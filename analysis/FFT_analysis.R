# Install packages if not already installed, then load them
if (!require(devtools)) install.packages('devtools')
if (!require(smisc)) devtools::install_github("stevenworthington/smisc")
smisc::ipak(c("utils", "tidyr", "dplyr", "ggplot2", "DescTools", "bspec",
              "pracma", "gridExtra", "data.table", "tables", "zoo", "parallel",
              "scales", "lazyeval", "stats", "gdata", "viridis", "gginnards"))

# Main function begins (encompasses other functions)
main_function <- function(ext_objects, display, iso_sides, sbtr, samp_per,
                          clumps, dep_var, pval, shuff, trends, smooth_method,
                          win_func, xaxisvals, duration, attn_filter,
                          catch_floor, side_bias, pre_range, post_range,
                          block_floor, miniblock_range, CTI_range){
  
  grouping_cnsts <- quos(participant, Trials_filtered_out, Acc_prefilter,       # Columns that are frequently used for grouping, variable means
                         Acc_postfilter, CatchAcc)                              # don't have to type them out every time we use them for grouping
  
  pcpts <- if(ext_objects == 2) 301:324 else 401:427 
  
  dep_var_abbr <- as.name(ifelse(dep_var == "Accuracy", "Acc", "RT"))                                                   # Converts character to name/symbol so we can refer to it as a column using tidy
  
  dem_df <- fread(file.path("data", "Demographics.csv")) %>%
    mutate_at(vars(SubjID), as.numeric)
  
  # For each participant, function reads in data, filters it, and transforms it to prepare for interpolation and FFT'ing
  pcpts_combine <- function(pcpt){
    fread(file.path("data", pcpt, paste0(pcpt, ".csv"))) %>%# Reads in participant data (the `select` part is because PsychoPy created an empty column at the end of the data frame for the first few participants, which meant those dataframes had different dimensions than subsequent data frames, which `do.call` doesn't like)
      filter(Trial > 0) %>%                                                     # Filters out practice trials
      mutate(CatchAcc = mean(ifelse(Opacity > 0, NA, Acc), na.rm = TRUE)) %>%   # Creates column indicating mean accuracy for catch trials
      group_by(CorrSide) %>%
      mutate(Side_Acc = mean(Acc, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(Side_Diff = max(Side_Acc) - min(Side_Acc)) %>%
      left_join(dem_df, by = c("participant" = "SubjID")) %>%
      filter(CatchAcc >= catch_floor,                                          # Filters out participants whose catch accuracy is below desired threshold
             grepl(ifelse(attn_filter, "fully alert" , ""), Q9),
             Opacity > 0,
             Side_Diff <= side_bias) %>%
      mutate(Acc_prefilter = mean(Acc, na.rm = TRUE),                           # Creates column indicating mean accuracy before we've filtered for `block_floor`, unlike `Acc_postfilter`
             CTI = RoundTo(RoundTo(lilsquareStartTime - flash_circleEndTime,
                                   1 / 60), samp_per)) %>% 
      filter(between(Acc_prefilter, pre_range[1], pre_range[2], incbounds = TRUE),                       # Filters out participants whose non-catch, pre-block-filtering accuracy is outside of desired range
             between(CTI, min(CTI_range), max(CTI_range))) %>%
      mutate(block = RoundTo(Trial, 54, ceiling) / 54,                          # Creates column indicating trial's block
             RT = ifelse(Acc == 1 & ButtonPressTime - lilsquareStartTime > .1,  #                indicating RT after target appeared on screen, only for correct trials with an RT < 100 ms
                         ButtonPressTime - lilsquareStartTime, NA),
             Stim_Sides = as.character(
               ifelse(CorrSide == FlashSide, "Valid", "Invalid")),              #                indicating whether cue was valid or invalid
             CorrSide = case_when(CorrSide == 1 ~ "Right",                      #                indicating which side the target appeared on
                                  CorrSide == -1 ~ "Left", TRUE ~ "Bottom"),
             Stim_Sides = case_when(iso_sides ~ paste(CorrSide, Stim_Sides, sep = "_"), # Overwrites `Stim_Sides` column if `sep_vis_fidels` parameter == `TRUE` to include which side of screen target was on,
                                    TRUE ~ Stim_Sides)) %>%                      # as well as whether it was valid with cue; if `iso_sides` parameter == `FALSE`, leaves `StimSides` unchanged
      group_by(block) %>%
      mutate(block_acc = mean(Acc)) %>%                                         # Creates column indicating block's mean accuracy
      ungroup() %>%
      mutate(rown = row_number(),
             miniblock = RoundTo(rown, 16, ceiling) / 16) %>%                   # Creates column indicating trial's mini-block (every 16 trials the opacity was readjusted)
      group_by(miniblock) %>%
      mutate(miniblock_avg = mean(Acc)) %>%
      ungroup() %>%
      mutate_at(vars(Acc, RT),
                list(~ifelse(block_acc <= block_floor | !between(miniblock_avg,  # Changes `Acc` and `RT` column values to NA if trial's block accuracy < `block_floor`
                                                                miniblock_range[1], miniblock_range[2]), NA, .))) %>%    # or miniblock was not in the desired range
      mutate(Trials_filtered_out = sum(is.na(Acc)) / n(),
             Acc_postfilter = mean(Acc, na.rm = TRUE)) %>%                      # Creates column indicating mean accuracy for non-catch trials; note this is after before we've filtered for `block_floor`, unlike `Acc_prefilter`
      filter(between(Acc_postfilter, post_range[1], post_range[2])) %>%                # Filters out participants whose non-catch, post-block-filtering accuracy is outside of desired range 
      group_by(CTI, Stim_Sides, !!!grouping_cnsts) %>%
      summarise_at(vars(Acc, RT), list(~mean(., na.rm = TRUE))) %>%              # Overwrites `Acc` and `RT` columns according to mean of each combination of `CTI` and `Stim_Sides`
      arrange(CTI) %>%
      group_by(Stim_Sides) %>%
      mutate_at(vars(Acc, RT), list(~na.approx(., na.rm = FALSE, rule = 2))) %>%
      mutate_at(vars(Acc, RT), list(~rollapply(., clumps + 1, mean, partial = TRUE)))
  }
  
  cmbd <- do.call(rbind, lapply(pcpts, pcpts_combine)) %>%             # Calls `pcpts_combine` function for argumenet `pcpts`; then combines each participant's dataframe into one
    arrange(Acc_prefilter, participant, Stim_Sides, CTI)
  CTIs <- unique(cmbd$CTI)
  if (win_func == "Tukey"){ win <- tukeywindow(length(CTIs), .5)} else { # Creates window, which if `tukey` will add the parameter `r` == `.5` —so 'only' half the data length will be non-flat
    win <- match.fun(paste0(tolower(win_func), "window"))(length(CTIs))}
  locations <- unique(cmbd$Stim_Sides)                                   # Creates vector of column names representing sides locations of target in reference to cue (and also potentially side of screen)
  pcpts <- unique(cmbd$participant)                                       # Creates vector of remaining participant numbers after `pcpts_combine` filtering
  cmbd_w <- cmbd %>%
    pivot_wide(CTI:CatchAcc, Stim_Sides, !!dep_var_abbr) %>%
    arrange(Acc_prefilter, participant, CTI) %>%
    group_by(participant) %>%
    mutate_at(vars(locations), list(~na.approx(., na.rm = FALSE, rule = 2)))
  
  # Analyzes Invalid - Valid instead of them separetely, if sbtr == `TRUE`
  if (sbtr){
    s <- tail(1:ncol(cmbd_w), length(locations)) [c(TRUE, FALSE)]
    cmbd_w[paste0(names(cmbd_w[s]), "_minus_",
                  names(cmbd_w)[s + 1])] <- cmbd_w[s] - cmbd_w[s + 1]
    locations = tail(colnames(cmbd_w), length(locations) / 2)
  }
  
  # Determines confidence intervals
  conf_int <- function(x, ...){ x %>%
      group_by_(.dots = lazy_dots(...)) %>%
      summarise_at(vars(locations), list(~qnorm(.975) * std_err(.)))
  }
  
  # Transforms from Time to Frequency Domain
  amplitude <- function(x, z){
    pre_pad <- length(pcpts) * (length(unique(cmbd_w$CTI))) * z
    x %>%
      group_by(participant) %>%
      mutate_at(vars(locations), list(~ case_when("Detrending" %in% trends ~ . - polyval(polyfit(CTI, ., 2), CTI), TRUE ~ .))) %>%
      mutate_at(vars(locations), list(~ case_when("Demeaning" %in% trends ~ . - mean(.), TRUE ~ .))) %>%
      mutate_at(vars(locations), list(~ . * win / Norm(win))) %>%
      ungroup() %>%
      add_row(participant = rep(pcpts, times = z * (ceiling((duration - diff(range(cmbd_w$CTI)))/samp_per)))) %>%
      head(-length(pcpts) * z) %>%
      mutate_at(vars(locations), list(~coalesce(., 0))) %>%
      mutate(samp_shuff = ifelse(row_number() <= pre_pad, 
                                 RoundTo(row_number(), pre_pad / z,
                                         ceiling) / (pre_pad / z), 
                                 RoundTo(row_number() - pre_pad, (n() - pre_pad) / (z),
                                         ceiling) / ((n() - pre_pad) /  z))) %>%
      group_by(participant, samp_shuff) %>%
      mutate_at(vars(locations), list(~Mod(sqrt(2/n()) * fft(.))^2)) %>%
      mutate(Hz = (row_number() - 1) / (n() * samp_per)) %>%
      ungroup() %>%
      select(-CTI)
  }
  
  amps <- amplitude(cmbd_w, 1)
  fft_x <- round(1 / (length(unique(amps$Hz)) * samp_per),1)
  xaxis_r <- RoundTo(xaxisvals, fft_x)
  
  
  # Set Up Graphing
  
  t_srs_g <- function(x){
    cmbd_w %>%
      gather(Location, !!dep_var_abbr, -c(CTI, !!!grouping_cnsts)) %>%
      right_join(gather(conf_int(cmbd_w, CTI), Location, Conf_Int, -CTI),
                 by = c("CTI", "Location")) %>%
      group_by(CTI, Location, Conf_Int, !!!head(grouping_cnsts, x)) %>%
      summarise(!!dep_var_abbr := mean(!!dep_var_abbr)) %>%                               # Keeps either `RT` or `Acc` column—depending on whether `dep_var_abbr` parameter == `RT` or `Acc`
      ggplot(aes(CTI, !!dep_var_abbr, group = Location, color = Location,
                 fill = Location, ymin = !!dep_var_abbr - Conf_Int,
                 ymax = !!dep_var_abbr + Conf_Int)) + 
      labs(title = paste0(dep_var_abbr, " by Cue-Target Interval ", ext_objects, "-object Task"),
           x = "Cue-Target Interval (ms)")
  }
  fft_g <- function(x) { x +
      labs(title = paste0("FFT of Target ", dep_var, ", ", ext_objects, "-object task"),
           col = "Target Location", y = "Spectral Power") +
      theme(panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank())
  }
  
  viridis_cols <- .7 + RoundTo(.0001 * RoundTo(length(locations),
                                               4, floor), .2, ceiling)
  graph <- function(y, x) {
    y(x) +
      theme_bw() +
      scale_color_viridis_d(option = "C",
                            end = viridis_cols,
                            labels = sapply(locations, function(x) gsub("_", " ", x),
                                            USE.NAMES = FALSE, simplify = TRUE)) +
      scale_fill_viridis_d(option = "C",
                           end = viridis_cols) +
      guides(colour = guide_legend(reverse = TRUE), fill = FALSE)
  }
  
  idvl_g <- function(y, x){
    graph(y, x) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)) +
      facet_wrap(~factor(participant, levels = pcpts),
                 ncol = round(length(pcpts) / 3), scales = 'free_x')                                   # `free` means the y_axis isn't fixed from participant to participant
  }
  cmbd_g <- function(y, x) {
    graph(y, x) + geom_line(size = 1.5) +
      theme(legend.key.size = unit(.55, "in")) +
      labs(subtitle = paste( "Data from", as.character(length(pcpts)),
                            "participants"))
  }
  
  sv_cmbd_g <- function(x) { x %>%
      ggsave(filename = file.path("plots", paste0(display, ".pdf")),
             width = xaxisvals)}
  
  # Conditionally returns data frame of prelim participant or post-FFT data, exits function
  if (display == "prelim_table") {
      write.csv(cmbd, file.path("plots", "Prelim_table.csv"))
  } else if (display == "fft_table") {
       write.csv(amps, file.path("plots", "FFT_table.csv"))
  } else if (display == "Time-Series + FFT by Individual"){
    
    # display Individuals' Plots ------------------------------------------------
    
    # Produces left half of final graph
    ts_facets <- idvl_g(t_srs_g, 1) +
      geom_line(alpha = I(2 / 10), color = "grey", show.legend = FALSE) +       # Graphs unsmoothed data in light gray
      stat_smooth(method = tolower(smooth_method), span = 0.2, se = FALSE,               # Smoothes data depending on `smooth_method` parameter
                  size = .5, show.legend = FALSE)
    
    # Produces label for each right-side graph
    plot_label <- amps %>%
      group_by(!!!grouping_cnsts) %>%
      summarise() %>%
      ungroup() %>%
      drop_na() %>%
      mutate_at(vars(!!!tail(grouping_cnsts, -1)),
                list(~paste(quo_name(quo(.)), "=", percent(.)))) %>%
      unite(lab, !!!tail(grouping_cnsts, -1), sep = "\n", remove = FALSE)

    # Produces right half of final graph
    fft_facets <- idvl_g(fft_g, amps %>%
                           group_by(participant, Hz) %>%
                           summarise_all(mean) %>%
                           gather(Flash_and_or_field, Power, -Hz, -samp_shuff, -c(!!!grouping_cnsts)) %>%
                           ggplot(aes(Hz, Power, color = Flash_and_or_field))) +
      geom_line() +
      scale_x_continuous(name = "Frequency (Hz)", limits = c(0, xaxisvals),
                         breaks = seq(0, xaxisvals, ifelse(xaxisvals > 10 | xaxisvals != xaxis_r, max(Closest(xaxis_r/ seq(fft_x, xaxis_r, fft_x), 5) /xaxis_r)^-1, fft_x))) +
      labs(caption = paste("Data from", as.character(length(pcpts)),
                           "participants")) +
      geom_text(data = as.data.frame(plot_label), inherit.aes = FALSE, size = 1.2,# Sets location for label overlayed onto graph
                aes(label = lab, x = Inf, y = Inf), vjust = 1.15, hjust = 1.05)
    
    side_by_side <- arrangeGrob(ts_facets, fft_facets, ncol = 2)                # Combines time series and FFT graphs into one plot
    ggsave(file.path("plots", "Indvls_Plots.pdf"), width = 25, side_by_side)
    
  } else if (display == "Time-Series Across Participants") { # Graph combined Time Series -----------
    
    (move_layers(cmbd_g(t_srs_g, 0) +
                   theme(panel.grid = element_blank()) +
                   geom_ribbon(alpha = 0.15, aes(color = NULL)), "GeomRibbon", position = "bottom")) %>%
      sv_cmbd_g
    
  } else { # Graph combined FFT ------------------------------------------------

    set.seed(123)
    fft_x <- 1 / (length(unique(amps$Hz)) * samp_per)
    
    # Produces 'shuff' # of null hypothesis permutations
    shuffle <- function(x){
      cmbd_w %>%
        group_by(participant) %>%
        sample_n(length(CTIs), weight = CTI) %>%
        mutate_at(vars(CTI), list(~seq(min(CTIs), max(CTIs), samp_per))) %>%
        mutate(samp_shuff = x)
    }
    
    amps_shuff <- do.call(rbind, lapply(1:shuff, shuffle)) %>%
      amplitude(shuff) %>%
      group_by(Hz, samp_shuff) %>%
      summarise_at(vars(locations), mean) %>%
      group_by(Hz) %>%
      summarise_at(vars(locations), list(~quantile(., probs = 1 - pval))) %>%
      combine(amps %>% group_by(Hz) %>% summarise_at(vars(locations), mean),
              names = (c("Significance Cutoff", "Observed Data"))) %>%
      gather(Location, Power, -c(Hz, source)) %>%
      right_join(gather(conf_int(amps, Hz), Location, Conf_Int, -Hz),
                 by = c("Hz", "Location"))
    (move_layers(cmbd_g(fft_g, ggplot(amps_shuff, aes(Hz, Power, col = Location, linetype = source, 
                                                      ymin = Power - Conf_Int, ymax = Power + Conf_Int, fill = Location))) +
                   scale_linetype_manual(values = c("solid", "dashed")) +
                   scale_x_continuous(name = "Frequency (Hz)", limits = c(0, xaxisvals),
                                      breaks = seq(0, xaxisvals, ifelse(fft_x > .5, round(fft_x, 2), 1))) +
                   labs(linetype = "",
                        caption = paste("Significance threshold at p < ",
                                         as.character(pval))) +
                   geom_ribbon(data = filter(amps_shuff, source == "Observed Data"), alpha = 0.15, aes(color = NULL)) +
                   geom_point(size = 3, data = amps_shuff %>% spread(source, Power) %>%
                                filter(`Observed Data` > `Significance Cutoff`) %>%
                                select(-c(`Significance Cutoff`, Conf_Int)) %>%
                                gather(source, Power, -Hz, -Location), 
                              aes(ymin = NULL, ymax = NULL)), "GeomRibbon", position = "bottom")) %>%
      sv_cmbd_g
    
  }
}



# Sets inputs for the 'main_function' function
main_function(ext_objects = 2,                                                  # `2` corresponds to the two-object task, `3` to the three-object task
              display = "FFT Across Participants",                              # Either `FFT Across Participants`, `Time-Series Across Participants`, `Time-Series + FFT by Individual`, `prelim_table` (lightly analyzed data), and `fft_table` (semi-ready for graphing data)
              iso_sides = FALSE,                                                 # Either `FALSE` or `TRUE`, which groups by not only valid and invalid but also by the side of the screen for each trial (i.e. going from `Valid` and `Invalid` to `Right Valid`, `Left Valid`, `Right Invalid`, `Left Invalid`)
              sbtr = FALSE,                                                     # Either `FALSE` or `TRUE`, which subtracts the dependent variable values at each CTI (valid - invalid) before performing analyses rather than analyzing valid and invalid trials independently
              samp_per = 1 / 60,                                                # Spacing between CTI intevals (in seconds); the data was originally sampled at 1 / 60, but one could re-sample at a different rate, which would just clump neighboring CTI's together (whereas the 'clump' variable groups neighbors but doesn't combine them, keeping the same total number of bins)
              clumps = 0,                                                       # Number of points to average at each CTI; `1` means this function does nothing, `3` means each CTI is the average of that CTI and its neighboring CTI's on each sides, etc...
              dep_var = "Accuracy",                                             # Either `Accuracy`or `Response Time`
              pval = .05,                                                       # The p-value to use for drawing the significance cutoff on the graphs
              shuff = 50,                                                       # The number of surrogate shuffles to use to determine the null hypothesis; NOTE: increasing this number slows down the run time
              trends = c("Detrending", "Demeaning"),                            # Either `Detrending`, `Demeaning`, both, or an empty vector
              smooth_method = "Loess",                                          # Either `Loess`, `LM`, `GLM`, or `GAM` smoothing methods for the time-series data of individuals
              win_func = "Tukey",                                               # Choose between the following types of windowing functions: `Tukey`, `Square`, `Hann`, `Welch`, `Triangle`, `Hamming`, `Cosine`, or `Kaiser`
              xaxisvals = 25,                                                   # Greatest x-axis value included in graph
              duration = 1,                                                     # Duration (Seconds) Analyzed Including Padding
              attn_filter = FALSE,                                              # Either `FALSE` or `TRUE`, which filters out any participants who indicated in the post-task questionnaire that they either dozed off at one point or were not fully focused on at least two of the eight blocks (the other options were that they were fully alert on all blocks or fully alert on all but one block)
              catch_floor = .85,                                                # Filter participants who perform below this hit rate on catch trials, which featured no target
              side_bias = .3,                                                   # Filter out participants whose hit rate at one visual field - another visual field is > `side_bias`
              pre_range = c(.45, .75),                                          # Filter out participants whose unfiltered/scrutinized data is outside of the selected range
              post_range = c(.45, .75),                                                 # Filter out participants whose filtered data is outside of the selected range
              block_floor = .65,                                                # Interpolate over trials if the average hit rate in that block, every 48 trials, was below this value
              miniblock_range = c(.45, .85),                                           # Interpolate over trials if the average hit rate in that mini-block, every 16 trials which is how often the task difficulty was adjusted to titrate to 65%, is outside this range
              CTI_range = c(.5, 1.29)                                          # Filters out trials outside this range of CTI bins
)
