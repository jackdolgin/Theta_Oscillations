# Install packages if not already installed, then load them
if (!require(devtools)) install.packages('devtools')
if (!require(smisc)) devtools::install_github("stevenworthington/smisc")
smisc::ipak(c("utils", "tidyr", "dplyr", "ggplot2", "DescTools", "bspec",
              "pracma", "gridExtra", "data.table", "tables", "zoo", "parallel",
              "scales", "purrr", "lazyeval", "stats", "gdata", "viridis",
              "gginnards"))

# Main function begins (encompasses other functions)
data_graphing <- function(catch_cutoff, block_floor, mini_block_floor,
                          mini_block_ceil, side_bias, attn_filter, catch,
                          pre_floor, pre_ceil, post_floor, post_ceil,
                          iso_sides, smooth_method, win_func, dep_var, sbtr,
                          samp_per, zeropads, latestart, earlyend, ext_objects,
                          output, clumps, xaxisvals, shuff, pval){
  
  grouping_cnsts <- quos(participant, Trials_filtered_out, Acc_prefilter,       # Columns that are frequently used for grouping, variable means
                         Acc_postfilter, CatchAcc)                          # don't have to type them out every time we use them for grouping
  
  pcpts <- if(ext_objects == 2) 301:324 else 401:427 
  
  dep_var <- as.name(dep_var)                                                   # Converts character to name/symbol so we can refer to it as a column using tidy
  
  dem_df <- fread(file.path("data", "Demographics.csv")) %>%
    mutate_at(vars(SubjID), as.numeric)
  
  # For each participant, function reads in data, filters it, and transforms it to prepare for interpolation and FFT'ing
  pcpts_combine <- function(pcpt){
    fread(file.path("data", pcpt, paste0(pcpt, ".csv")), select = c(1:23)) %>%# Reads in participant data (the 'select' part is because PsychoPy created an empty column at the end of the data frame for the first few participants, which meant those dataframes had different dimensions than subsequent data frames, which 'do.call' doesn't like)
      filter(Trial > 0) %>%                                                     # Filters out practice trials
      mutate(CatchAcc = mean(ifelse(Opacity > 0, NA, Acc), na.rm = TRUE)) %>%   # Creates column indicating mean accuracy for catch trials
      group_by(CorrSide) %>%
      mutate(Side_Acc = mean(Acc, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(Side_Diff = max(Side_Acc) - min(Side_Acc)) %>%
      left_join(dem_df, by = c("participant" = "SubjID")) %>%
      filter(CatchAcc >= catch_cutoff,                                          # Filters out participants whose catch accuracy is below desired threshold
             grepl(ifelse(attn_filter == "On", "fully alert" , ""), Q9),
             Opacity > catch,                                                   #             catch trials if 'catch' parameter is assigned to 0; if it's assigned to -1, this line does nothing)
             Side_Diff <= side_bias) %>%
      mutate(Acc_prefilter = mean(Acc, na.rm = TRUE),                           # Creates column indicating mean accuracy before we've filtered for 'block_floor', unlike 'Acc_postfilter'
             CTI = RoundTo(RoundTo(lilsquareStartTime - flash_circleEndTime,
                                   1 / 60), samp_per)) %>% 
      filter(between(Acc_prefilter, pre_floor, pre_ceil, incbounds = TRUE),                       # Filters out participants whose non-catch, pre-block-filtering accuracy is outside of desired range
             between(CTI, min(CTI) + latestart, max(CTI) - earlyend)) %>%
      mutate(block = RoundTo(Trial, 54, ceiling) / 54,                          # Creates column indicating trial's block
             RT = ifelse(Acc == 1 & ButtonPressTime - lilsquareStartTime > .1,  #                indicating RT after target appeared on screen, only for correct trials with an RT < 100 ms
                         ButtonPressTime - lilsquareStartTime, NA),
             Stim_Sides = as.character(
               ifelse(CorrSide == FlashSide, "Valid", "Invalid")),              #                indicating whether cue was valid or invalid
             CorrSide = case_when(CorrSide == 1 ~ "Right",                      #                indicating which side the target appeared on
                                  CorrSide == -1 ~ "Left", TRUE ~ "Bottom"),
             Stim_Sides = case_when(iso_sides == "No" ~ Stim_Sides,        # Overwrites 'Stim_Sides' column if 'sep_vis_fidels' parameter == 'Yes' to include which side of screen target was on,
                                    TRUE ~ paste(CorrSide, Stim_Sides, sep = "_"))) %>%             # as well as whether it was valid with cue; if 'iso_sides' parameter == 'No', leaves 'StimSides' unchanged
      group_by(block) %>%
      mutate(block_acc = mean(Acc)) %>%                                         # Creates column indicating block's mean accuracy
      group_by(Opacity > 0) %>%
      mutate(rown = row_number(),
             miniblock = ifelse(Opacity > 0,                                    # Creates column indicating trial's mini-block (every 16 trials the opacity was readjusted)
                                RoundTo(rown, 16, ceiling) / 16, NA)) %>%
      group_by(miniblock) %>%
      mutate(miniblock_avg = mean(Acc)) %>%
      ungroup() %>%
      mutate_at(vars(Acc, RT),
                funs(ifelse(block_acc <= block_floor | !between(miniblock_avg,  # Changes 'Acc' and 'RT' column values to NA if trial's block accuracy < 'block_floor'
                                                                mini_block_floor, mini_block_ceil), NA, .))) %>%    # or mini_block was not in the desired range
      mutate(Trials_filtered_out = sum(is.na(Acc)) / n(),
             Acc_postfilter = mean(Acc, na.rm = TRUE)) %>%                      # Creates column indicating mean accuracy for non-catch trials; note this is after before we've filtered for block_floor, unlike 'Acc_prefilter'
      filter(between(Acc_postfilter, post_floor, post_ceil)) %>%                # Filters out participants whose non-catch, post-block-filtering accuracy is outside of desired range 
      group_by(CTI, Stim_Sides, !!!grouping_cnsts) %>%
      summarise_at(vars(Acc, RT), funs(mean(., na.rm = TRUE))) %>%              # Overwrites 'Acc' and 'RT' columns according to mean of each combination of 'CTI' and 'Stim_Sides'
      arrange(CTI) %>%
      group_by(Stim_Sides) %>%
      mutate_at(vars(Acc, RT), funs(na.approx(., na.rm = FALSE, rule = 2))) %>%
      mutate_at(vars(Acc, RT), funs(rollapply(., clumps, mean, partial = TRUE)))
  }
  
  cmbd <- do.call(rbind, mclapply(pcpts, pcpts_combine)) %>%             # Calls 'pcpts_combine' function for argumenet 'pcpts'; then combines each participant's dataframe into one
    arrange(Acc_prefilter, participant, Stim_Sides, CTI)
  if (win_func == "tukey"){ win <- tukeywindow(length(unique(cmbd$CTI)), .5)} else { # Creates window, which if 'tukey' will add the parameter 'r == .5'—so 'only' half the data length will be non-flat
    win <- match.fun(paste0(win_func, "window"))(length(unique(cmbd$CTI)))}
  locations <- unique(cmbd$Stim_Sides)                                   # Creates vector of column names representing sides locations of target in reference to cue (and also potentially side of screen)
  pcpts <- unique(cmbd$participant)                                       # Creates vector of remaining participant numbers after 'pcpts_combine' filtering
  cmbd_w <- cmbd %>%
    pivot_wide(CTI:CatchAcc, Stim_Sides, !!dep_var) %>%
    arrange(Acc_prefilter, participant, CTI) %>%
    group_by(participant) %>%
    mutate_at(vars(locations), funs(na.approx(., na.rm = FALSE, rule = 2)))
  
  # Analyzes Invalid - Valid instead of them separetely, if sbtr == "Yes"
  if (sbtr == "Yes"){
    s <- tail(1:ncol(cmbd_w), length(locations)) [c(TRUE, FALSE)]
    cmbd_w[paste0(names(cmbd_w[s]), "_minus_",
                  names(cmbd_w)[s + 1])] <- cmbd_w[s] - cmbd_w[s + 1]
    locations = tail(colnames(cmbd_w), length(locations) / 2)
  }
  
  # Determines confidence intervals
  conf_int <- function(x, ...){ x %>%
      group_by_(.dots = lazy_dots(...)) %>%
      summarise_at(vars(locations), funs(qnorm(.975) * std_err(.)))
  }
  
  # Transforms from Time to Frequency Domain
  amplitude <- function(x, z){
    pre_pad <- length(pcpts) * (length(unique(cmbd_w$CTI))) * z
    x %>%
      group_by(participant) %>%
      mutate_at(vars(locations), funs(detrend(.) * win)) %>%
      ungroup() %>%
      add_row(participant = rep(pcpts, times = z * (zeropads + 1))) %>%
      head(-length(pcpts) * z) %>%
      mutate_at(vars(locations), funs(coalesce(., 0))) %>%
      mutate(samp_shuff = ifelse(row_number() <= pre_pad, 
                                RoundTo(row_number(), pre_pad / z,
                                        ceiling) / (pre_pad / z), 
                                RoundTo(row_number() - pre_pad, (n() - pre_pad) / (z),
                                        ceiling) / ((n() - pre_pad) /  z))) %>%
      group_by(participant, samp_shuff) %>%
      mutate_at(vars(locations), funs(Mod(fft(.)))) %>%          # Detrends, multiplies by window, applies FFT, and then takes magnitude
      mutate(Hz = (row_number() - 1) / (n() * samp_per)) %>%
      ungroup() %>%
      select(-CTI)
  }
  
  amps <- amplitude(cmbd_w, 1)
  
  
  # Set Up Graphing
  
  t_srs_g <- function(x){
    cmbd_w %>%
      gather(Location, !!dep_var, -c(CTI, !!!grouping_cnsts)) %>%
      right_join(gather(conf_int(cmbd_w, CTI), Location, Conf_Int, -CTI),
                 by = c("CTI", "Location")) %>%
      group_by(CTI, Location, Conf_Int, !!!head(grouping_cnsts, x)) %>%
      summarise(!!dep_var := mean(!!dep_var)) %>%                               # Keeps either RT or Acc column—depending on whether 'dep_var' parameter == 'RT' or 'Acc'
      ggplot(aes(CTI, !!dep_var, group = Location, color = Location,
                 fill = Location, ymin = !!dep_var - Conf_Int,
                 ymax = !!dep_var + Conf_Int)) + 
      labs(title = paste(dep_var, "by Cue-Target Interval", ", ", ext_objects, "-object task"),
           x = "Cue-Target Interval (ms)")
  }
  fft_x <- 1 / (length(unique(amps$Hz)) * samp_per)
  fft_g <- function(x) { x +
      labs(title = paste0("FFT of Target ", dep_var, ", ", ext_objects, "-object task"),
           col = "Target Location") +
      theme(panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank()) +
      scale_x_continuous(name = "Frequency (Hz)", limits = c(0, xaxisvals), 
                         breaks = seq(0, xaxisvals, ifelse(fft_x > .5, round(fft_x, 2), 1)))
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
                 ncol = 4, scales = 'free_x')                                   # 'free' means the y_axis isn't fixed from participant to participant
  }
  cmbd_g <- function(y, x) {
    graph(y, x) + geom_line(size = 1.5) +
      theme(legend.key.size = unit(.55, "in")) +
      labs(subtitle = paste( "Data from", as.character(length(pcpts)),
                            "participants"))
  }
  
  sv_cmbd_g <- function(x) { x %>%
      ggsave(filename = file.path("analysis", paste0(output, ".pdf")),
             width = xaxisvals)}
  
  # Conditionally returns data frame of prelim participant or post-FFT data, exits function
  if (output == "prelim_table") {
      write.csv(cmbd, file.path("analysis", "Prelim_table.csv"))
  } else if (output == "fft_table") {
       write.csv(amps, file.path("analysis", "FFT_table.csv"))
  } else if (output == "individuals"){
    
    # Output Individuals' Plots ------------------------------------------------
    
    # Produces left half of final graph
    ts_facets <- idvl_g(t_srs_g, 1) +
      geom_line(alpha = I(2 / 10), color = "grey", show.legend = FALSE) +       # Graphs unsmoothed data in light gray
      stat_smooth(method = smooth_method, span = 0.2, se = FALSE,               # Smoothes data depending on 'smooth_method' parameter
                  size = .5, show.legend = FALSE)
    
    # Produces label for each right-side graph
    plot_label <- amps %>%
      group_by(!!!grouping_cnsts) %>%
      summarise() %>%
      ungroup() %>%
      mutate_at(vars(!!!tail(grouping_cnsts, -1)),
                funs(paste(quo_name(quo(.)), "=", percent(.)))) %>%
      unite(lab, !!!tail(grouping_cnsts, -1), sep = "\n", remove = FALSE)
    
    # Produces right half of final graph
    fft_facets <- idvl_g(fft_g, amps %>%
                           group_by(participant, Hz) %>%
                           summarise_all(mean) %>%
                           gather(Flash_and_or_field, Magnitude, -Hz, -c(!!!grouping_cnsts)) %>%
                           ggplot(aes(Hz, Magnitude, color = Flash_and_or_field))) +
      geom_line() +
      labs(caption = paste("Data from", as.character(length(pcpts)),
                           "participants")) +
      geom_text(data = as.data.frame(plot_label), inherit.aes = FALSE, size = 2,# Sets location for label overlayed onto graph
                aes(label = lab, x = Inf, y = Inf), vjust = 1.15, hjust = 1.05)
    
    side_by_side <- arrangeGrob(ts_facets, fft_facets, ncol = 2)                # Combines time series and FFT graphs into one plot
    ggsave(file.path("analysis", "Indvls_Plots.pdf"), width = 25, side_by_side)
    
  } else if (output == "combined_ts") { # Graph combined Time Series -----------
    
    (move_layers(cmbd_g(t_srs_g, 0) +
                   theme(panel.grid = element_blank()) +
                   geom_ribbon(alpha = 0.15, aes(color = NULL)), "GeomRibbon", position = "bottom")) %>%
      sv_cmbd_g
    
  } else { # Graph combined FFT ------------------------------------------------
    
    # Produces 'shuff' # of null hypothesis permutations
    shuffle <- function(.data, n, perm_cols){
      cols_ids <- match(perm_cols, colnames(.data))
      ids <- seq_len(nrow(.data))
      n_ids <- rerun(n, sample(ids))
      
      map_dfr(n_ids, function(x){
        .data[ids, cols_ids] <- .data[x, cols_ids]
        .data
      })
    }
    
    set.seed(123)
    
    # Produces and save graph
    amps_shuff <- shuffle(.data = cmbd_w, n = shuff,
                          perm_cols = c("participant", locations)) %>%
      amplitude(shuff) %>%
      group_by(Hz, samp_shuff) %>%
      summarise_at(vars(locations), mean) %>%
      group_by(Hz) %>%
      summarise_at(vars(locations), funs(quantile(., probs = 1 - pval))) %>%
      combine(amps %>% group_by(Hz) %>% summarise_at(vars(locations), mean),
              names = (c("Significance Cutoff", "Observed Data"))) %>%
      gather(Location, Magnitude, -c(Hz, source)) %>%
      right_join(gather(conf_int(amps, Hz), Location, Conf_Int, -Hz),
                 by = c("Hz", "Location"))
    (move_layers(cmbd_g(fft_g, ggplot(amps_shuff, aes(Hz, Magnitude, col = Location, linetype = source, 
                                                      ymin = Magnitude - Conf_Int, ymax = Magnitude + Conf_Int, fill = Location))) +
                   scale_linetype_manual(values = c("solid", "dashed")) +
                   labs(linetype = "",
                        caption = paste0(zeropads, " zero pads added","\nSignificance threshold at p < ",
                                         as.character(pval))) +
                   geom_ribbon(data = filter(amps_shuff, source == "Observed Data"), alpha = 0.15, aes(color = NULL)) +
                   geom_point(size = 3, data = amps_shuff %>% spread(source, Magnitude) %>%
                                filter(`Observed Data` > `Significance Cutoff`) %>%
                                select(-c(`Significance Cutoff`, Conf_Int)) %>%
                                gather(source, Magnitude, -Hz, -Location), 
                              aes(ymin = NULL, ymax = NULL)), "GeomRibbon", position = "bottom")) %>%
      sv_cmbd_g
    
  }
}



# Sets inputs for the 'data_graphing' function
data_graphing(catch_cutoff = .85,                                               # Filters participants by their accuracy on catch trials (no target on the screen), which is used as a manipualtion check 
              block_floor = .65,                                                # Interpolates over trials inside a block (every 48 trials) averaging less than this hit rate
              mini_block_floor = .45,                                           # Interpolates over trials inside a mini-block (every 16 trials, after which stimulus difficulty re-tritates to 65%) averaging outside this range of hit rates
              mini_block_ceil = .85,                                            # Line above is the floor of the filter, this line is the ceiling
              side_bias = .3,                                                   # Filter out participants whose hit rate at one visual field (i.e. right or left side) was > `side_bias` better than another visual field
              attn_filter = "Off",                                              # Either `Off` or `On`, which filters out any participants who indicated in the post-task questionnaire that they either dozed off at one point or were not fully focused on at least two of the eight blocks (the other options were that they were fully alert on all blocks or fully alert on all but one block)
              catch = 0,                                                        # Either `-1` or `0` to include or exclude catch trials from the analyses (besides the participant meeting the catch trial cutoff); for example `-1` would mean analyzing hit rates or response times across all trials, mixing catch and non-catch trials together (response times will be messier since a correct catch trial response is not responding for a full second, therefore an inverse relationship between response time and accuracy)
              pre_floor = .45,                                                  # Filter out participants whose unfiltered/scrutinized data is outside of the selected range
              pre_ceil = .85,                                                   # Line above is floor, this line the ceiling of the range
              post_floor = .35,                                                 # Filter out participants whose filtered data is outside of the selected range
              post_ceil = .75,                                                  # Line above is the floor of the filter, this line is the ceiling
              iso_sides = "No",                                                 # Either `No` or `Yes`, which groups by not only valid and invalid but also by the side of the screen for each trial (i.e. going from `Valid` and `Invalid` to `Right Valid`, `Left Valid`, `Right Invalid`, `Left Invalid`)
              smooth_method = "loess",                                          # Either `loess`, `lm`, `glm`, or `gam` smoothing methods for the time-series data of individuals
              win_func = "tukey",                                               # Choose between the following types of windowing functions: `tukey`, `square`, `hann`, `welch`, `triangle`, `hamming`, `cosine`, or `kaiser`
              dep_var = "Acc",                                                  # Either `Acc` (accuracy) or `RT` (response time)
              sbtr = "No",                                                      # Either `No` or `Yes`, which subtracts the dependent variable values at each CTI (valid - invalid) before performing analyses rather than analyzing valid and invalid trials independently)
              samp_per = 1 / 60,                                                # Spacing between CTI intevals (in seconds); the data was originally sampled at 1 / 60, but one could re-sample at a different rate, which would just clump neighboring CTI's together (whereas the 'clump' variable groups neighbors but doesn't combine them, keeping the same total number of bins)
              zeropads = 0,                                                     # The number of zero pads to use; can be set from 0 upwards
              latestart = 0, earlyend = 0,                                      # Filters out trials within the first `latestart` seconds of CTI bins or the last `earlyend` seconds of CTI bins
              ext_objects = 2,                                                  # `2` corresponds to the two-object task, `3` to the three-object task
              output = "combined_fft",                                          # Either `combined_fft`, `combined_ts`, `individuals`, `prelim_table` (lightly analyzed data), and `fft_table` (semi-ready for graphing data)
              clumps = 1,                                                       # Number of points to average at each CTI; `1` means this function does nothing, `3` means each CTI is the average of that CTI and its neighboring CTI's on each sides, etc...
              xaxisvals = 10,                                                   # Greatest x-axis value included in graph
              shuff = 50,                                                       # The number of surrogate shuffles to use to determine the null hypothesis; NOTE: increasing this number slows down the run time
              pval = .05)                                                       # The p-value to use for drawing the significance cutoff on the graphs
