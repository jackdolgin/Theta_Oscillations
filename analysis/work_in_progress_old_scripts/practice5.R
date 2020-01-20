# Install packages if not already installed, then load them
if (!require(devtools)) install.packages("pacman")
pacman::p_load(
  utils, tidyr, dplyr, ggplot2, DescTools, bspec, pracma, gridExtra, data.table,
  tables, zoo, scales, stats, viridis, gginnards, purrr, stingr, rlang, here)

# Main function begins (encompasses other functions)
main_function <- function(display, dset, wm_exp, iso_sides, sbtr, samp_per,
                          clumps, dep_var, α, shuff, mult_correcs, trends,
                          smooth_method, win_func, xmax, duration, attn_filter,
                          catch_floor, side_bias, wm_floor, invalid_floor,
                          pre_range, post_range, filtered_cap, block_range,
                          blocks_desired, miniblock_range, CTI_range){
  
  grouping_cnsts <- quos(                                                       # Columns that are frequently used for grouping, variable means...
    participant, Trials_filtered_out, Acc_prefilter, Acc_postfilter, CatchAcc)  # ... don't have to type them out every time we use them for grouping
  
  dvcol <- (ifelse(dep_var == "Accuracy", "Acc", "RT")) %>% as.name             # Converts character to name/symbol so we can refer to it as a column using tidycom
  
  dv_mutate <- function(x, s){
    x %>%
      mutate_at(vars(!!dvcol), list(~!!s))
  }
  
  if (batch == "3a"){
    blocksize <- 79
  }
  
  cmbd <-
    dir(path = here("data", dset),
      pattern = "*.csv",
      full.names = TRUE,
      recursive = TRUE) %>%
    map_df(fread) %>%
    filter(Trial > 0) %>%                                                     # Prunes practice trials
    mutate(
      Stim_Sides =
        ifelse(CorrSide == FlashSide, "Valid", "Invalid") %>% as.character,     # Creates column indicating whether cue was valid or invalid
      CTI =
        (lilsquareStartTime - flash_circleEndTime) %>%
        RoundTo(1 / 60) %>%
        RoundTo(samp_per),
      block =
        RoundTo(Trial, blocksize, ceiling) / blocksize,                         # Creates column indicating trial's block
      RT =
        ifelse(Acc == 1 & ButtonPressTime - lilsquareStartTime > .1,       #                           RT after target appeared on screen, only for correct trials with an RT > 100 ms
               ButtonPressTime - lilsquareStartTime,
               NA)) %>%
    dv_mutate(expr(ifelse(block %in% blocks_desired, .,  NA))) %>%              # Converts trials' Acc and RT to NA if they are not in a desired block
    mutate(
      wm_Acc =
        ifelse(session == "4",
               mean(Acc_wmarith, na.rm = TRUE),
               1),
      CatchAcc =
        ifelse(Opacity != 0, NA, Acc) %>% mean(na.rm = TRUE)) %>%  #whats up with taking this mean?
    filter(
      CatchAcc >= catch_floor,                                           # Prunes participants whose catch accuracy is below desired threshold
      Opacity != 0) %>%                                                  #        catch trials
    filter(mean(Acc[Stim_Sides == "Invalid"], na.rm = TRUE) >= 
             invalid_floor) %>%
    mutate(
      CorrSide =
        case_when(CorrSide == 1 ~ "Right",
                  CorrSide == -1 ~ "Left",
                  TRUE ~ "Bottom"),
      Stim_Sides =
        case_when(iso_sides ~ paste(CorrSide, Stim_Sides,  sep = "_"),
                  TRUE ~ Stim_Sides)) %>%                          # ... as well as whether it was valid with cue; if `iso_sides` == `FALSE`, leaves `StimSides` unchanged
    group_by(CorrSide) %>%
    mutate(Side_Acc = mean(Acc, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(Side_Diff = max(Side_Acc) - min(Side_Acc)) %>%
    left_join(dem_df, by = c("participant" = "SubjID")) %>%
    filter(
      ifelse(attn_filter, "task" , "") %>% grepl(Attentiveness),            # Prunes participants reporting lack of alertness on at least two blocks, when `attn_filter` == TRUE
      Side_Diff <= side_bias,                                            #                     whose hit rate at one visual field - another visual field is > `side_bias`
      wm_Acc >= wm_floor) %>%                                            #                           working memory task accuracy is < `wm_floor`
    mutate(Acc_prefilter = mean(Acc, na.rm = TRUE)) %>%                       # Creates column indicating mean accuracy before we've filtered for `block_range`, unlike `Acc_postfilter`
    filter(
      Acc_prefilter %>% between( pre_range[1], pre_range[2], incbounds = TRUE), # Prunes participants whose non-catch, pre-block-filtering accuracy is outside of desired range
      between(CTI, min(CTI_range), max(CTI_range))) %>%                  #        trials outside of desired CTI range
    group_by(block) %>%
    mutate(block_acc = mean(Acc)) %>%                                         # Creates column indicating block's mean accuracy
    ungroup() %>%
    mutate(
      rown = row_number(),
      miniblock = RoundTo(rown, 16, ceiling) %>% `/` (16)) %>%                   #                           trial's mini-block (every 16 non-catch trials the opacity was readjusted)
    group_by(miniblock) %>%
    mutate(miniblock_avg = mean(Acc)) %>%
    ungroup() %>%
    dv_mutate(expr(ifelse(                                    # Changes `Acc` and `RT` (dv) column values to NA if...
      between(block_acc, block_range[1], block_range[2]) & between(           #  ... trial's block accuracy outside of `block_range`...
        miniblock_avg, miniblock_range[1], miniblock_range[2]),
      .,
      NA))) %>%  # ... or miniblock was not in the desired range or block was not in `blocks_desired`
    mutate(
      Trials_filtered_out = sum(is.na(Acc)) / n(),
      Acc_postfilter = mean(Acc, na.rm = TRUE)) %>%                      
    filter(
      between(Acc_postfilter, post_range[1], post_range[2]),             # Prunes participants whose non-catch, post-block-filtering accuracy is outside of desired range
      Trials_filtered_out <= filtered_cap)
  
  CTIs <- unique(cmbd$CTI)                                    # Calls `pcpts_combine` function for argument `pcpts` + combines each participant's dataframe into one
  if (win_func == "Tukey"){ win <- tukeywindow(length(CTIs), .5)} else {        
    win <- match.fun(paste0(tolower(win_func), "window"))(length(CTIs))}
  locations <- unique(cmbd$Stim_Sides)                                         
  pcpts <- unique(cmbd$participant)  
  
  empty_rows_per_pcpt <- 
    cmbd$CTI %>% range %>% diff %>% - duration %>%                         # Add empty rows (other than subject ID) as additional CTI's needed to reach desired padded...
    `/` (-samp_per) %>% + .00001 %>% floor %>% - 1                            # ... `duration` of intervals for each participant for each shuffle (`y` represents each shuffle)...
  
  # Determines confidence intervals
  conf_int <- function(x, y){ x %>%
      summarise_at(vars(y), list(~qnorm(.975) * std_err(.)))
  }
  
  amplitude <- function(x, shuf){
    set.seed(x)
    cmbd %>%     
      group_by(Stim_Sides, !!!grouping_cnsts) %>%
      dv_mutate(expr(case_when(shuf ~ sample(.), TRUE ~ .))) %>%
      group_by(CTI, Stim_Sides, !!!grouping_cnsts) %>%
      summarise_at(vars(!!dvcol), list(~mean(., na.rm = TRUE))) %>%               # Overwrites dv column according to mean of each combination of `CTI`, `Stim_Sides`, and `participant`
      dv_mutate(expr(na.approx(., na.rm = FALSE, rule = 2))) %>%
      dv_mutate(expr(rollapply(., clumps + 1, mean, partial = TRUE))) %>%
      group_by(participant, Stim_Sides) %>%
      dv_mutate(expr(na.approx(., na.rm = FALSE, rule = 2))) %>%
      ungroup() %>%
      mutate_at(vars(Stim_Sides), list(~case_when(!sbtr ~ ., TRUE ~ str_replace(
        ., regex("_?(In)?valid", ignore_case = T), "Invalid - Valid")))) %>%
      group_by(Stim_Sides, CTI, !!!grouping_cnsts) %>%
      summarise_at(vars(!!dvcol), list(~ifelse(sbtr, first(.) - last(.), .))) %>%
      group_by(participant, Stim_Sides) %>%
      dv_mutate(expr(case_when("Detrending" %in% trends ~
                                    . - (polyfit(CTI, ., 2) %>% polyval(CTI)),
                                  TRUE ~ .))) %>%
      dv_mutate(expr(case_when("Demeaning" %in% trends ~                    # Works just like the detrending except for demeaning
                                     . - mean(.), TRUE ~ .))) %>%
      dv_mutate(expr(. * win / Norm(win))) %>%               # Apply window
      ungroup() %>%
      add_row(
        participant = rep(pcpts, empty_rows_per_pcpt * (length(locations))),
        Stim_Sides = rep(locations, each = empty_rows_per_pcpt * length(pcpts)),
        !!dvcol := 0) %>%
      group_by(participant, Stim_Sides) %>%
      mutate(Power = fft(!!dvcol) %>%
               `*`(sqrt(2 / n())) %>%
               Mod %>% `^` (2)) %>%
      mutate(Hz = (rank(CTI) - 1) / (n() * samp_per)) %>%
      filter(
        rank(Hz) - 1 <= floor(n_distinct(Hz) / 2),                   # Remove alias frequencies above Nyquist
        Hz < xmax) %>%
      group_by(Hz, Stim_Sides) %>%
      summarise_at(vars(Power), mean)
  }
  
  tasco <- 1:shuff %>%
    map_dfr(~amplitude(.x, TRUE)) %>%
    group_by(Hz, Stim_Sides) %>%
    mutate(ntile = PercentRank(desc(Power)))
  
  fillerfunc <- function(Hz0, Stim_Sides0, Power0, xray0, yray0, zray0){
    tasco %>%                                                        # ... finds row in `amps_shuff` that...
      filter(
        Hz0 == Hz,
        Stim_Sides0 == Stim_Sides,
        !!parse_expr(yray0)) %>%                                            # ....3) with a smaller Power0 than `amps`' Power...
      top_n(1, Power) %>%                                                   # ... and 4) with largest Power0 among the remaining rows (we are being conservative and ntile reference...
      pull(xray0)  
  }
  
  amped_up <- map_dfr(1, ~amplitude(.x, FALSE)) %>%
    group_by(Stim_Sides) %>%
    mutate(
      ntile2 = 
        pmap_dbl(list(Hz0 = Hz, Stim_Sides0 = Stim_Sides, Power0 = Power,
                      xray0 = "ntile",
                      yray0 = "Power0 > Power | Power == min(Power)"),
                 fillerfunc),
      ntile2 =
        p.adjust(ntile2, method = mult_correcs) / ntile2,
      cutoff =
        map_dbl(list(Hz0 = Hz, Stim_Sides0 = Stim_Sides, Power0 = Power,
                     xray0 = "Power", yray0 = "ntile <  zray0 * (1 - α)",
                     zray0 = ntile2),
                fillerfunc))
  
  
  # Set Up Graphing
  
  t_srs_g <- function(x){
    cmbd_w %>%
      gather(Location, !!dvcol, -c(CTI, !!!grouping_cnsts)) %>%
      right_join(gather(conf_int(group_by(cmbd_w, CTI), locations),
                        Location, Conf_Int, -CTI),
                 by = c("CTI", "Location")) %>%
      group_by(CTI, Location, Conf_Int, !!!head(grouping_cnsts, x)) %>%
      summarise(!!dvcol := mean(!!dvcol)) %>%                     # Keeps either `RT` or `Acc` column—depending on whether `dvcol` parameter == `RT` or `Acc`
      ggplot(aes(CTI, !!dvcol, group = Location, color = Location,
                 fill = Location, ymin = !!dvcol - Conf_Int,
                 ymax = !!dvcol + Conf_Int)) + 
      labs(title = paste(dvcol, "by Cue-Target Interval, Exp.", dset),
           x = "Cue-Target Interval (ms)")
  }
  fft_g <- function(x) { x +
      labs(title = paste0("FFT of Target ", dep_var, ", Exp., ", dset),
           col = "Target Location", y = "Spectral Power") +
      theme(panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank())
  }
  
  graph <- function(y, x) {
    viridis_cols <- .7 +
      RoundTo(.0001 * RoundTo(length(locations), 4, floor), .2, ceiling)
    
    y(x) +
      theme_bw() +
      scale_color_viridis_d(option = "C",
                            end = viridis_cols,
                            labels = sapply(locations, simplify = TRUE,
                                            function(x) gsub("_", " ", x),
                                            USE.NAMES = FALSE)) +
      scale_fill_viridis_d(option = "C",
                           end = viridis_cols) +
      guides(colour = guide_legend(reverse = TRUE), fill = FALSE)
  }
  
  idvl_g <- function(x, y){
    graph(y, x) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)) +
      facet_wrap(~factor(participant, levels = pcpts),
                 ncol = round(length(pcpts) / 3), scales = "free_x")            # `free` means the y_axis isn't fixed from participant to participant
  }
  cmbd_g <- function(x, y) {
    graph(y, x) + geom_line(size = 1.5) +
      theme(legend.key.size = unit(.55, "in")) +
      labs(subtitle = paste("Data from", as.character(length(pcpts)),
                            "participants"))
  }
  
  sv_cmbd_g <- function(x) { x %>%
      ggsave(filename = file.path("plots", paste0(paste(display, dep_var, dset),
                                                  ".pdf")),
             width = xmax)}
  
  # Conditionally returns data frame of prelim participant or post-FFT data,
  # exits function
  if (display == "prelim_table") {
    write.csv(cmbd, file.path("plots", "Prelim_table.csv"))
  } else if (display == "fft_table") {
    write.csv(amps("Location", NULL), file.path("plots", "FFT_table.csv"))
  } else if (display == "Time-Series + FFT by Individual"){
    
    # Display Individuals' Plots -----------------------------------------------
    
    # Produces left half of final graph
    ts_facets <- idvl_g(1, t_srs_g) +
      geom_line(alpha = I(2 / 10), color = "grey", show.legend = FALSE) +       # Graphs unsmoothed data in light gray
      stat_smooth(method = tolower(smooth_method), span = 0.2, se = FALSE,      # Smoothes data depending on `smooth_method` parameter
                  size = .5, show.legend = FALSE)
    
    # Produces label for each right-side graph
    plot_label <- amps("Power0", NULL) %>%
      group_by(!!!grouping_cnsts) %>%
      summarise() %>%
      ungroup() %>%
      drop_na() %>%
      mutate_at(vars(!!!tail(grouping_cnsts, -1)),
                funs(paste(quo_name(quo(.)), "=", percent(.)))) %>%
      unite(lab, !!!tail(grouping_cnsts, -1), sep = "\n", remove = FALSE)
    
    # Produces right half of final graph
    fft_x <- round(1 / duration, 1)
    xaxis_r <- RoundTo(xmax, fft_x)
    
    fft_facets <- amps("Power0", c("participant", "Hz")) %>%
      summarise_all(mean) %>%
      gather(Flash_and_or_field, Power, -Hz, -samp_shuff,
             -c(!!!grouping_cnsts)) %>%
      ggplot(aes(Hz, Power, color = Flash_and_or_field)) %>%
      idvl_g(fft_g) +
      geom_line() +
      scale_x_continuous(
        name = "Frequency (Hz)", limits = c(0, xmax), breaks = seq(
          0, xmax, ifelse(xmax > 10 | xmax != xaxis_r,
                          xaxis_r %>% `/`(seq(fft_x, xaxis_r, fft_x)) %>%
                            Closest(5) %>% `/` (xaxis_r) %>% max %>% `^` (-1),
                          fft_x))) +
      labs(caption = paste("Data from", as.character(length(pcpts)),
                           "participants")) +
      geom_text(data = as.data.frame(plot_label), inherit.aes = FALSE,          # Sets location for label overlayed onto graph
                size = 1.2,
                aes(label = lab, x = Inf, y = Inf), vjust = 1.15, hjust = 1.05)
    
    side_by_side <- arrangeGrob(ts_facets, fft_facets, ncol = 2)                # Combines time series and FFT graphs into one plot
    ggsave(file.path("plots", "Indvls_Plots.pdf"), width = 25, side_by_side)
    
  } else if (display == "Time-Series Across Participants") {                    # Graph combined Time Series -----------
    
    (cmbd_g(0, t_srs_g) +
       theme(panel.grid = element_blank()) +
       geom_ribbon(alpha = 0.15, aes(color = NULL))) %>%
      move_layers("GeomRibbon", position = "bottom") %>%
      sv_cmbd_g
    
  } else { # Graph combined FFT ------------------------------------------------
    
    fft_x <- 1 / duration
    
    # Produces `shuff` # of null hypothesis permutations
    amps_shuff <- map_dfr(1:shuff, function(x){                                 # 1:shuff creates a set of CTI's*(length(pcpts)) rows for each number in this vector (before zero-padding)
      set.seed(x)                                                               # Set seed inside loop so each run of `shuffle` using a repeatable seed
      cmbd_w %>%
        group_by(participant) %>%                                               # The three lines randomize row order of CTI's for each participant; `group_by`, inherited from cmbd_w creation, is for  shuffling only...
        sample_n(n()) %>%                                                       # ... within rows with the same particpant number; `mutate_at` below allows shuffling only for CTI...
        mutate_at(vars(CTI), list(~seq(min(CTIs), max(CTIs), samp_per)))        # ... column; specifically, keeps rows intact except resets CTI's in descending order, even though...
    }) %>%                                                                      # ... previous line was just randomizing, thereby randomizing CTI vs. performance
      amplitude(shuff, "Location", "Power",                                     # Tabulates amplitude for each participant's data for `shuff` number of shuffles
                c("Hz", "samp_shuff", "Location")) %>%
      summarise_at(vars(Power), mean) %>%                                       # Finds average amplitude at each CTI across all participants per shuffle
      group_by(Hz, Location) %>%
      mutate_at(vars(Power), list(ntile = ~1.02 - ecdf(.)(.)))                  # Add columns converting null amplitudes -> a null distribution/rankings in the form of p-values
    
    
    amp_and_shf <- amps("Power0", c("Hz", "Location")) %>%
      summarise_at(vars(Power0), mean) %>%
      ungroup() %>%
      mutate(ntile = pmap_dbl(., function(Hz, Location, Power0, ...){           # For each row in `amps`...
        amps_shuff %>%                                                          # ... finds row in `amps_shuff` that...
          rename(Hz0 = Hz, Location0 = Location) %>%
          filter(Hz0 == Hz,                                                     # ... 1) matches in the Hz column...
                 Location0 == Location,                                         # ... 2) matches sides locations (e.g. `Invalid`)...
                 Power < Power0) %>%                                            # ....3) with a smaller Power0 than `amps`' Power...
          ungroup() %>%                                                         # ... and 4) with largest Power0 among the remaining rows (we are being conservative and ntile reference...
          top_n(1, Power0) %>%                                                  # .... Power0 is always a smaller power than would be Power)
          select(ntile) %>%
          as.double()
      })) %>%
      group_by(Location) %>%
      mutate_at(vars(ntile), list(~p.adjust(., method = mult_correcs) / .)) %>%
      select(Hz, Location, ntile) %>%
      safe_right_join(amps_shuff, by = c("Hz", "Location"),
                      conflict = `*`) %>%
      group_by(Hz, Location) %>%
      filter(ntile > 1 - α) %>% # at least one value will always be true - the last percentrank is equal to 1.00, which will always be greater than 1 - α
      top_n(1, desc(ntile)) # the smallest remaining ntile becomes the cutoff for the p-value test
    
    filter(ntile < α) %>% # new line # finds the largest percentile from null distribution that still is smaller than α
      top_n(1, ntile) %>% # new line (check that I don't need to use `desc` for how `ntile is getting sorted`)
      # filter(abs(ntile - α) == min(abs(ntile - α))) %>% # commented out line from before that is getting replaced with two lines above
      ungroup() %>%
      select(-c(samp_shuff, ntile)) %>%
      combine(amps("Power", c("Location", "Hz")) %>%                            # Finds average amplitude at each CTI for real (not surrogate) data, then merges that data with surrogate
                summarise(Power = mean(Power)) %>%
                ungroup(),
              names = c("Significance Cutoff", "Observed Data")) %>%
      right_join(conf_int(amps("Conf_Int", c("Hz", "Location")), "Conf_Int"),
                 by = c("Hz", "Location"))
    
    
    amp_and_shf2 <- amp_and_shf %>%
      ggplot(aes(Hz, Power, col = Location,
                 linetype = source, fill = Location,
                 ymin = Power - Conf_Int,
                 ymax = Power + Conf_Int)) %>%
      cmbd_g(fft_g) +
      scale_linetype_manual(values = c("solid", "dashed")) +
      scale_x_continuous(name = "Frequency (Hz)",
                         limits = c(0, xmax),
                         breaks = seq(0, xmax, ifelse(fft_x > .5,
                                                      round(fft_x, 2), 1))) +
      labs(linetype = "",
           caption = paste("Significance threshold at p < ",
                           as.character(α))) +
      geom_ribbon(data = filter(amp_and_shf, source == "Observed Data"),
                  alpha = 0.15, aes(color = NULL)) +
      geom_point(size = 3,
                 data = amp_and_shf %>%
                   spread(source, Power) %>%
                   filter(`Observed Data` > `Significance Cutoff`) %>%
                   select(-c(`Significance Cutoff`, Conf_Int)) %>%
                   gather(source, Power, -Hz, -Location), 
                 aes(ymin = NULL, ymax = NULL))
    move_layers(amp_and_shf2, "GeomRibbon", position = "bottom") %>%
      sv_cmbd_g
  }
}



# Sets inputs for the `main_function` function
main_function(display = "Time-Series Across Participants",                              # Either `FFT Across Participants`, `Time-Series Across Participants`, `Time-Series + FFT by...
              # ... Individual`, `prelim_table` (lightly analyzed data), and `fft_table` (semi-ready for graphing data)
              
              dset = "1a",                                                      # Either `1a`, `1b`, `2a`, `2b`, `2c`, or `3a`
              
              wm_exp = FALSE,                                                   # Either `FALSE` or `TRUE` for working memory participants
              
              iso_sides = FALSE,                                                # Either `FALSE` or `TRUE`, which groups by not only valid and invalid but also by the side of the...
              # ... screen for each trial (i.e. going from `Valid` and `Invalid` to `Right Valid`, `Left Valid`,... 
              # ... `Right Invalid`, and `Left Invalid`)
              
              sbtr = FALSE,                                                     # Either `FALSE` or `TRUE`, which subtracts the dependent variable values at each CTI (valid - invalid)...
              # ... before performing analyses rather than analyzing valid and invalid trials independently
              
              samp_per = 3/60,                                                  # Spacing between CTI intevals (in seconds); the data was originally sampled at 1 / 60, but one could...
              # ... re-sample at a different rate, which would just clump neighboring CTI's together (whereas the ...
              # ... `clump` variable groups neighbors but doesn't combine them, keeping the same total number of bins)
              
              clumps = 0,                                                       # Number of points to average at each CTI; `1` means this function does nothing, `3` means each CTI is...
              # ... the average of that CTI and its neighboring CTI's on each sides, etc...
              
              dep_var = "Accuracy",                                             # Either `Accuracy`or `Response Time`
              
              α = .05,                                                          # The alpha threshold to use for drawing the significance cutoff on the graphs
              
              shuff = 50,                                                       # The number of surrogate shuffles to use to determine the null hypothesis; NOTE: increasing this number...
              # ... slows down the run time
              
              mult_correcs = "BH",                                              # Choose between the following types of multiple corrections: `holm`, `hochberg`, `hommel`, ...  
              # ... `bonferroni`, `BH`, `BY`, `fdr`, and `none`
              
              trends = c("Detrending", "Demeaning"),                            # Either `Detrending`, `Demeaning`, both, or an empty vector
              
              smooth_method = "Loess",                                          # Either `Loess`, `LM`, `GLM`, or `GAM` smoothing methods for the time-series data of individuals
              
              win_func = "Tukey",                                               # Choose between the following types of windowing functions: `Tukey`, `Square`, `Hann`, `Welch`, ...
              # ... `Triangle`, `Hamming`, `Cosine`, or `Kaiser`
              
              xmax = 15,                                                        # Greatest x-axis value included in graph
              
              duration = 1,                                                     # Duration (Seconds) Analyzed Including Padding
              
              attn_filter = FALSE,                                              # Either `FALSE` or `TRUE`, which prunes any participants who indicated in the post-task questionnaire...
              # ... that they either dozed off at one point or were not fully focused on at least two of the eight...
              # ... blocks (the other options were that they were fully alert on all blocks or fully alert on all but...
              # ... one block)
              
              catch_floor = .85,                                                # Prunes participants who perform below this hit rate on catch trials, which featured no target
              
              side_bias = .2,                                                   #                     whose hit rate at one visual field - another visual field is > `side_bias`
              
              wm_floor = .7,                                                    #                           working memory task accuracy is < `wm_floor`
              
              invalid_floor = .02,                                              #                           invalid trial accuracy is < `invalid_floor`
              
              pre_range = c(.45, .85),                                          #                           unfiltered/scrutinized data is outside of the selected range
              
              post_range = c(.45, .85),                                         #                           filtered data is outside of the selected range
              
              filtered_cap = .4,                                                #                     with at least `filtered_cap` % of trials filtered out
              
              block_range = c(.40, .80),                                        # Interpolates over trials if the average hit rate in that block, every 48 trials, is outside of the...
              # ... select range
              
              blocks_desired = 1:8,                                             #                                 block is not in the `blocks_desired`
              
              miniblock_range = c(.20, .80),                                    #                                 average hit rate in that mini-block, every 16 trials which is how...
              # ... often the task difficulty was adjusted to titrate to 65%, is outside this range
              
              CTI_range = c(.3, 1.09)                                           # Filters trials outside this range of CTI bins
)
