# Install packages if not already installed, then load them
if (!require(devtools)) install.packages('devtools')
if (!require(smisc)) devtools::install_github("stevenworthington/smisc")
smisc::ipak(c("utils", "tidyr", "dplyr", "ggplot2", "DescTools", "e1071",
              "pracma", "gridExtra", "data.table", "tables", "zoo",
              "parallel", "scales"))

# Main function begins (encompasses other functions)
data_graphing <- function(catch_trial_cutoff, block_floor, mini_block_floor,
                          mini_block_ceil, catch, lowacc_prefilter,
                          highacc_prefilter, lowacc_postfilter,
                          highacc_postfilter, sep_vis_fields, smooth_method,
                          dep_var, sampling_freq, latestart, earlyend,
                          ptcpts, output, save_output, clumping){
  
  grouping_constants <- quos(participant, Trials_filtered_out, Acc_prefilter,   # Columns that are frequently used for grouping, variable means
                             Acc_postfilter, CatchAcc)                          # don't have to type them out every time we use them for grouping

  dep_var <- as.name(dep_var)                                                   # Converts character to name/symbol so we can refer to it as a column using tidy

  # For each participant, function reads in data, filters it, and transforms it to prepare for interpolation and FFT'ing
  ptcpt_tabulate <- function(participant_num){
    ptcpt_path <- file.path("data", participant_num,                            # Creates file path for reading in participant data,
                            paste0(participant_num, ".csv"))                    # accomodating all operating systems
    fread(ptcpt_path, select = c(1:23)) %>%                                     # Reads in participant data (the 'select' part is because PsychoPy created an empty column at the end of the data frame for the first few participants, which meant those dataframes had different dimensions than subsequent data frames, which 'do.call' doesn't like)
      filter(Trial > 0) %>%                                                     # Filters out practice trials
      mutate(CatchAcc = mean(ifelse(Opacity > 0, NA, Acc), na.rm = TRUE)) %>%   # Creates column indicating mean accuracy for catch trials
      filter(CatchAcc >= catch_trial_cutoff,                                    # Filters out participants whose catch accuracy is below desired threshold
             Opacity > catch) %>%                                               #             catch trials if 'catch' parameter is assigned to 0; if it's assigned to -1, this line does nothing)
      mutate(Acc_prefilter = mean(Acc, na.rm = TRUE)) %>%                       # Creates column indicating mean accuracy before we've filtered for 'block_floor', unlike 'Acc_postfilter'
      filter(between(Acc_prefilter, lowacc_prefilter, highacc_prefilter)) %>%   # Filters out participants whose non-catch, pre-block-filtering accuracy is outside of desired range                                            #             catch trials if 'catch' parameter is assigned to 0; if it's assigned to -1, this line does nothing)
      mutate(block = RoundTo(Trial, 54, ceiling) / 54,                          # Creates column indicating trial's block
             CTI = RoundTo(RoundTo(lilsquareStartTime - flash_circleEndTime,
                                   1 / 60), sampling_freq),
             RT = ifelse(Acc == 1 & ButtonPressTime - lilsquareStartTime > .1,  #                indicating RT after target appeared on screen, only for correct trials with an RT < 100 ms
                         ButtonPressTime - lilsquareStartTime, NA),
             Stim_Sides = as.character(
               ifelse(CorrSide == FlashSide, "Congruent", "Incongruent")),      #                indicating whether cue was congruent or incongruent
             CorrSide = if_else(CorrSide == 1, "Right", "Left"),                #                indicating which side the target appeared on
             Stim_Sides = case_when(sep_vis_fields == "No" ~ Stim_Sides,        # Overwrites 'Stim_Sides' column if 'sep_vis_fidels' parameter == 'Yes' to include which side of screen target was on,
                TRUE ~ paste(Stim_Sides, CorrSide, sep = "_"))) %>%             # as well as whether it was congruent with cue; if 'sep_vis_fields' parameter == 'No', leaves 'StimSides' unchanged
      group_by(block) %>%
      mutate(block_acc = mean(Acc)) %>%                                         # Creates column indicating block's mean accuracy
      group_by(Opacity > 0) %>%
      mutate(rown = row_number(),
             miniblock = ifelse(Opacity > 0,                                    # Creates column indicating trial's mini-block (every 16 trials the opacity was readjusted)
                                RoundTo(rown, 16, ceiling) / 16, NA)) %>%
      group_by(miniblock) %>%
      mutate(miniblock_avg = mean(Acc)) %>%
      ungroup() %>%
      mutate_at(vars(Acc, RT), funs(ifelse(block_acc <= block_floor |           # Changes 'Acc' and 'RT' column values to NA if trial's block accuracy < 'block_floor'
                                           !between(miniblock_avg,
                                                    mini_block_floor,           # or mini_block was not in the desired range
                                                    mini_block_ceil),
                                           NA, .))) %>%              
      mutate(Trials_filtered_out = sum(is.na(Acc)) / n(),
             Acc_postfilter = mean(Acc, na.rm = TRUE)) %>%                      # Creates column indicating mean accuracy for non-catch trials; note this is after before we've filtered for block_floor, unlike Acc_prefilter
      filter(between(Acc_postfilter, lowacc_postfilter, highacc_postfilter)) %>%# Filters out participants whose non-catch, post-block-filtering accuracy is outside of desired range 
      group_by(CTI, Stim_Sides, !!!grouping_constants) %>%
      summarise_at(vars(Acc, RT), funs(mean(., na.rm = TRUE))) %>%              # Overwrites 'Acc' and 'RT' columns according to mean of each combination of 'CTI' and 'Stim_Sides'
      arrange(Stim_Sides, CTI) %>%
      ungroup() %>%
      mutate_at(vars(Acc, RT), funs(na.approx(., na.rm = FALSE, rule = 2))) %>%
      mutate_at(vars(Acc, RT), funs(rollapply(.,
                                              clumping, mean, partial = TRUE)))
    }
  
  ptcpts_tabulated <- do.call(rbind, mclapply(ptcpts, ptcpt_tabulate)) %>%      # Calls 'ptcpt_tabulate' function for argumenet 'ptcpts'; then combines each participant's dataframe into one
                        arrange(Acc_prefilter)
  ptcpts_remaining <- unique(ptcpts_tabulated$participant)                      # Creates vector of remaining participant numbers after 'ptcpt_tabulate' filtering
  hann <- hanning.window(length(unique(ptcpts_tabulated$CTI)))                  # Creates Hanning window, gets applied later
  Sides <- unique(ptcpts_tabulated$Stim_Sides)
  
  amplitudes <- ptcpts_tabulated %>%
    pivot_wide(CTI:CatchAcc, Stim_Sides, !!dep_var) %>%
    group_by(participant) %>%
    mutate_at(vars(!!Sides), funs(Mod(fft(detrend(.) * hann)))) %>%             # Detrends, multiplies by Hanning window, applies FFT, and then takes magnitude   
    mutate(Hz = (row_number()-1) / (n() * sampling_freq)) %>%
    ungroup() %>%
    select(-CTI)
  
 ### Outputs ####
  
  # Conditionally returns data frame of prelim participant or post-FFT data, exits function
  if (output == "prelim_table") {
    if (save_output == "Yes") {                                                 # Conditionally saves data frame as .csv in analysis folder
      write.csv(ptcpts_tabulated, file.path("analysis", "Prelim_table.csv"))
    }
    return(ptcpts_tabulated)
  } else if (output == "fft_table") {
    if (save_output == "Yes") {                                                 # Conditionally saves data frame as .csv in analysis folder
      write.csv(amplitudes, file.path("analysis", "FFT_table.csv"))
    }
    return(amplitudes)
  }          
  
  #Graphing
  
  # Produces left half of the final graph
  time_series_facets <- rbind((ptcpts_tabulated %>%                             # Copies each row in 'ptcpts_tabulated' except 'participant' = 'All', which creates the top graph on the left labeled 'All'
                                 mutate(participant = "All")),
                              ptcpts_tabulated) %>%
    group_by(CTI, Stim_Sides, !!!grouping_constants) %>%
    summarise(!!dep_var := mean(!!dep_var)) %>%                                 # Keeps either RT or Acc column—depending on whether 'dep_var' parameter == 'RT' or 'Acc'
    ggplot(aes(CTI, !!dep_var, group = Stim_Sides, color = Stim_Sides)) +       # CTI is x axis, dep_var is y axis
    geom_line(alpha = I(2 / 10), color = "grey", show.legend = FALSE) +         # Graphs unsmoothed data in light gray
    stat_smooth(method = smooth_method, span = 0.2, se = FALSE,                 # Smoothes data depending on 'smooth_method' parameter
                size = .5, show.legend = FALSE) +
    labs(title = paste0(dep_var, " by Cue-Target Interval, ", smooth_method,
                        "-Smoothed"),
         x = "Cue-Target Interval (ms)",
         caption = paste("Not fft'ed;", 120 * sampling_freq,                    # 2 * 60 = 120; the more clumping, the more we multiply the 2 data points by
                         "data points/bin/target side/participant\n",
                         as.character(length(ptcpts_remaining)),
                         "participants averaged")) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    facet_wrap(~ factor(participant, levels = c("All", ptcpts_remaining)),
               ncol = 4, scales = 'free_x')
  
  # Creates label each right-side graph
  graph_label <- amplitudes %>%
    group_by(!!!grouping_constants) %>%
    summarise() %>%
    ungroup() %>%
    mutate_at(vars(!!!tail(grouping_constants,-1)),
              funs(paste(quo_name(quo(.)),"=", percent(.)))) %>%
    unite(lab, !!!tail(grouping_constants,-1), sep = "\n", remove = FALSE)
  
  # Produces right half of final graph
  fft_facets <- rbind((amplitudes %>% mutate(participant = "All")),             # Copies each row in 'amplitudes' except 'participant' = 'All', which creates the top graph on the right labeled 'All'
                      amplitudes) %>%
    group_by(participant, Hz) %>%
    summarise_all(mean) %>%
    gather(Flash_and_or_field, Magnitude, -Hz, -c(!!!grouping_constants)) %>%
    ggplot(aes(Hz, Magnitude, color = Flash_and_or_field)) +
    geom_line() +
    scale_x_continuous(breaks = c(4, 8),
                       limits = c(0, max(amplitudes$Hz) / 2)) +                 # Creates thick vertical lines at 4 Hz and 8 Hz, and sets x axis range to 0-16
    labs(title = paste("FFT of", as.character(dep_var)),
         caption = paste("Data from", as.character(length(ptcpts_remaining)),
                         "participants")) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid.major = element_line(colour = "white", size = 3)) +
    facet_wrap(~ factor(participant, levels = c("All", ptcpts_remaining)),
               ncol = 4, scales = 'free') +                                     # 'free' means the y_axis isn't fixed from participant to participant
    geom_text(data = as.data.frame(graph_label), vjust = 1, hjust = 1,          # Sets location for label overlayed onto graph
              aes(label = lab, x = Inf, y = Inf), inherit.aes = FALSE, size = 1)

  
  side_by_side <- arrangeGrob(time_series_facets, fft_facets, ncol = 2)         # Combines time series and FFT graphs into one plot
  ggsave(file.path("analysis", "Time-Series_+_FFT_Plots.pdf"), width = 25,
         side_by_side)

}



# Sets inputs for the 'data_graphing' function
data_graphing(catch_trial_cutoff = .85,                                         # Filters out participants with a catch trial accuracy below this value 
              block_floor = .25,                                                # Converts a trial's accuracy to NA if its block's accuracy below this value
              mini_block_floor = .35,                                           # Converts a trial's accuracy to NA if its mini-block's accuracy below this value
              mini_block_ceil = .9,                                             # Converts a trial's accuracy to NA if its mini-block's accuracy above this value
              catch = 0,                                                        # -1/0 to include/exclude catch trials in analyses
              lowacc_prefilter = .35,                                           # Filters participants by their accuracy before their blocks below 'block_floor' have been interpolated over
              highacc_prefilter = .9,                                           # Line above is the floor of the filter, this line is the ceiling
              lowacc_postfilter = .35,                                          # Filters participants by their accuracy after their blocks below 'block_floor' have been interpolated over 
              highacc_postfilter = .9,                                          # Line above is the floor of the filter, this line is the ceiling
              sep_vis_fields = "No",                                            # Use 'No' and 'Yes'; 'No' means incongruent cue and target trials are analyzed differently at each CTI than trials with congruent cue and target; 'Yes' means we further divide analyses by which side of screen target appeared, yielding four conditions (congruent/incongruent and left/right)
              smooth_method = "loess",                                          # See geom_smooth documentation for available smoothing methods
              dep_var = "Acc",                                                  # Use 'Acc' or 'RT'
              sampling_freq = 1 / 60,                                           # Equals spacing between CTI intevals (in seconds)
              latestart = 0,  earlyend = 0,                                     # Filters out, for FFT analysis, trials with a CTI < 'latestart' or a CTI > ((largest CTI [so 1.3]) - 'earlyend'); units = seconds
              ptcpts = 201:230,
              output = "graph",                                                 # Use "graph", "prelim_table" (before interpolation and FFT'ing), and "fft_table"
              save_output = "Yes",                                              # Use "Yes" and "No" (if 'no', table outputs will still be visible in R)
              clumping = 3)                                                     # Use 1 (no clumping) and 3 (each bin is the average of itself and its two neighbors)
  