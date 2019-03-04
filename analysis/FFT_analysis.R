# Install packages if not already installed, then load them
if (!require(devtools)) install.packages('devtools')
devtools::install_github("stevenworthington/smisc")
smisc::ipak(c("utils", "tidyr", "dplyr", "ggplot2", "DescTools",
              "e1071", "pracma", "gridExtra", "data.table",
              "tables", "zoo", "tidyverse", "parallel", "scales"))

# Main function begins (encompasses other functions)                            # Can take several few minutes to run, especially on Windows
data_graphing <- function(catch_trial_cutoff,  block_acc_cutoff, catch,
                          lowacc_ptcpts_pre_block_filter,
                          highacc_ptcpts_pre_block_filter,
                          lowacc_ptcpts_post_block_filter,
                          highacc_ptcpts_post_block_filter,
                          jitter, sep_vis_fields, smooth_method,  dep_var,
                          sampling_freq, win_freq, win_size, latestart,
                          earlyend, ptcpts, output, save_output){
  
  grouping_constants <- quos(participant, Trials_below_block_cutoff,            # Columns that are frequently used for grouping,
                             Acc_blocks_unfiltered, Acc_blocks_filtered,        # variable means don't have to type them out every
                             CatchAcc)                                          # time we use them for grouping
  dep_var <- as.name(dep_var)                                                   # Converts character to name/symbol so we can refer to it as a column using tidy

  # For each participant, function reads in data, filters it,
  # and transforms it to prepare for interpolation and FFT'ing
  ptcpt_tabulate <- function(participant_num){
    ptcpt_path <- file.path("data", participant_num,                            # Creates file path for reading in participant data,
                            paste0(participant_num, ".csv"))                    # accomodating all operating systems
    rawdata <- fread(ptcpt_path, select = c(1:23))                              # Reads in participant data (the 'select' part is because PsychoPy created an empty column at the end of the data frame for the first few participants, which meant those dataframes had different dimensions than subsequent data frames, which 'do.call' doesn't like)
    rawfactors <- sort(unique((rawdata %>% filter(Trial > 0))$Opacity))
    rawdata %>%
      mutate(ExpAcc = ifelse(Opacity > 0, Acc, NA),                             # Creates column indicating accuracy for non-catch trials, and indicating NA for catch trials
             CatchAcc = ifelse(Opacity > 0, NA, Acc),                           #                indicating accuracy for catch trials, and indicating NA for non-catch trials
             Acc_blocks_unfiltered = mean(ExpAcc, na.rm = TRUE),                #                indicating mean accuracy for non-catch trials; note this is created before we've filtered for 'block_acc_cutoff', unlike 'Acc_blocks_filtered'
             CatchAcc = mean(CatchAcc, na.rm = TRUE)) %>%                       #                indicating mean accuracy for catch trials
      filter(between(Acc_blocks_unfiltered, lowacc_ptcpts_pre_block_filter,     # Filters out participants whose non-catch, pre-block-filtering accuracy is outside of desired range
                     highacc_ptcpts_pre_block_filter),
             CatchAcc > catch_trial_cutoff,                                     #             participants whose catch accuracy is below desired threshold
             Trial > 0,                                                         #             practice trials
             Opacity > catch) %>%                                               #             catch trials if 'catch' parameter is assigned to 0; if it's assigned to -1, this line does nothing)
      mutate(Acc = ifelse(((Opacity > 0 &                                       # Overwrites 'Acc' column, since original PsychoPy code was written to count any trial with
                              ((Key == 'l' & CorrSide == 1) |                   # a false press as incorrect; this code won't change anything for newest code
                                Key == 'a' & CorrSide == -1)) |
                           (Opacity == 0 & Key == 'nope')), 1, 0),
             block = RoundTo(Trial, 44, ceiling)/44,                            # Creates column indicating trial's block
             CTI = RoundTo(lilsquareStartTime - flash_circleEndTime,
                                   sampling_freq),
             RT = ifelse(Acc == 1 | ButtonPressTime - lilsquareStartTime > .1,  #                indicating RT after target appeared on screen, only for correct trials with an RT < 100 ms 
                         ButtonPressTime - lilsquareStartTime, NA),   
             Stim_Sides = as.character(                                         
               ifelse(CorrSide == FlashSide, "Congruent", "Incongruent")),      #                indicating whether cue was congruent or incongruent
             CorrSide = if_else(CorrSide == 1, "Right", "Left"),                #                indicating which side the target appeared on
             Stim_Sides = case_when(sep_vis_fields == "No" ~ Stim_Sides,        # Overwrites 'Stim_Sides' column if 'sep_vis_fidels' parameter == 'No' to include which side of screen target was on,
                TRUE ~ paste(Stim_Sides, CorrSide, sep = "_"))) %>%             # as well as whether it was congruent with cue; if 'sep_vis_fidels' parameter == 'Yes', leaves 'StimSides' unchanged
      group_by(block) %>%
      mutate(block_acc = mean(Acc)) %>%                                         # Creates column indicating block's mean accuracy
      ungroup() %>%
      mutate_at(vars(Acc, RT),                                                  # Changes 'Acc' and 'RT' column values to NA if trial's block accuracy < 'block_acc_cutoff'
             funs(ifelse(block_acc < block_acc_cutoff, NA, .))) %>%
      mutate(Trials_below_block_cutoff = sum(is.na(Acc))/n()) %>%               # Creates column indicating proportion of all trials in a block below 'block_acc_cutoff'
      mutate_at(vars(Acc, RT),
                funs(ifelse(Opacity %in% rawfactors[jitter + 4], ., NA))) %>%   # Changes 'Acc' and 'RT' column values to NA if trial's jitter level != 
      mutate(Acc_blocks_filtered = mean(Acc, na.rm = TRUE)) %>%                 # Creates column indicating mean accuracy for non-catch trials; note this is after before we've filtered for block_acc_cutoff, unlike Acc_blocks_unfiltered
      filter(between(Acc_blocks_filtered, lowacc_ptcpts_post_block_filter,      # Filters out participants whose non-catch, 
                     highacc_ptcpts_post_block_filter)) %>%                     # post-block-filtering accuracy is outside of desired range
      group_by(CTI, Stim_Sides, !!!grouping_constants) %>%
      summarise_at(vars(Acc, RT), funs(mean(., na.rm = TRUE))) %>%              # Overwrites 'Acc' and 'RT' columns according to mean of each combination of 'CTI' and 'Stim_Sides'
      arrange(Stim_Sides, CTI) %>%
      ungroup() %>%
      mutate_at(vars(Acc, RT), funs(na.approx(., na.rm = FALSE, rule = 2)))     # 'CTI's with an avg of NaN take the mean of the 'CTI' accuracies before and after it
    }
  
  ptcpts_tabulated <- do.call(rbind,lapply(ptcpts, ptcpt_tabulate))             # Calls 'ptcpt_tabulate' function for argumenet 'ptcpts'; then combines each participant's dataframe into one
  
  # Conditionally returns data frame of prelim participant data, exits function
  if (output == "prelim_table") {
    if (save_output == "Yes") {                                                 # Conditionally saves data frame as .csv in analysis folder
       write.csv(ptcpts_tabulated, file.path("analysis", "Prelim_table.csv"))
      }
    return(ptcpts_tabulated)
    }                      

  ptcpts_remaining <- unique(ptcpts_tabulated$participant)                      # Creates vector of remaining participant numbers after 'ptcpt_tabulate' filtering
  ptcpts_tabulated <- rbind((ptcpts_tabulated %>%                               # Copies each row in ptcpts_tabulated, but then changes 'participant' column to 'All;
                               mutate(participant = "All")), ptcpts_tabulated)  # this is going to be helpful for the graph
  
  # Produces left half of the final graph
  time_series_facets <- rbind((ptcpts_tabulated %>%                             # Copies each row in 'ptcpts_tabulated' except 'participant' = 'All', which creates the top graph on the left labeled 'All'
                                mutate(participant = "All")),
                              ptcpts_tabulated) %>%
    group_by(CTI, Stim_Sides, !!!grouping_constants) %>%
    summarise(!!dep_var := mean(!!dep_var)) %>%                                 # Keeps either RT or Acc column—depending on whether 'dep_var' parameter == 'RT' or 'Acc'
    ggplot(aes(CTI, !!dep_var, group = Stim_Sides, color = Stim_Sides)) +       # CTI is x axis, dep_var is y axis
      geom_line(alpha = I(2/10), color = "grey", show.legend = FALSE) +         # Graphs unsmoothed data in light gray
      stat_smooth(size = 1, method = smooth_method, span = 0.2, se = FALSE,     # Smoothes data depending on 'smooth_method' parameter
                  show.legend = FALSE) +
      labs(title = paste0(dep_var, " by Cue-Target Interval, ", smooth_method,
                          "-Smoothed"),
           x = "Cue-Target Interval (ms)",
           caption = paste("Not fft'ed; 2 data points per x-value per target side per participant\n",
                           as.character(length(ptcpts_remaining)),
                           "participants averaged")) +
        theme(plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5)) +
    facet_wrap(~ factor(participant, levels = c("All", ptcpts_remaining)),
               ncol = 1, scales = 'free_x')
  
  
  # Interpolates + Detrends + Hanning Window + FFT
  
  # Creates a sliding window for interpolation
  win_size <- win_size - .001                                                   # -.001 is because the between function between is inclusive
  first_win <- min(ptcpts_tabulated$CTI) + latestart
  last_win <- max(ptcpts_tabulated$CTI) - earlyend
  wins <- seq(first_win, last_win, win_freq)
  hann <- hanning.window(length(wins))                                          # Creates Hanning window, gets applied later
  
  # Interpolates and returns amplitude for each participant (assuming they haven't already been filtered out above)
  amplitude_func <- function(participant_num){
   exp_trials_noncatch <- ptcpts_tabulated %>%
     filter(participant == participant_num)
   Sides <- levels(factor(exp_trials_noncatch$Stim_Sides))
   
    # Interpolates as a sliding window
    slidingwindow <- function(win_start){
      win_end <- win_start + win_size
      exp_trials_noncatch %>%
        filter(between(CTI, win_start, win_end)) %>%
        group_by(Stim_Sides, !!!grouping_constants) %>%
        summarise(!!dep_var := mean(!!dep_var)) %>%
        spread(Stim_Sides, !!dep_var)
    }
    as.data.frame(t(mcmapply(slidingwindow, wins))) %>%
      modify(as.numeric) %>%
  	  mutate_at(vars(!!Sides), funs(Mod(fft(detrend(.)*hann)))) %>%             # Detrends, multiplies by Hanning window, applies FFT, and then takes magnitude
      modify(as.numeric) %>%
      mutate(Hz = (row_number()-1)/(n()*win_freq))
  }
  
  # Calls 'amplitude_func' for all remaining participants, then combine the data into one data frame; the interpolating is what makes this line slow the whole script
  amplitudes <- do.call(rbind, mclapply(ptcpts_remaining, amplitude_func))
  
  # Creates label each graph
  graph_label <- amplitudes %>%
    group_by(!!!grouping_constants) %>%
    summarise() %>%
    ungroup() %>%
    mutate_at(vars(!!!tail(grouping_constants,-1)),
              funs(paste(quo_name(quo(.)),"=", percent(.)))) %>%
    unite(lab, !!!tail(grouping_constants,-1), sep = "\n", remove = FALSE)
  
  # Returns data frame of post-FFT participant data,
  # exits function if "output == "fft_table""
  if (output == "fft_table") {
    if (save_output == "Yes") {                                                 # Conditionally saves data frame as .csv in analysis folder
      write.csv(amplitudes, file.path("analysis", "FFT_table.csv"))
    }
    return(amplitudes)
  }          
  
  # Produces right half of final graph
  fft_facets <- rbind((amplitudes %>% mutate(participant = "All")),             # Copies each row in 'amplitudes' except 'participant' = 'All', which creates the top graph on the right labeled 'All'
                      amplitudes) %>%
    group_by(participant, Hz) %>%
    summarise_all(mean) %>%
    gather(Flash_and_or_field, Magnitude, -Hz, -c(!!!grouping_constants)) %>%
    ggplot(aes(Hz, Magnitude, color = Flash_and_or_field)) +
    geom_line() +
    scale_x_continuous(breaks = c(4, 8), limits = c(0, 16)) +                   # Creates thick vertical lines at 4 Hz and 8 Hz, and sets x axis range to 0-16
    labs(title = paste("FFT of", as.character(dep_var)),
         caption = paste("Data from", as.character(length(ptcpts_remaining)),
                         "participants")) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid.major = element_line(colour = "white", size = 3)) +
    facet_wrap(~ factor(participant, levels = c("All", ptcpts_remaining)),
               ncol = 1, scales = 'free') +                                     # 'free' means the y_axis isn't fixed from participant to participant
    geom_text(data = as.data.frame(graph_label), vjust = 1, hjust = 1,          # Sets location for label overlayed onto graph
              aes(label = lab, x = Inf, y = Inf), inherit.aes = FALSE)
 
  # Saves graph to 'analysis' folder
  side_by_side <- arrangeGrob(time_series_facets, fft_facets, ncol = 2)     # Combines time series and FFT graphs into one plot
  ggsave(file.path("analysis", "Time-Series_+_FFT_Plots_.pdf"),
         side_by_side, limitsize = FALSE, width = 25,
         height = length(ptcpts_remaining) * 3 + 5)
}



# Sets inputs for the 'data_graphing' function
data_graphing(catch_trial_cutoff = .85,                                         # Filters out participants with a catch trial accuracy below this value 
              block_acc_cutoff = .25,                                           # Converts a trial's accuracy to NA if block's accuracy below this value
              catch = 0,                                                        # -1/0 to include/exclude catch trials in analyses
              lowacc_ptcpts_pre_block_filter = .35,                             # Filters participants by their accuracy before their blocks below 'block_acc_cutoff' have been interpolated over
              highacc_ptcpts_pre_block_filter = .9,                              # Line above is the floor of the filter, this line is the ceiling
              lowacc_ptcpts_post_block_filter = .35,                            # Filters participants by their accuracy after their blocks below 'block_acc_cutoff' have been interpolated over 
              highacc_ptcpts_post_block_filter = .9,                            # Line above is the floor of the filter, this line is the ceiling
              jitter = c(-2:2),                                                 # Use "-2:2" to keep all jitter levels; or -2, -1, 0, 1, and/or 2 to refer to jitter level, separate non-consecutive jitters by commas (negative means below staircase threshold, 0 means at staircase threshold, positive means above)
              sep_vis_fields = "No",                                            # Use 'No' and 'Yes'; 'No' means incongruent cue and target trials are analyzed differently at each CTI than trials with congruent cue and target; 'Yes' means we further divide analyses by which side of screen target appeared, yielding four conditions (congruent/incongruent and left/right)
              smooth_method = "loess",                                          # See geom_smooth documentation for available smoothing methods
              dep_var = "Acc",                                                  # Use 'Acc' or 'RT'
              sampling_freq = 1/60,                                             # Equals spacing between CTI intevals (in seconds)
              win_freq = .001, win_size = .05,                                  # Set parameters for sliding window (in seconds)
              latestart = 0,  earlyend = 0,                                     # Filters out, for FFT analysis, trials with a CTI < 'latestart' or a CTI > ((largest CTI [so 1.836]) - 'earlyend'); units = seconds
              ptcpts = 201:230,                                                 # Error w/ PsychoPy coding for ptcpts 201-205
              output = "graph",                                                 # Use "graph", "prelim_table" (before interpolation and FFT'ing), and "fft_table"
              save_output = "Yes")                                              # Use "Yes" and "No" (if 'no', table outputs will still be visible in R)
