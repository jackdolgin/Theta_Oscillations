if (!require(devtools)) install.packages('devtools')
devtools::install_github("stevenworthington/smisc")
smisc::ipak(c("utils", "tidyr", "dplyr", "ggplot2", "DescTools", "e1071", "pracma", "gridExtra", "data.table", "tables", "zoo", "tidyverse", "parallel"))

data_graphing <- function(lowacc_ptcpts, highacc_ptcpts, catch_trial_cutoff,
                          block_acc_cutoff, catch, sep_vis_fields, dep_var,
                          sampling_freq, win_freq, win_size, latestart, earlyend, ptcpts, save){
  
  grouping_constants <- quos(participant, Trials_below_block_cutoff, ExpAccmeanbefore, ExpAccmeanafter, CatchAccmean)
  label_constants <- quos(Trials_below_block_cutoff, ExpAccmeanbefore, ExpAccmeanafter, CatchAccmean)
  dep_var <- as.name(dep_var)
  
  ptcpt_tabulate <- function(participant_num){
    ptcpt_path <- file.path("data", participant_num, paste0(participant_num, ".csv"))
    fread(ptcpt_path, select = c(1:23)) %>%
      mutate(ExpAcc = ifelse(Opacity > 0, Acc, NA),
             CatchAcc = ifelse(Opacity > 0, NA, Acc),
             ExpAccmeanbefore = mean(ExpAcc, na.rm = TRUE),
             CatchAccmean = mean(CatchAcc, na.rm = TRUE)) %>%
      filter(Trial > 0,
             between(ExpAccmeanbefore, lowacc_ptcpts, highacc_ptcpts),
             CatchAccmean > catch_trial_cutoff,
             Opacity > catch) %>%
      mutate(Acc = ifelse(((Opacity > 0 &
                              ((Key == 'l' & CorrSide == 1) |
                                Key == 'a' & CorrSide == -1)) |
                           (Opacity == 0 & Key == 'nope')), 1, 0),
             block = RoundTo(Trial, 44, ceiling)/44,
             SquareOnset = RoundTo(lilsquareStartTime - flash_circleEndTime, sampling_freq),
             RT = ifelse(Acc == 1, ButtonPressTime - lilsquareEndTime, NA),
             Stim_Sides = as.character(
               ifelse(CorrSide == FlashSide, "Same", "Opposite")),
             CorrSide = if_else(CorrSide == 1, "Right", "Left"),
             Stim_Sides = case_when(sep_vis_fields == "No" ~ Stim_Sides,
                TRUE ~ paste(Stim_Sides, CorrSide, sep = "_"))) %>%
      group_by(block) %>%
      mutate(block_acc = mean(Acc)) %>%
      ungroup() %>%
      mutate_at(vars(Acc, RT),
             funs(ifelse(block_acc < block_acc_cutoff, NA, .))) %>%
      mutate(Trials_below_block_cutoff = sum(block_acc < block_acc_cutoff),
             ExpAccmeanafter = mean(Acc, na.rm = TRUE)) %>%
      filter(between(ExpAccmeanafter, lowacc_ptcpts, highacc_ptcpts)) %>%
      group_by(SquareOnset, Stim_Sides, !!!grouping_constants) %>%
      summarise_at(vars(Acc, RT), funs(mean(., na.rm = TRUE))) %>%
      arrange(Stim_Sides, SquareOnset) %>%
      ungroup() %>%
      mutate_at(vars(Acc, RT), funs(na.approx(., na.rm = FALSE, rule = 2))) #square onsets with an avg of NaN take the mean of the square onset accuracies before and after it ; if it's the last row in the data or after the last NA, take just the most recent NA value
    }
  
  ptcpts_tabulated <- do.call(rbind,lapply(ptcpts, ptcpt_tabulate))

  ptcpts_remaining <- unique(ptcpts_tabulated$participant)
  ptcpts_tabulated <- rbind((ptcpts_tabulated %>%
                               mutate(participant = "All")), ptcpts_tabulated)
  
  time_course_graph_facets <- ptcpts_tabulated %>%
    group_by(SquareOnset, Stim_Sides, !!!grouping_constants) %>%
    summarise(!!dep_var := mean(!!dep_var)) %>%
    ggplot(aes(SquareOnset, !!dep_var, group=Stim_Sides, color=Stim_Sides)) +
      geom_line(alpha = I(2/10), color="grey", show.legend = FALSE) +
      stat_smooth(size=1, span=0.2, se=FALSE, show.legend = FALSE) +
      labs(title = "Accuracy by Square Onset Time, Loess-Smoothed",
           x = "Square Onset (ms)",
           caption = paste("Not fft'ed; 2 data points per x-value per target side per participant\n", as.character(length(ptcpts_remaining)), "participants averaged")) +
        theme(plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5)) +
    facet_wrap(~participant, ncol = 1)
  
  
  #---------#

  
  win_size <- win_size - .001 #because the between function is inclusive
  first_win <- min(ptcpts_tabulated$SquareOnset) + latestart
  last_win <- max(ptcpts_tabulated$SquareOnset) - earlyend
  wins <- seq(first_win, last_win, win_freq)
  hann <- hanning.window(length(wins))
  
  
  amplitude_func <- function(participant_num){
   exp_trials_noncatch <- ptcpts_tabulated %>%
     filter(participant == participant_num)
   
    slidingwindow <- function(win_start){
      win_end <- win_start + win_size
      exp_trials_noncatch %>%
        filter(between(SquareOnset, win_start, win_end)) %>%
        group_by(Stim_Sides, !!!grouping_constants) %>%
        summarise(!!dep_var := mean(!!dep_var)) %>%
        spread(Stim_Sides, !!dep_var)
    }
    as.data.frame(t(mcmapply(slidingwindow, wins))) %>%
      modify(as.numeric) %>%
  	  mutate_at(vars(Same, Opposite), funs(Mod(fft(detrend(.)*hann)))) %>%
      modify(as.numeric) %>%
      mutate(Hz = (row_number()-1)/(n()*win_freq))
  }
  
  amplitudes <- do.call(rbind,mclapply(ptcpts_remaining, amplitude_func)) #may take several minutes to run depending on number of ptcpts
  
  graph_label <- amplitudes %>%
    group_by(!!!grouping_constants) %>%
    summarise() %>%
    ungroup() %>%
    mutate_at(vars(!!!label_constants), funs(paste(quo_name(quo(.)),"=", .))) %>%
    unite(lab, !!!label_constants, sep = "\n", remove = FALSE)
  
  fft_facets <-
    rbind((amplitudes %>% mutate(participant = "All")), amplitudes) %>%
    group_by(participant, Hz) %>%
    summarise_all(mean) %>%
    gather(Flash_and_or_field, Magnitude, -Hz, -c(!!!grouping_constants)) %>%
    ggplot(aes(Hz, Magnitude, color = Flash_and_or_field)) +
    geom_line() +
    scale_x_continuous(breaks = c(4, 8), limits = c(0, 16)) +
    labs(title = paste("FFT of", as.character(dep_var)),
         caption = paste("Data from", as.character(length(ptcpts_remaining)),
                         "participants")) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid.major = element_line(colour="white", size = 3)) +
    facet_wrap(~ factor(participant, levels = c("All", ptcpts_remaining)), ncol = 1) +
    geom_text(data = as.data.frame(graph_label), aes(label = lab, x = 16, y = 70), inherit.aes = FALSE)
  
  grid.arrange(time_course_graph_facets, fft_facets, ncol = 2)
  if (save == "Yes"){
      side_by_side <- arrangeGrob(time_course_graph_facets, fft_facets, ncol = 2)
      ggsave(paste0("data/Time-Series_+_FFT_Plots_", dep_var, "_low_", lowacc_ptcpts,
                    "_high_", highacc_ptcpts, ".pdf"),
                    side_by_side, limitsize = FALSE, width = 25, height = 100)
  }
}

data_graphing(lowacc_ptcpts = .1, highacc_ptcpts = 1.1,
           catch_trial_cutoff = .85, block_acc_cutoff = .0,
           catch = 0, #-1/0 to include/exclude catch trials in analyses
           sep_vis_fields = "No", #Use "Yes" and "No"
           dep_var = "RT", #Use "Acc" or "RT"
           sampling_freq = 1/60,
           win_freq = .001, win_size = .05,
           latestart = 0,  earlyend = 0, #unit = seconds
           ptcpts = 201:230,
           save = "No") #Use "Yes" and "No"
