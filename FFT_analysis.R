library(tidyr)
library(dplyr)
library(ggplot2)
library(DescTools)
library(e1071)
library(pracma)
library(gridExtra)
library(data.table)
library(tables)
library(zoo)
library(tidyverse)

data_graphing <- function(lowacc_ptcpts, highacc_ptcpts, catch_trial_cutoff,
                          block_acc_cutoff, catch, sep_vis_fields, dep_var,
                          win_freq, win_size, latestart, earlyend, ptcpts, save){
  
  ptcpt_tabulate <- function(participant_num){
    dep_var <- as.name(dep_var)
    ptcpt_path <- paste0("data/", participant_num, "/", participant_num, ".csv")
    fread(ptcpt_path, select = c(1:23)) %>%
      filter(Trial > 0) %>%
      mutate(ExpAcc = ifelse(Opacity > 0, Acc, NA),
             CatchAcc = ifelse(Opacity > 0, NA, Acc)) %>%
      filter(between(mean(ExpAcc, na.rm = TRUE),
                     lowacc_ptcpts, highacc_ptcpts),
             mean(CatchAcc, na.rm = TRUE) > catch_trial_cutoff,
             Opacity > catch) %>%
      mutate(Acc = ifelse(((Opacity > 0 &
                              ((Key == 'l' & CorrSide == 1) |
                                Key == 'a' & CorrSide == -1)) |
                           (Opacity == 0 & Key == 'nope')), 1, 0),
             block = RoundTo(Trial, 44, ceiling)/44,
             SquareOnset = RoundTo(lilsquareStartTime - flash_circleEndTime, 1/60),
             RT = ifelse(Acc == 1, ButtonPressTime - lilsquareEndTime, NA),
             Stim_Sides = as.character(
               ifelse(CorrSide == FlashSide, "Same", "Opposite")),
             CorrSide = if_else(CorrSide == 1, "Right", "Left"),
             Stim_Sides = case_when(sep_vis_fields == "No" ~ Stim_Sides,
                TRUE ~  paste(Stim_Sides, CorrSide, sep = "_"))) %>%
      group_by(block) %>%
      mutate(block_acc = mean(Acc),
             !!dep_var := ifelse(block_acc < block_acc_cutoff, NA,
                                 !!dep_var)) %>%
      group_by(participant, SquareOnset, Stim_Sides) %>%
      summarise(!!dep_var := mean(!!dep_var, na.rm = TRUE)) %>%
      arrange(Stim_Sides, SquareOnset) %>%
      ungroup() %>%
      mutate(!!dep_var := na.approx(!!dep_var, na.rm = FALSE),
             !!dep_var := na.locf.default(!!dep_var),
             !!dep_var := na.locf.default(!!dep_var, fromLast = TRUE)) #square onsets with an avg of NaN take the mean of the square onset accuracies before and after it ; if it's the last row in the data or after the last NA, take just the most recent NA value
  }
  
  ptcpts_tabulated <- do.call(rbind,lapply(ptcpts, ptcpt_tabulate))
  ptcpts_remaining <- unique(ptcpts_tabulated$participant)
  ptcpts_tabulated <- rbind((ptcpts_tabulated %>%
                               mutate(participant = 1)), ptcpts_tabulated)
  
  time_course_graph_facets <- ptcpts_tabulated %>%
    group_by(participant, SquareOnset, Stim_Sides) %>%
    summarise(Acc = mean(Acc)) %>%
    ggplot(aes(SquareOnset, Acc, group=Stim_Sides, color=Stim_Sides)) +
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
        group_by(Stim_Sides) %>%
        summarise(Acc = mean(Acc)) %>%
        spread(Stim_Sides, Acc)
    }
    as.data.frame(t(mapply(slidingwindow, wins))) %>%
      modify(as.numeric) %>%
  	  mutate_all(funs(Mod(fft(detrend(.)*hann)))) %>%
      modify(as.numeric) %>%
      mutate(Hz = (row_number()-1)/(n()*win_freq),
             Participant = participant_num)
  }
  
  amplitudes <- do.call(rbind,lapply(ptcpts_remaining, amplitude_func))
  fft_facets <- rbind((amplitudes %>% mutate(Participant = 1)), amplitudes) %>%
    group_by(Participant, Hz) %>%
    summarise_all(mean) %>%
    gather(Flash_and_or_field, Magnitude, -Participant, -Hz) %>%
    ggplot(aes(Hz, Magnitude, color = Flash_and_or_field)) +
    geom_line() +
    scale_x_continuous(breaks = c(4, 8), limits = c(0, 16)) +
    labs(title = paste("FFT of", as.character(dep_var)),
         caption = paste("Data from", as.character(length(ptcpts_remaining)),
                         "participants")) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid.major = element_line(colour="white", size = 3)) +
    facet_wrap(~Participant, ncol = 1)
  
  
  grid.arrange(time_course_graph_facets, fft_facets, ncol = 2)
  if (save == "Yes"){
      side_by_side <- arrangeGrob(time_course_graph_facets,
                                  fft_facets, ncol = 2)
      ggsave("Time-Series_+_FFT_Plotsjbc.pdf", side_by_side,
             limitsize = FALSE, width = 15, height = 8)
  }
}

data_graphing(lowacc_ptcpts = .7, highacc_ptcpts = .9,
           catch_trial_cutoff = .85, block_acc_cutoff = .4,
           catch = 0, #-1/0 to include/exclude catch trials in analyses
           sep_vis_fields = "No", #Use "Yes" and "No"
           dep_var = "Acc",
           win_freq = .001,  win_size = .05,
           latestart = 0,  earlyend = 0, #unit = seconds
           ptcpts = 201:230,
           save = "No") #Use "Yes" and "No"
