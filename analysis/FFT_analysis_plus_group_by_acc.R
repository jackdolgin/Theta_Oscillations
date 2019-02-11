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
                          win_freq, win_size, latestart, earlyend, ptcpts,
                          purpose, table_int, table_skip, output, save){
  
      if (purpose != "ffting"){
          lowacc_ptcpts <- .1
          highacc_ptcpts <- 1.1}
  
      df0 <- data.frame(Participant = numeric(), Acc = numeric())
        
      ptcpt_tabulate <- function(participant_num){
        dep_var <- as.name(dep_var)
        ptcpt_path <- paste0("data/", participant_num, "/", participant_num, ".csv")
        pholder <- fread(ptcpt_path, select = c(1:23)) %>%
          filter(Trial > 0) %>%
          mutate(ExpAcc = ifelse(Opacity > 0, Acc, NA),
                 CatchAcc = ifelse(Opacity > 0, NA, Acc)) %>%
          filter(between(mean(ExpAcc, na.rm = TRUE),
                         lowacc_ptcpts, highacc_ptcpts),
                 mean(CatchAcc, na.rm = TRUE) > catch_trial_cutoff,
                 Opacity > catch)
        df0[nrow(df0) + 1, ] <<- c(participant_num, mean(pholder$ExpAcc))
        pholder %>%
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
          mutate(!!dep_var := na.approx(!!dep_var, na.rm = FALSE, rule = 2)) #square onsets with an avg of NaN take the mean of the square onset accuracies before and after it ; if it's the last row in the data or after the last NA, take just the most recent NA value
      }
      
      ptcpts_tabulated <- do.call(rbind,lapply(ptcpts, ptcpt_tabulate))
      ptcpts_remaining <- unique(ptcpts_tabulated$participant)
      ptcpts_tabulated <- rbind((ptcpts_tabulated %>%
                                   mutate(participant = 1)), ptcpts_tabulated)
      
      
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
      
      if (purpose == "ffting"){
        
          time_course_graph_facets <- ptcpts_tabulated %>%
            ggplot(aes(SquareOnset, Acc, group=Stim_Sides, color=Stim_Sides)) +
              geom_line(alpha = I(2/10), color="grey", show.legend = FALSE) +
              stat_smooth(size=1, span=0.2, method = "loess", se=FALSE, show.legend = FALSE) +
              labs(title = "Accuracy by Square Onset Time, Loess-Smoothed",
                   x = "Square Onset (ms)",
                   caption = paste("Not fft'ed; 2 data points per x-value per target side per participant\n",
                                   as.character(length(ptcpts_remaining)), "participants averaged")) +
                theme(plot.title = element_text(hjust = 0.5),
                      plot.subtitle = element_text(hjust = 0.5)) +
            facet_wrap(~participant, ncol = 1)
          
          
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
      else if (output == "table") {
      
        df <- data.frame(low_acc = numeric(),  high_acc = numeric(), count = numeric(),
             HzSame1 = numeric(),  HzSame2 = numeric(),  HzSame3 = numeric(),
             HzOpp1 = numeric(),   HzOpp2 = numeric(),   HzOpp3 = numeric(),
             MagSame1 = numeric(), MagSame2 = numeric(), MagSame3 = numeric(),
             MagOpp1 = numeric(),  MagOpp2 = numeric(),  MagOpp3 = numeric())
        
        for (low_acc in seq(.3, .9, table_skip)){
          
          df01 <- df0 %>%
            filter(between(Acc, low_acc, low_acc + table_int))
        
          fft_prep <- amplitudes %>%
            filter(Participant %in% df01$Participant) %>%
            group_by(Hz) %>%
            summarise_all(mean) %>%
            filter(Hz < 20)
  
            Same_freqs <- fft_prep %>%
              select(-Opposite) %>%
              arrange(-Same)
            Opp_freqs <- fft_prep %>%
              select(-Same) %>%
              arrange(-Opposite)
          
            df[nrow(df) + 1,] <- c(low_acc, low_acc + table_int, nrow(df01),
                        Same_freqs[1,"Hz"],      Same_freqs[2,"Hz"],      Same_freqs[3,"Hz"],
                        Opp_freqs[1,"Hz"],       Opp_freqs[2,"Hz"],       Opp_freqs[3,"Hz"],
                        Same_freqs[1,"Same"],    Same_freqs[2,"Same"],    Same_freqs[3,"Same"],
                        Opp_freqs[1,"Opposite"], Opp_freqs[2,"Opposite"], Opp_freqs[3,"Opposite"]) 
        }
        
        df
      } 
      
      else {
        
        Same_freqs <- amplitudes %>%
          filter(Hz < 20) %>%
          group_by(Participant) %>%
          filter(Same == max(Same)) %>%
          select(-Opposite, -Same) %>%
          rename(HzSame = Hz)
        Opp_freqs <- amplitudes %>%
          filter(Hz < 20) %>%
          group_by(Participant) %>%
          filter(Opposite == max(Opposite)) %>%
          select(-Opposite, -Same) %>%
          rename(HzOpp = Hz)
        
        tograph <- left_join(Same_freqs, Opp_freqs, by = "Participant") %>%
          left_join(., df0, by = "Participant") 
        bigdot <- geom_point(size = 3)
        thicken <- scale_y_continuous(breaks = c(4, 8), limits = c(0, 12))
        xlabs <- xlab("Mean Accuracy Across Non-Catch Trials")
        ylabs <- ylab("Freq w/ Greatest Magnitude (Hz)")
        themes <- theme(plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              text = element_text(size=20),
              panel.grid.major = element_line(colour="white", size = 3))
        
        same_facets <- tograph %>%
          ggplot(aes(Acc, HzSame, color = Participant)) +
          bigdot + thicken +
          labs(title = "Trial Accuracy vs. Same Side Target Top Frequency",
                 caption = paste("Data from", as.character(nrow(unique(Same_freqs))),
                                 "participants")) + xlabs + ylabs + themes
        
        opp_facets <- tograph %>%
          ggplot(aes(Acc, HzOpp, color = Participant)) +
          bigdot + thicken +
          labs(title = "Trial Accuracy vs. Opposite Side Target Top Frequency",
                 caption = paste("Data from", as.character(nrow(unique(Same_freqs))),
                                 "participants")) + xlabs + ylabs + themes
        
        
        grid.arrange(same_facets, opp_facets, ncol = 1)
        
    }   
}

data_graphing(lowacc_ptcpts = .3, highacc_ptcpts = .6,
           catch_trial_cutoff = .85, block_acc_cutoff = .4,
           catch = 0, #-1/0 to include/exclude catch trials in analyses
           sep_vis_fields = "No", #Use "Yes" and "No"
           dep_var = "Acc",
           win_freq = .001,  win_size = .05,
           latestart = 0,  earlyend = 0, #unit = seconds
           ptcpts = 201:230,
           purpose = "else", #use 'ffting' or anything else
           table_int = .3, table_skip = .3,
           output = "graph", #Use "graph" or "table"
           save = "No") #Use "Yes" and "No"
