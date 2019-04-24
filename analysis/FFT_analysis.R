# Install packages if not already installed, then load them
if (!require(devtools)) install.packages('devtools')
if (!require(smisc)) devtools::install_github("stevenworthington/smisc")
smisc::ipak(c("utils", "tidyr", "dplyr", "ggplot2", "DescTools", "e1071",
              "pracma", "gridExtra", "data.table", "tables", "zoo", "parallel",
              "scales", "purrr", "lazyeval", "stats", "gdata", "viridis"))

# Main function begins (encompasses other functions)
data_graphing <- function(catch_cutoff, block_floor, mini_block_floor,
                          mini_block_ceil, side_bias, catch, pre_floor, pre_ceil,
                          post_floor, post_ceil, sep_vis_fields, smooth_method,
                          dep_var, sbtr, samp_freq, latestart, earlyend, pcpts,
                          output, save_output, clumps, xaxisvals, shuff, pval){

  grouping_cnsts <- quos(participant, Trials_filtered_out, Acc_prefilter,       # Columns that are frequently used for grouping, variable means
                             Acc_postfilter, CatchAcc)                          # don't have to type them out every time we use them for grouping
  
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
             grepl("fully alert", Q9),
             Opacity > catch,                                                   #             catch trials if 'catch' parameter is assigned to 0; if it's assigned to -1, this line does nothing)
             Side_Diff < side_bias) %>%
      mutate(Acc_prefilter = mean(Acc, na.rm = TRUE),                           # Creates column indicating mean accuracy before we've filtered for 'block_floor', unlike 'Acc_postfilter'
             CTI = RoundTo(RoundTo(lilsquareStartTime - flash_circleEndTime,
                                   1 / 60), samp_freq)) %>% 
      filter(between(Acc_prefilter, pre_floor, pre_ceil, incbounds = TRUE),                       # Filters out participants whose non-catch, pre-block-filtering accuracy is outside of desired range
             between(CTI, min(CTI) + latestart, max(CTI) - earlyend)) %>%
      mutate(block = RoundTo(Trial, 54, ceiling) / 54,                          # Creates column indicating trial's block
             RT = ifelse(Acc == 1 & ButtonPressTime - lilsquareStartTime > .1,  #                indicating RT after target appeared on screen, only for correct trials with an RT < 100 ms
                         ButtonPressTime - lilsquareStartTime, NA),
             Stim_Sides = as.character(
               ifelse(CorrSide == FlashSide, "Valid", "Invalid")),              #                indicating whether cue was valid or invalid
             CorrSide = case_when(CorrSide == 1 ~ "Right",                      #                indicating which side the target appeared on
                                  CorrSide == -1 ~ "Left", TRUE ~ "Bottom"),
             Stim_Sides = case_when(sep_vis_fields == "No" ~ Stim_Sides,        # Overwrites 'Stim_Sides' column if 'sep_vis_fidels' parameter == 'Yes' to include which side of screen target was on,
                TRUE ~ paste(CorrSide, Stim_Sides, sep = "_"))) %>%             # as well as whether it was valid with cue; if 'sep_vis_fields' parameter == 'No', leaves 'StimSides' unchanged
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
  hann <- hanning.window(length(unique(cmbd$CTI)))                       # Creates Hanning window, gets applied later
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
  amplitude <- function(x, ...){ x %>%
      group_by_(.dots = lazy_dots(...)) %>%
      mutate_at(vars(locations), funs(Mod(fft(detrend(.) * hann)))) %>%         # Detrends, multiplies by Hanning window, applies FFT, and then takes magnitude   
      mutate(Hz = (row_number() - 1) / (n() * samp_freq)) %>%
      ungroup() %>%
      select(-CTI)
  }
  amps <- amplitude(cmbd_w, participant) 

  
 # Set Up Graphing
 
  graph <- function(y, x) {
    y(x) +
      theme_bw() +
      scale_color_viridis_d(option = "C",
                            end = .7 + RoundTo(.0001 * RoundTo(length(locations),
                                                               4, floor), .2, ceiling),
                            labels = sapply(locations, function(x) gsub("_", " ", x),
                                            USE.NAMES = FALSE, simplify = TRUE)) +
      guides(colour = guide_legend(reverse = TRUE), fill = "none")
  }
  
  avgd <- function(y, x) {
    graph(y, x) + geom_line(size = 1.5) +
      theme(legend.key.size = unit(.55, "in")) +
      labs(subtitle = paste("Data from", as.character(length(pcpts)),
                       "participants"))
  }
  
  save_avgd <- function(x) { x %>%
      ggsave(filename = file.path("analysis", paste0(output, ".pdf")),
             width = xaxisvals)}
  
  idvl <- function(y, x){
    graph(y, x) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)) +
      facet_wrap(~factor(participant, levels = pcpts),
                 ncol = 4, scales = 'free_x')                                   # 'free' means the y_axis isn't fixed from participant to participant
  }
  
  t_srs <- function(x){ 
    cmbd_w %>%
      gather(Location, !!dep_var, -c(CTI, !!!grouping_cnsts)) %>%
      right_join(gather(conf_int(cmbd_w, CTI), Location, Conf_Int, -CTI),
                 by = c("CTI", "Location")) %>%
      group_by(CTI, Location, Conf_Int, !!!head(grouping_cnsts, x)) %>%
      summarise(!!dep_var := mean(!!dep_var)) %>%                               # Keeps either RT or Acc column—depending on whether 'dep_var' parameter == 'RT' or 'Acc'
      ggplot(aes(CTI, !!dep_var, group = Location, color = Location,
                 fill = Location, ymin = !!dep_var - Conf_Int,
                 ymax = !!dep_var + Conf_Int)) + 
      labs(title = paste(dep_var, "by Cue-Target Interval"),
           x = "Cue-Target Interval (ms)")
  }
  
  fft_g <- function(x) { x +
      labs(title = paste("FFT of Target", dep_var),
           col = "Target Location") +
      theme(panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank()) +
      scale_x_continuous(name = "Frequency (Hz)", limits = c(0, xaxisvals), 
         breaks = seq(0, xaxisvals, 1 / (length(unique(amps$Hz)) * samp_freq)))
    }
  

  # Conditionally returns data frame of prelim participant or post-FFT data, exits function
  if (output == "prelim_table") {
    
    # Output Prelim Data Table -------------------------------------------------

    if (save_output == "Yes") {                                             # Conditionally saves data frame as .csv in analysis folder
      write.csv(cmbd, file.path("analysis", "Prelim_table.csv"))
    }
    return(cmbd_w)
  } else if (output == "fft_table") {
    # Output FFT Data Table ----------------------------------------------------
      
    if (save_output == "Yes") {                                                 # Conditionally saves data frame as .csv in analysis folder
      write.csv(amps, file.path("analysis", "FFT_table.csv"))
    }
    return(amps)
  } else if (output == "individuals"){
    
    # Output Individuals' Plots ------------------------------------------------
    
    # Produces left half of final graph
    ts_facets <- idvl(t_srs, 1) +
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
    fft_facets <- idvl(fft_g, amps %>%
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
    
    (avgd(t_srs, 0) +
      theme(panel.grid = element_blank()) +
      geom_ribbon(alpha = 0.15, aes(color = NULL))) %>%
      save_avgd
    
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
    amps <- shuffle(.data = cmbd_w, n = shuff,
                    perm_cols = c("participant", locations)) %>%
      mutate(randsamp = RoundTo(row_number(), n() / shuff, ceiling)) %>%
      amplitude(participant, randsamp) %>%
      group_by(Hz, randsamp) %>%
      summarise_at(vars(locations), mean) %>%
      group_by(Hz) %>%
      summarise_at(vars(locations), funs(quantile(., probs = 1 - pval))) %>%
      combine(amps %>% group_by(Hz) %>% summarise_at(vars(locations), mean),
              names = (c("Significance Cutoff", "Observed Data"))) %>%
      gather(Location, Magnitude, -c(Hz, source)) %>%
      right_join(gather(conf_int(amps, Hz), Location, Conf_Int, -Hz),
                 by = c("Hz", "Location"))
    (avgd(fft_g, ggplot(amps, aes(Hz, Magnitude, col = Location, linetype = source, 
                 ymin = Magnitude - Conf_Int, ymax = Magnitude + Conf_Int))) +
      scale_linetype_manual(values = c("dashed", "solid")) +
      labs(linetype = "",
           caption = paste("Significance threshold at p <",
                           as.character(pval))) +
      geom_errorbar(data = filter(amps, source == "Observed Data"),
                        position = position_dodge(width = 0.9), width = 0.2) +
      geom_point(size = 3, data = amps %>% spread(source, Magnitude) %>%
                   filter(`Observed Data` > `Significance Cutoff`) %>%
                   select(-c(`Significance Cutoff`, Conf_Int)) %>%
                   gather(source, Magnitude, -Hz, -Location), 
                 aes(ymin = NULL, ymax = NULL))) %>%
      save_avgd

  }
}



# Sets inputs for the 'data_graphing' function
data_graphing(catch_cutoff = .85,                                               # Filters out participants with a catch trial accuracy below this value 
              block_floor = .25,                                                # Converts a trial's accuracy to NA if its block's accuracy below this value
              mini_block_floor = .35,                                           # Converts a trial's accuracy to NA if its mini-block's accuracy below this value
              mini_block_ceil = .9,                                             # Converts a trial's accuracy to NA if its mini-block's accuracy above this value
              side_bias = .3,
              catch = 0,                                                        # -1/0 to include/exclude catch trials in analyses
              pre_floor = .35,                                                  # Filters participants by their accuracy before their blocks below 'block_floor' have been interpolated over
              pre_ceil = .75,                                                   # Line above is the floor of the filter, this line is the ceiling
              post_floor = .35,                                                 # Filters participants by their accuracy after their blocks below 'block_floor' have been interpolated over 
              post_ceil = .75,                                                  # Line above is the floor of the filter, this line is the ceiling
              sep_vis_fields = "No",                                            # Use 'No' and 'Yes'; 'No' means invalid cue and target trials are analyzed differently at each CTI than trials with valid cue and target; 'Yes' means we further divide analyses by which side of screen target appeared, yielding four conditions (valid/invalid and left/right)
              smooth_method = "loess",                                          # See geom_smooth documentation for available smoothing methods
              dep_var = "Acc",                                                  # Use 'Acc' or 'RT'
              sbtr = "No",
              samp_freq = 1 / 60,                                               # Equals spacing between CTI intevals (in seconds)
              latestart = 0, earlyend = 0,                                      # Filters out, for FFT analysis, trials with a CTI < 'latestart' or a CTI > ((largest CTI [so 1.3]) - 'earlyend'); units = seconds
              pcpts = 401:424,#304:322,
              output = "combined_fft",                                           # Use "combined_fft", "combined_ts", "individuals", "prelim_table" (before interpolation and FFT'ing), and "fft_table"
              save_output = "Yes",                                              # Use "Yes" and "No" (if 'no', table outputs will still be visible in R)
              clumps = 3,                                                       # Use 1 (no clumping) and 3 (each bin is the average of itself and its two neighbors)
              xaxisvals = 10,
              shuff = 5000,
              pval = .001)
