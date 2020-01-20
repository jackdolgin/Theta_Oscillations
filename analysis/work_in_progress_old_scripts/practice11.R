# Install packages if not already installed, then load them
if (!require(devtools)) install.packages("pacman")
pacman::p_load(
  utils, tidyr, dplyr, ggplot2, DescTools, bspec, pracma, gridExtra, data.table,
  zoo, scales, stats, viridis, gginnards, purrr, furrr, stingr, rlang, here, binom)

# Main function begins (encompasses other functions)
main_function <- function(averaging, visualize, batch, wm_exp, iso_sides, sbtr, samp_per,
                          clumps, dep_var, α, shuff, mult_correcs, trends,
                          smooth_method, show_CI, display, win_func, xmax, duration, attn_filter,
                          catch_floor, side_bias, wm_floor, invalid_floor,
                          pre_range, post_range, filtered_cap, block_range,
                          blocks_desired, miniblock_range, CTI_range){
  
  dvcol <- (ifelse(dep_var == "Accuracy", "Acc", "RT")) %>% as.name             # Converts character to name/symbol so we can refer to it as a column using tidycom
  
  dv_mutate <- function(x, s){ mutate_at(x, vars(!!dvcol), ~!!s)}
  
  if (batch == "3a"){
    blocksize <- 79
  }
  
  dem_df <- fread(here("analysis", "data", "demographics.csv"))
  
  cmbd <- dir(path = here("analysis", "data", batch),
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
    group_by(participant) %>%
    mutate(
      wm_Acc =
        ifelse(session == "4",
               mean(Acc_wmarith, na.rm = TRUE),
               1),
      CatchAcc =
        ifelse(Opacity != 0, NA, Acc) %>% mean(na.rm = TRUE)) %>%
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
    group_by(CorrSide, participant) %>%
    mutate(Side_Acc = mean(Acc, na.rm = TRUE)) %>%
    group_by(participant) %>%
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
    group_by(block, participant) %>%
    mutate(block_acc = mean(Acc)) %>%                                         # Creates column indicating block's mean accuracy
    group_by(participant) %>%
    mutate(
      rown = row_number(),
      miniblock = RoundTo(rown, 16, ceiling) %>% `/` (16)) %>%                   #                           trial's mini-block (every 16 non-catch trials the opacity was readjusted)
    group_by(participant, miniblock) %>%
    mutate(miniblock_avg = mean(Acc)) %>%
    ungroup() %>%
    dv_mutate(expr(ifelse(                                    # Changes `Acc` and `RT` (dv) column values to NA if...
      between(block_acc, block_range[1], block_range[2]) & between(           #  ... trial's block accuracy outside of `block_range`...
        miniblock_avg, miniblock_range[1], miniblock_range[2]),
      .,
      NA))) %>%  # ... or miniblock was not in the desired range or block was not in `blocks_desired`
    group_by(participant) %>%
    mutate(
      Trials_filtered_out = sum(is.na(Acc)) / n(),
      Acc_postfilter = mean(Acc, na.rm = TRUE)) %>%                      
    filter(
      between(Acc_postfilter, post_range[1], post_range[2]),             # Prunes participants whose non-catch, post-block-filtering accuracy is outside of desired range
      Trials_filtered_out <= filtered_cap)
  
  pcpts_original <- unique(cmbd$participant)
  
  if(averaging == "Individual") {
    grouping_cnsts <- quos(participant, Trials_filtered_out, Acc_prefilter,     
                           Acc_postfilter, CatchAcc)
  } else if(averaging == "Aggregate") {
    grouping_cnsts <- quos(participant)
    cmbd$participant <- 1
  }
  
  pcpts <- unique(cmbd$participant)
  
  cmbd_summarize <- function(to_shuf){
    cmbd %>%
      group_by(Stim_Sides, !!!grouping_cnsts) %>%
      dv_mutate(expr(case_when(!!to_shuf ~ sample(.), TRUE ~ .))) %>%
      group_by(CTI, Stim_Sides, !!!grouping_cnsts) %>%
      summarise(
        Conf_Int = ifelse(
          !!to_shuf | (averaging == "Individual" & visualize == "Individual"),
          0,
          binom.exact(sum(!!dvcol, na.rm = TRUE), sum(!is.na(!!dvcol))) %>%
            summarise((upper - lower) / 2) %>%
            pull),
        !!dvcol := mean(!!dvcol, na.rm = TRUE)) %>%
      dv_mutate(expr(na.approx(., na.rm = FALSE, rule = 2))) %>%
      dv_mutate(expr(rollapply(., clumps + 1, mean, partial = TRUE))) %>%
      group_by(Stim_Sides, !!!grouping_cnsts) %>%
      dv_mutate(expr(na.approx(., na.rm = FALSE, rule = 2))) %>%
      ungroup() %>%
      mutate_at(vars(Stim_Sides), ~case_when(!sbtr ~ ., TRUE ~ str_replace(
        ., regex("_?(In)?valid", ignore_case = T), "Invalid - Valid"))) %>%
      group_by(CTI, Stim_Sides, !!!grouping_cnsts, Conf_Int) %>%
      summarise_at(vars(!!dvcol), ~ifelse(sbtr, first(.) - last(.), .)) %>%
      group_by(Stim_Sides, !!!grouping_cnsts) %>%
      dv_mutate(expr(case_when("Detrending" %in% trends ~
                                 . - (polyfit(CTI, ., 2) %>% polyval(CTI)),
                               TRUE ~ .))) %>%
      dv_mutate(expr(case_when("Demeaning" %in% trends ~                    # Works just like the detrending except for demeaning
                                 . - mean(.), TRUE ~ .)))
  }
  
  grouping_viz <- if (visualize == "Individual"){
    grouping_cnsts} else if (visualize == "Aggregate") {
      quo()
    }
  
  
  line_size <- ifelse(averaging == "Individual" & visualize == "Individual",
                      .5, 1.5)
  
  # TS Graph
  time_srs_graph <- cmbd_summarize(FALSE) %>%
    group_by(Stim_Sides, CTI, !!!grouping_viz) %>%
    summarise(
      Conf_Int = case_when(
        averaging == "Individual" & visualize == "Aggregate" ~
          qnorm(.975) * std_err(!!dvcol),
        averaging == "Aggregate"  ~ mean(Conf_Int),
        TRUE ~ 0),
      !!dvcol := mean(!!dvcol)) %>%
    ggplot(aes(CTI, !!dvcol, group = Stim_Sides, color = Stim_Sides,
               fill = Stim_Sides, ymin = !!dvcol - Conf_Int,
               ymax = !!dvcol + Conf_Int)) +
    labs(title = paste(dvcol, "by Cue-Target Interval, Exp.", batch),
         x = "Cue-Target Interval (ms)") +
    geom_ribbon(alpha = ifelse(show_CI, 0.15, 0), aes(color = NULL)) +
    geom_line(alpha = I(2 / 10), color = "grey", show.legend = FALSE, size = line_size) +
    stat_smooth(
      method = tolower(smooth_method), span = 0.2, se = FALSE, size = line_size,
      show.legend = ifelse("time_series" %in% display & "FFT" %in% display, FALSE, TRUE))
  

  
  win <- if (win_func == "Tukey"){
    tukeywindow(n_distinct(cmbd$CTI), .5)} else {
      match.fun(paste0(tolower(win_func), "window")) (n_distinct(cmbd$CTI))}
  
  locations <- unique(cmbd$Stim_Sides)                                         
  
  empty_rows_per_pcpt <- cmbd$CTI %>% range %>% diff %>% - duration %>%         # Add empty rows (other than subject ID) as additional CTI's needed to reach desired padded...
    `/` (-samp_per) %>% + .00001 %>% floor %>% - 1                            # ... `duration` of intervals for each participant for each shuffle (`y` represents each shuffle)...
  
  amplitude <- function(x, to_shuf){
    set.seed(x)
    cmbd_summarize(to_shuf) %>%
      dv_mutate(expr(. * win / Norm(win))) %>%               # Apply window
      ungroup() %>%
      add_row(
        participant = rep(pcpts, empty_rows_per_pcpt * (length(locations))),
        Stim_Sides = rep(locations, each = empty_rows_per_pcpt * length(pcpts)),
        !!dvcol := 0) %>%
      group_by(Stim_Sides, participant) %>%
      mutate(
        Observed = fft(!!dvcol) %>% `*`(sqrt(2 / n())) %>% Mod %>% `^` (2),
        Hz = (rank(row_number()) - 1) / (n() * samp_per)) %>%
      filter(
        rank(Hz) - 1 <= floor(n_distinct(Hz) / 2),                  
        Hz < xmax) %>%
      group_by(Hz, Stim_Sides, !!!grouping_viz) %>%
      summarise(
        Conf_Int = ifelse(averaging == "Individual" &
                            visualize == "Aggregate" & to_shuf == FALSE,
                          qnorm(.975) * std_err(Observed),
                          0),
        Observed = mean(Observed))
  }
  
  tasco <- 1:shuff %>%
    future_map_dfr(~amplitude(.x, TRUE)) %>%
    group_by(Hz, Stim_Sides) %>%
    mutate(ntile = PercentRank(desc(Observed)))
  
  fillerfunc <- function(Hz0, Stim_Sides0, Observed0, xray0, yray0, zray0){
    tasco %>%                                                        # ... finds row in `amps_shuff` that...
      filter(
        Hz0 == Hz,
        Stim_Sides0 == Stim_Sides,
        !!parse_expr(yray0)) %>%                                            # ....3) with a smaller Observed0 than `amps`' Observed...
      top_n(1, Observed) %>%                                               # ... and 4) with largest Observed0 among the remaining rows (we are being conservative and ntile reference...
      pull(xray0)  
  }
  
  amped_up <- future_map_dfr(1, ~amplitude(.x, FALSE)) %>%
    mutate(ntile2 = future_pmap_dbl(
      list(Hz0 = Hz, Stim_Sides0 = Stim_Sides, Observed0 = Observed, xray0 = "ntile",
           yray0 = "Observed0 > Observed | Observed == min(Observed)"), fillerfunc)) %>%
    group_by(Stim_Sides) %>%
    mutate_at(vars(ntile2), list(~ . / p.adjust(., method = mult_correcs))) %>%
    ungroup() %>%
    mutate(`Significance Cutoff` = future_pmap_dbl(
      list(Hz0 = Hz, Stim_Sides0 = Stim_Sides, Observed0 = Observed, xray0 = "Observed",
           yray0 = "ntile >=  zray0 * α", zray0 = ntile2), fillerfunc))
  
  if(averaging == "Individual" & visualize == "Individual"){
    
    amped_for_graph <- amped_up %>%
      pivot_longer(c(Observed, `Significance Cutoff`),
                   "source", values_to = "Power") %>%
      filter(source == "Observed")
    
    line_col <- NULL
    sig_points <- NULL
    
  } else{
    amped_for_graph <- amped_up %>%
      pivot_longer(c(Observed, `Significance Cutoff`),
                   "source", values_to = "Power") %>%
      mutate_at(vars(Conf_Int), ~ifelse(source == "Observed", ., 0))
    
    line_col <- as.symbol("source")
    sig_points <- list(geom_point(
      size = 3,
      data = amped_up %>%
        filter(Observed > `Significance Cutoff`) %>%
        select(-c(`Significance Cutoff`, ntile2)) %>%
        pivot_longer(Observed, "source", values_to = "Power"),
      aes(x = Hz,
          y = Power,
          linetype = NULL,
          fill = NULL,
          ymin = NULL,
          ymax = NULL)))
  }
  
  fft_aggregate_graph <- ggplot(amped_for_graph,
                                aes(Hz,
                                    Power,
                                    col = Stim_Sides,
                                    linetype = !!line_col,
                                    fill = Stim_Sides,
                                    ymin = Power - Conf_Int,
                                    ymax = Power + Conf_Int)) +
    labs(title = paste0("FFT of Target ", dep_var, ", Exp. ", batch),
         col = "Target Location", y = "Spectral Power") +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank())
  
  fft_x <- round(1 / duration, 1)
  xaxis_r <- RoundTo(xmax, fft_x)
  
  fft_aggregate_graph2 <- list(
    geom_line(size = line_size),
    scale_linetype_manual(values = c("solid", "dashed")),
    scale_x_continuous(name = "Frequency (Hz)",
                       limits = c(0, xmax),
                       breaks = seq(0, xmax, ifelse(fft_x > .5,
                                                    round(fft_x, 2), 1))),
    labs(linetype = "",
         caption = paste("Significance threshold at p < ",
                         as.character(α))),
    geom_ribbon(alpha = ifelse(show_CI, 0.15, 0), aes(color = NULL)),
    sig_points)
  
  
  
  
  clustering_theme <-
    if (averaging == "Individual" & visualize == "Individual"){ list(
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)),
      facet_wrap(~factor(participant, levels = pcpts),
                 ncol = round(length(pcpts) / 3), scales = "free_x")
      )} else { list(
        theme(legend.key.size = unit(.55, "in")),
        labs(subtitle = paste("Data from", as.character(length(pcpts_original)),
                            "participants")))
    }
  
  viridis_cols <- .7 +
    RoundTo(.0001 * RoundTo(length(locations), 4, floor), .2, ceiling)
  
  common_theme <- list(
    theme_bw(),
    scale_color_viridis_d(option = "C",
                          end = viridis_cols,
                          labels = str_replace(sort(locations), "_", " ")),
    scale_fill_viridis_d(option = "C",
                         end = viridis_cols),
    guides(fill = FALSE),
    clustering_theme
  )
  
  

  time_srs_graph +
    common_theme


  fft_aggregate_graph + 
    common_theme +
    fft_aggregate_graph2
}



# Sets inputs for the `main_function` function
main_function(averaging = "Aggregate", #individual or grouped
              
              visualize = "Aggregate",
              
              batch = "3a",                                                      # Either `1a`, `1b`, `2a`, `2b`, `2c`, or `3a`
              
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
              
              show_CI = TRUE,
              
              display = c("time_series", "FFT"),
              
              win_func = "Tukey",                                               # Choose between the following types of windowing functions: `Tukey`, `Square`, `Hann`, `Welch`, ...
              # ... `Triangle`, `Hamming`, `Cosine`, or `Kaiser`
              
              xmax = 11,                                                        # Greatest x-axis value included in graph
              
              duration = 1,                                                     # Duration (Seconds) Analyzed Including Padding
              
              attn_filter = FALSE,                                              # Either `FALSE` or `TRUE`, which prunes any participants who indicated in the post-task questionnaire...
              # ... that they either dozed off at one point or were not fully focused on at least two of the eight...
              # ... blocks (the other options were that they were fully alert on all blocks or fully alert on all but...
              # ... one block)
              
              catch_floor = .85,                                                # Prunes participants who perform below this hit rate on catch trials, which featured no target
              
              side_bias = .5,                                                   #                     whose hit rate at one visual field - another visual field is > `side_bias`
              
              wm_floor = .7,                                                    #                           working memory task accuracy is < `wm_floor`
              
              invalid_floor = .02,                                              #                           invalid trial accuracy is < `invalid_floor`
              
              pre_range = c(.25, .85),                                          #                           unfiltered/scrutinized data is outside of the selected range
              
              post_range = c(.25, .85),                                         #                           filtered data is outside of the selected range
              
              filtered_cap = .5,                                                #                     with at least `filtered_cap` % of trials filtered out
              
              block_range = c(.20, .80),                                        # Interpolates over trials if the average hit rate in that block, every 48 trials, is outside of the...
              # ... select range
              
              blocks_desired = 1:8,                                             #                                 block is not in the `blocks_desired`
              
              miniblock_range = c(.20, .80),                                    #                                 average hit rate in that mini-block, every 16 trials which is how...
              # ... often the task difficulty was adjusted to titrate to 65%, is outside this range
              
              CTI_range = c(.3, 1.09)                                           # Filters trials outside this range of CTI bins
)
