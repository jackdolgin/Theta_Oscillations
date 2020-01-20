batch <- "2a"

dvcol <- (ifelse(dep_var == "Accuracy", "Acc", "RT")) %>% as.name             # Converts character to name/symbol so we can refer to it as a column using tidycom

dv_mutate <- function(x, s){ mutate_at(x, vars(!!dvcol), ~!!s)}

blocksize <-
  if (batch == "3a"){
    79 } else if (batch %in% c("1a", "1b")){
      54 } else if (batch %in% c("2a", "2b", "2c")){
        80}

dem_df <- fread(here("analysis", "data", "demographics.csv"))

cmbd2a <- 
  dir(path = here("analysis", "data", batch),
      pattern = "*.csv",
      full.names = TRUE,
      recursive = TRUE) %>%
  map_df(fread) %>%
  filter(Trial > 0) %>%                                                       # Prunes practice trials
  mutate(
    Stim_Sides = ifelse(CorrSide == FlashSide, "Valid", "Invalid"),           # Creates column indicating whether cue was valid or invalid
    CTI = RoundTo(lilsquareStartTime - flash_circleEndTime, samp_per),
    block = RoundTo(Trial, blocksize, ceiling) / blocksize,                   # Creates column indicating trial's block
    RT = ifelse(Acc == 1, ButtonPressTime - lilsquareStartTime, NA)) %>%      # RT only gets calculated for trials with correct responses
  dv_mutate(expr(ifelse(                                                      # DV becomes NA if...
    (is.na(RT) | RT > .1) &                                                   #                 ... RT <= 100 ms
      block %in% blocks_desired &                                             #                 ... not in a desired block
      lilsquareEndFrame - lilsquareStartFrame == 2 &                          #                 ... lil squares on screen for longer than 2 frames; affects both these trials and...
      lag(lilsquareEndFrame) - lag(lilsquareStartFrame) == 2,                 # ... subsequent trials since lilsquare never gets removed from screen until end of next trial
    .,
    NA))) %>%
  group_by(participant) %>%
  mutate(
    wm_Acc = ifelse(session == "4", mean(Acc_wmarith, na.rm = TRUE), 1),
    CatchAcc = ifelse(Opacity != 0, NA, Acc) %>% mean(na.rm = TRUE)) %>%
  filter(
    CatchAcc >= catch_floor,                                           # Prunes participants whose catch accuracy is below desired threshold
    Opacity != 0) %>%                                                     #        catch trials
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
    Acc_prefilter %>% between(min(pre_range), max(pre_range), incbounds = TRUE), # Prunes participants whose non-catch, pre-block-filtering accuracy is outside of desired range
    between(CTI, min(CTI_range), max(CTI_range))) %>%                  #        trials outside of desired CTI range
  group_by(block, participant) %>%
  mutate(block_acc = mean(Acc, na.rm = TRUE)) %>%                                         # Creates column indicating block's mean accuracy
  group_by(participant) %>%
  mutate(miniblock = RoundTo(row_number(), 16, ceiling) %>% `/` (16)) %>%                   #                           trial's mini-block (every 16 non-catch trials the opacity was readjusted)
  group_by(participant, miniblock) %>%
  mutate(miniblock_avg = mean(Acc, na.rm = TRUE)) %>%
  ungroup() %>%
  dv_mutate(expr(ifelse(                                    # Changes `Acc` and `RT` (dv) column values to NA if...
    between(block_acc, min(block_range), max(block_range)) & between(           #  ... trial's block accuracy outside of `block_range`...
      miniblock_avg, min(miniblock_range), max(miniblock_range)),
    .,
    NA))) %>%  # ... or miniblock was not in the desired range or block was not in `blocks_desired`
  group_by(participant) %>%
  mutate(
    Trials_filtered_out = sum(is.na(Acc)) / n(),
    Acc_postfilter = mean(Acc, na.rm = TRUE)) %>%                      
  filter(
    between(Acc_postfilter, min(post_range), max(post_range)),             # Prunes participants whose non-catch, post-block-filtering accuracy is outside of desired range
    Trials_filtered_out <= filtered_cap | averaging == "Aggregate")






batch <- "3a"

dvcol <- (ifelse(dep_var == "Accuracy", "Acc", "RT")) %>% as.name             # Converts character to name/symbol so we can refer to it as a column using tidycom

dv_mutate <- function(x, s){ mutate_at(x, vars(!!dvcol), ~!!s)}

blocksize <-
  if (batch == "3a"){
    79 } else if (batch %in% c("1a", "1b")){
      54 } else if (batch %in% c("2a", "2b", "2c")){
        80}

dem_df <- fread(here("analysis", "data", "demographics.csv"))

cmbd3a <- 
  dir(path = here("analysis", "data", batch),
      pattern = "*.csv",
      full.names = TRUE,
      recursive = TRUE) %>%
  map_df(fread) %>%
  filter(Trial > 0) %>%                                                       # Prunes practice trials
  mutate(
    Stim_Sides = ifelse(CorrSide == FlashSide, "Valid", "Invalid"),           # Creates column indicating whether cue was valid or invalid
    CTI = RoundTo(lilsquareStartTime - flash_circleEndTime, samp_per),
    block = RoundTo(Trial, blocksize, ceiling) / blocksize,                   # Creates column indicating trial's block
    RT = ifelse(Acc == 1, ButtonPressTime - lilsquareStartTime, NA)) %>%      # RT only gets calculated for trials with correct responses
  dv_mutate(expr(ifelse(                                                      # DV becomes NA if...
    (is.na(RT) | RT > .1) &                                                   #                 ... RT <= 100 ms
      block %in% blocks_desired &                                             #                 ... not in a desired block
      lilsquareEndFrame - lilsquareStartFrame == 2 &                          #                 ... lil squares on screen for longer than 2 frames; affects both these trials and...
      lag(lilsquareEndFrame) - lag(lilsquareStartFrame) == 2,                 # ... subsequent trials since lilsquare never gets removed from screen until end of next trial
    .,
    NA))) %>%
  group_by(participant) %>%
  mutate(
    wm_Acc = ifelse(session == "4", mean(Acc_wmarith, na.rm = TRUE), 1),
    CatchAcc = ifelse(Opacity != 0, NA, Acc) %>% mean(na.rm = TRUE)) %>%
  filter(
    CatchAcc >= catch_floor,                                           # Prunes participants whose catch accuracy is below desired threshold
    Opacity != 0) %>%                                                     #        catch trials
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
    Acc_prefilter %>% between(min(pre_range), max(pre_range), incbounds = TRUE), # Prunes participants whose non-catch, pre-block-filtering accuracy is outside of desired range
    between(CTI, min(CTI_range), max(CTI_range))) %>%                  #        trials outside of desired CTI range
  group_by(block, participant) %>%
  mutate(block_acc = mean(Acc, na.rm = TRUE)) %>%                                         # Creates column indicating block's mean accuracy
  group_by(participant) %>%
  mutate(miniblock = RoundTo(row_number(), 16, ceiling) %>% `/` (16)) %>%                   #                           trial's mini-block (every 16 non-catch trials the opacity was readjusted)
  group_by(participant, miniblock) %>%
  mutate(miniblock_avg = mean(Acc, na.rm = TRUE)) %>%
  ungroup() %>%
  dv_mutate(expr(ifelse(                                    # Changes `Acc` and `RT` (dv) column values to NA if...
    between(block_acc, min(block_range), max(block_range)) & between(           #  ... trial's block accuracy outside of `block_range`...
      miniblock_avg, min(miniblock_range), max(miniblock_range)),
    .,
    NA))) %>%  # ... or miniblock was not in the desired range or block was not in `blocks_desired`
  group_by(participant) %>%
  mutate(
    Trials_filtered_out = sum(is.na(Acc)) / n(),
    Acc_postfilter = mean(Acc, na.rm = TRUE)) %>%                      
  filter(
    between(Acc_postfilter, min(post_range), max(post_range)),             # Prunes participants whose non-catch, post-block-filtering accuracy is outside of desired range
    Trials_filtered_out <= filtered_cap | averaging == "Aggregate")













combine(cmbd2a, cmbd3a) %>%
  filter(!is.na(Acc)) %>%
  group_by(source) %>%
  mutate(source_count = n()) %>%
  group_by(source, source_count, miniblock_avg) %>%
  summarise(miniblock_count = 100 * n()/ mean(source_count)) %>%
  ggplot(aes(x = miniblock_avg, y = miniblock_count, fill = source)) +
  geom_col(color="#e9ecef", alpha = .5, position = 'identity') +
  # geom_histogram(color="#e9ecef", alpha = .5, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080"))


cmbd2a %>%
  filter(!is.na(Acc)) %>%
  pull(miniblock_avg) %>%
  sd

  
cmbd3a %>%
  filter(!is.na(Acc)) %>%
  pull(miniblock_avg) %>%
  sd



cmbd2a %>%
  filter(!is.na(Acc)) %>%
  group_by(CTI) %>%
  count %>%
  pull(n) %>%
  median


cmbd3a %>%
  filter(!is.na(Acc)) %>%
  group_by(CTI) %>%
  count %>%
  pull(n) %>%
  median
