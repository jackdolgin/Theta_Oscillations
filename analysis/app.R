# Install packages if not already installed, then load them
# if (!require(devtools)) install.packages("pacman")
# pacman::p_load(utils, tidyr, dplyr, ggplot2, DescTools, bspec, pracma,
#                gridExtra, data.table, tables, zoo, scales, stats, gdata,
#                viridis, gginnards, purrr, shiny, shinyWidgets)
# pacman::p_load_gh("moodymudskipper/safejoin")

library(utils); library(tidyr); library(dplyr); library(ggplot2)
library(DescTools); library(bspec); library(pracma); library(gridExtra)
library(data.table); library(tables); library(zoo); library(scales)
library(stats); library(gdata); library(viridis); library(gginnards);
library(purrr); library(shiny); library(shinyWidgets); library(safejoin)

ui <- fluidPage(
  title = "Theta Oscillations Analysis",
  plotOutput("mygraph",  height = 800, width = 1350),
  hr(),       
  fluidRow(class = "text-center",
           column(4, h3( "Data Set, Task, and Graph Choice"), offset = 3)), br(), br(),
  fluidRow(class = "text-center",
           column(2, radioGroupButtons("dset", choices = c("1a", "1b", "2a", "2b", "2c", "3a" ), selected = "1a", status = "primary")),
           # column(2, radioGroupButtons("ext_objects", choices = c("2-object Task", "3-object Task"), status = "primary")),
           column(2, switchInput("wm_exp", "Working Memory Experiment", labelWidth = 300)),
           column(4, radioGroupButtons("display", choices = c("Time-Series Across Participants", "FFT Across Participants", "Time-Series + FFT by Individual"), selected = "FFT Across Participants", status = "primary"))
           ),
  br(), br(), br(), br(), 
  fluidRow(class = "text-center", 
           column(4, h3("Quantification and Statistical Analyses"), offset = 3)), br(), br(),
  fluidRow(column(4, br(),
                  switchInput("iso_sides", "Separate Hemifields", labelWidth = 150), br(),
                  switchInput("sbtr", "Analyze Invalid - Valid", labelWidth = 150), helpText("Subtract the dependent variable values at each CTI (valid - invalid) before performing analyses rather than analyzing valid and invalid trials independently"), br(),
                  numericInput("samp_per", "Sampling Period", min = round(1 / 60, 4), max = round(30 / 60, 4), value = round(3 / 60, 4), step = round(1 / 60, 4)), helpText(paste0("Spacing between CTI intevals (in seconds); the data was originally sampled at ", round(1 / 60, 4), ", but one could re-sample at a different rate, which would just clump neighboring CTI's together (whereas the below field groups neighbors but doesn't combine them, maintaining the total number of bins)")), br(),
                  # numericInput("samp_per", "Sampling Period", min = round(1 / 60, 4), max = round(30 / 60, 4), value = round(3 / 60, 4), step = round(1 / 60, 4)), helpText(paste0("Spacing between CTI intevals (in seconds); the data was originally sampled at ", round(1 / 60, 4), ", but one could re-sample at a different rate, which would just clump neighboring CTI's together (whereas the below field groups neighbors but doesn't combine them, maintaining the total number of bins)")), br(),
                  numericInput("clumps", "Neighbors to average at each CTI", min = 0, max = 14, value = 0, step = 2), helpText("`0` means this function does nothing, `2` means each CTI is the average of that CTI and its neighboring CTI's on each sides, etc...")
                  ),
            column(3,
                  radioGroupButtons("dep_var", "Dependent Variable", c("Accuracy", "Response Time"), selected = "Accuracy", status = "primary"), br(),
                  numericInput("α", "Alpha Value", max = .99, value = .05), helpText("The α to use for drawing the significance cutoff on the graphs"), br(),
                  numericInput("shuff", "Surrogate Shuffles for Null Hypothesis", min = 1, max = 10000, value = 50, step = 1), helpText("NOTE: increasing this number slows down the run time"), br(),
                  radioGroupButtons("mult_correcs", "Multiple Comparisons Correction", c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), selected = "BH", status = "primary")
            ),
           column(5, br(),
                  checkboxGroupButtons("trends", choices = c("Detrending", "Demeaning"), selected = c("Detrending", "Demeaning"), status = "primary"), br(),
                  radioGroupButtons("smooth_method", "Smoothing Individuals\" Time-Series", c("GLM", "GAM", "Loess", "LM"), selected = "Loess", status = "primary"), br(),
                  radioGroupButtons("win_func", "Windowing Function", c("Cosine", "Hamming", "Hann", "Kaiser", "Square", "Triangle", "Tukey", "Welch"), selected = "Tukey", status = "primary"), br(),
                  numericInput("xmax", "Max Hz Analyzed", min = 1, value = 15), br(),
                  numericInput("duration", "Duration (Seconds) Analyzed Including Padding", min = .8, max = 10, value = 1, step = .001)
           )),
  fluidRow(class = "text-center", 
          column(4, br(), br(), h3( "Filtering Participants"), offset = 3)), br(), br(),
  fluidRow(column(4, br(),
      switchInput("attn_filter", "Attention Filter", labelWidth = 100), helpText("Remove participants who reported either dozing off or not fully concentrating on multiple blocks"), br(),
      numericInput("catch_floor", "Catch-Trial Accuracy Floor", min = 0, max = 1, value = .85), helpText("Filter out participants whose filtered data is outside of the selected range")
  ),
  column(4, br(),
         numericInput("side_bias", "Side Bias Ceiling", min = .01, max = .99, value = .2), helpText("Filter out participants whose hit rate at one visual field - another visual field is > this ceiling"), br(),
         numericInput("wm_floor", "Working Memory Accuracy Floor", min = .01, max = .99, value = .7), helpText("Filter out participants whose working memory task accuracy is < this cutoff"), br(),
         numericInput("invalid_floor", "Invalid Trial Accuracy Floor", min = .01, max = .99, value = .02), helpText("Filter out participants whose accuracy on invalid trials is < this cutoff")
  ),
  column(4, br(),
         sliderInput("pre_range", "Pre-Filtered Accuracy Cutoffs", min = 0, max = 1, value = c(.45, .85)), helpText("Remove participants whose unfiltered accuracy is outside this range"), br(),
         sliderInput("post_range", "Post-Filtered Accuracy Cutoffs", min = 0, max = 1, value = c(.45, .85)), helpText("Remove participants whose filtered data is outside of the selected range"), br(),
         sliderInput("filtered_cap", "Maximum % Trials Filtered Out", min = 0, max = 1, value = .4), helpText("Remove participants whose trials are interpolated over more than percent of the time")
  )),
  fluidRow(class = "text-center", br(),
           column(4, h3( "Filtering Trials"), offset = 3)), br(), br(),
  fluidRow(column(4, br(),
                  sliderInput("block_range", "Block Accuracy Range", min = 0, max = 1, value = c(.20, .80)), helpText("Interpolate over trials if the average hit rate in that block, every 48 trials, was below this value"), br(),
                  checkboxGroupButtons("blocks_desired", "Blocks Included in Analysis", choices = 1:8, selected = 1:8, status = "primary")),
            column(4,
                   sliderInput("miniblock_range", "Mini-Block Accuracy Cutoffs", min = 0, max = 1, value = c(.40, .80)), helpText("Interpolate over trials if the average hit rate in that mini-block, every 16 trials which is how often the task difficulty was adjusted to titrate to 65%, is outside this range")),
           column(4, br(),
                  sliderInput("CTI_range", "Remove CTI\'s Outside This Range", min = .3, max = 1.29, value = c(.3, 1.09)))), br(), br())


server <- function(input, output, session) {
  
  grouping_cnsts <- quos(participant, Trials_filtered_out, Acc_prefilter, Acc_postfilter, CatchAcc)
  
  observe({
    batch <- substring(input$dset, 1, 1)
    batch_version <- substring(input$dset, 2, 2)
    if (batch == "1"){
      blocksize <- 54
      if (batch_version == "a"){
        pcpts <- 301:324
      } else if (batch_version == "b"){
        pcpts <- 401:427
      }
    } else if (batch == "2"){
      blocksize <- 80
      if (batch_version == "a"){
        pcpts <- c(501:522, 524:546)
      } else if (batch_version == "b"){
        pcpts <- 601:644
      } else if (batch_version == "c"){
        701:730
      }
    } else if(batch == "3"){
      blocksize <- 79
      pcpts <- 901:922
    }
    
    dep_var_abbr <- as.name(ifelse(input$dep_var == "Accuracy", "Acc", "RT"))
    
    dem_df <- fread(file.path("data", "Demographics.csv")) %>%
      mutate_at(vars(SubjID), as.numeric)
    
    # For each participant, function reads in data, filters it, and transforms it to prepare for interpolation and FFT'ing
    pcpts_combine <- function(pcpt){
      fread(file.path("data", pcpt, paste0(pcpt, ".csv"))) %>%# Reads in participant data (the `select` part is because PsychoPy created an empty column at the end of the data frame for the first few participants, which meant those dataframes had different dimensions than subsequent data frames, which `do.call` doesn't like)
        filter(Trial > 0) %>%                                                   # Prunes practice trials
        mutate(Stim_Sides = as.character(
                 ifelse(CorrSide == FlashSide, "Valid", "Invalid")),            # Creates column indicating whether cue was valid or invalid
               CTI = (lilsquareStartTime - flash_circleEndTime) %>%
                 RoundTo(1 / 60) %>% RoundTo(input$samp_per),
               block = RoundTo(Trial, blocksize, ceiling) / blocksize,          # Creates column indicating trial's block
               RT = ifelse(Acc == 1 & ButtonPressTime - lilsquareStartTime > .1,#                           RT after target appeared on screen, only for correct trials with an RT > 100 ms
                      ButtonPressTime - lilsquareStartTime, NA)) %>%
        mutate_at(vars(Acc, RT), list(~ifelse(block %in% input$blocks_desired,  # Converts trials' Acc and RT to NA if they are not in a desired block
                                              ., NA))) %>%
        mutate(wm_Acc = ifelse(session == "4",                                  # Creates column indicating wm accuracy; na.rm = TRUE because majority of trials lack wm probe, create NA
                               mean(Acc_wmarith, na.rm = TRUE), 1),
               CatchAcc = mean(ifelse(Opacity != 0, NA, Acc), na.rm = TRUE)) %>%#                           mean accuracy for catch trials
        filter(CatchAcc >= input$catch_floor,                                   # Prunes participants whose catch accuracy is below desired threshold
               Opacity != 0) %>%                                                #        catch trials
        filter(mean(Acc[Stim_Sides == "Invalid"],                                 #        participants whose invalid trial accuracy is...
                    na.rm = TRUE) >= input$invalid_floor) %>%                           #        ...below desired threshold
        mutate(CorrSide = case_when(                                              # Creates column indicating which side the target appeared on
          CorrSide == 1 ~ "Right", CorrSide == -1 ~ "Left", TRUE ~ "Bottom"),
               Stim_Sides = case_when(                                                 # Overwrites `Stim_Sides` column if `sep_vis_fidels` parameter == `TRUE`...
                 input$iso_sides ~ paste(CorrSide, Stim_Sides,  sep = "_"),                  # to include which side of screen target was on,...    
                 TRUE ~ Stim_Sides)) %>%                          # ... as well as whether it was valid with cue; if `iso_sides` == `FALSE`, leaves `StimSides` unchanged
        group_by(CorrSide) %>%
        mutate(Side_Acc = mean(Acc, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(Side_Diff = max(Side_Acc) - min(Side_Acc)) %>%
        # left_join(dem_df, by = c("participant" = "SubjID")) %>%
        filter(#grepl(ifelse(input$attn_filter, "task" , ""), Attentiveness),    # Prunes participants reporting lack of alertness on at least two blocks, when `attn_filter` == TRUE
               Side_Diff <= input$side_bias,                                    #                     whose hit rate at one visual field - another visual field is > `side_bias`
               wm_Acc >= input$wm_floor) %>%                                    #                           working memory task accuracy is < `wm_floor`
        mutate(Acc_prefilter = mean(Acc, na.rm = TRUE)) %>%                     # Creates column indicating mean accuracy before we've filtered for `block_range`, unlike `Acc_postfilter`
        filter(Acc_prefilter %>%                                                  # Prunes participants whose non-catch, pre-block-filtering accuracy is outside of desired range
                 between(input$pre_range[1], input$pre_range[2], incbounds = TRUE),
               between(CTI, min(input$CTI_range), max(input$CTI_range))) %>%    #        trials outside of desired CTI range
        group_by(block) %>%
        mutate(block_acc = mean(Acc)) %>%                                       # Creates column indicating block's mean accuracy
        ungroup() %>%
        mutate(rown = row_number(),
               miniblock = RoundTo(rown, 16, ceiling) / 16) %>%                 #                           trial's mini-block (every 16 non-catch trials the opacity was readjusted)
        group_by(miniblock) %>%
        mutate(miniblock_avg = mean(Acc)) %>%
        ungroup() %>%
        mutate_at(vars(Acc, RT), list(~ifelse(                                    # Changes `Acc` and `RT` column values to NA if...
          between(block_acc, input$block_range[1], input$block_range[2]) & between(           #  ... trial's block accuracy outside of `block_range`...
            miniblock_avg, input$miniblock_range[1], input$miniblock_range[2]), ., NA))) %>%
        # mutate_at(vars(Acc, RT), list(~ifelse(session != "4"  | Acc_wmarith == 1, ., NA))) %>%
        mutate(Trials_filtered_out = sum(is.na(Acc)) / n(),
               Acc_postfilter = mean(Acc, na.rm = TRUE)) %>%                    # Creates column indicating mean accuracy for non-catch trials; note this is after before we've filtered
        # for `block_range`, unlike `Acc_prefilter`
        filter(between(Acc_postfilter, input$post_range[1], input$post_range[2]),  # Prunes participants whose non-catch, post-block-filtering accuracy is outside of desired range 
               Trials_filtered_out <= input$filtered_cap) %>%
        group_by(CTI, Stim_Sides, !!!grouping_cnsts) %>%
        summarise_at(vars(Acc, RT), list(~mean(., na.rm = TRUE))) %>%           # Overwrites `Acc` and `RT` columns according to mean of each combination of `CTI` and `Stim_Sides`
        arrange(CTI) %>%
        group_by(Stim_Sides) %>%
        mutate_at(vars(Acc, RT), list(~na.approx(., na.rm = FALSE, rule = 2))) %>% # If any combination of `CTI` and `Stim_Sides` has only NA values, it takes on the average of its
        # neighboring CTI with same `Stim_Sides`
        mutate_at(vars(Acc, RT), list(~rollapply(., input$clumps + 1, mean,     # Averages each CTI with neighbors
                                                 partial = TRUE)))
    }
    
    cmbd <- map_dfr(pcpts, pcpts_combine) %>%                    # Calls `pcpts_combine` function for argumenet `pcpts`; then combines each participant's dataframe into one
      arrange(Acc_prefilter, participant, Stim_Sides, CTI)
    CTIs <- unique(cmbd$CTI)
    if (input$win_func == "Tukey"){ win <- tukeywindow(length(CTIs), .5)} else { # Creates window, which if `tukey` will add the parameter `r` == `.5` —so 'only' half the data length will be non-flat
      win <- match.fun(paste0(tolower(input$win_func), "window"))(length(CTIs))}
    locations <- unique(cmbd$Stim_Sides)                                   # Creates vector of column names representing sides locations of target in reference to cue (and also potentially side of screen)
    pcpts <- unique(cmbd$participant)                                       # Creates vector of remaining participant numbers after `pcpts_combine` filtering
    cmbd_w <- cmbd %>%
      pivot_wider(CTI:CatchAcc, Stim_Sides, values_from = !!dep_var_abbr) %>%
      arrange(Acc_prefilter, participant, CTI) %>%
      group_by(participant) %>%
      mutate_at(vars(locations), list(~na.approx(., na.rm = FALSE, rule = 2)))
    
    # Analyzes Invalid - Valid instead of them separetely, if sbtr == `Yes`
    if (input$sbtr){
      s <- tail(1:ncol(cmbd_w), length(locations)) [c(TRUE, FALSE)]
      cmbd_w[paste0(names(cmbd_w[s]), "_minus_",
                    names(cmbd_w)[s + 1])] <- cmbd_w[s] - cmbd_w[s + 1]
      locations = tail(colnames(cmbd_w), length(locations) / 2)
    }
    
    # Determines confidence intervals
    conf_int <- function(x, y){ x %>%
        summarise_at(vars(y), list(~qnorm(.975) * std_err(.)))
    }
    
    # Transforms from Time to Frequency Domain
    amplitude <- function(x, y, z, q, ...){
      pre_pad <- nrow(x)                                                          # Number of rows expected with one row per pcpt per CTI per shuffle
      x %>%
        ungroup() %>%
        mutate(samp_shuff = row_number() %>% RoundTo(pre_pad / y, ceiling) %>%    # Creates variable to track shuffle number which is then used for group_by
                 `*` (y / pre_pad)) %>%             
        group_by(participant, samp_shuff) %>%
        mutate_at(vars(locations),
                  list(~ case_when("Detrending" %in% input$trends ~                     # If `Detrending` selected...
                                     . - (polyfit(CTI, ., 2) %>% polyval(CTI)),        # detrend with this formula...
                                   TRUE ~ .))) %>%                                # otherwise ignore
        mutate_at(vars(locations), list(~case_when("Demeaning" %in% input$trends ~      # Works just like the detrending except for demeaning
                                                     . - mean(.), TRUE ~ .))) %>%
        mutate_at(vars(locations), list(~ . * win / Norm(win))) %>%               # Apply window
        ungroup() %>%
        add_row(participant = cmbd_w$CTI %>% range %>% diff %>% - input$duration %>%    # Add empty rows (other than subject ID) as additional CTI's needed to reach desired padded...
                  `/` (-samp_per) %>% + .00001 %>% floor %>% - 1 %>% `*` (y) %>%        # ... `duration` of intervals for each participant for each shuffle (`y` represents each shuffle)
                  rep(pcpts, .)) %>%
        mutate_at(vars(locations), list(~coalesce(., 0))) %>%
        mutate_at(vars(samp_shuff),                                               # Tags the padded rows with one of the shuffles created earlier with `samp_shuff` the non-padded rows...
                  list(~coalesce(., (row_number() - pre_pad) %>%                  # ... (<= pre_pad) appear first and have already been tagged, so we keep them as they are
                                   RoundTo((n() - pre_pad) / y, ceiling) %>%
                                   `/`((n() - pre_pad) / y)))) %>%
        group_by(participant, samp_shuff) %>%                                   # This group_by is critical so we're only taking the FFT of each shuffle
        mutate_at(vars(locations), list(~fft(.) %>%                               # Tabulate amplitude
                                          `*`(sqrt(2 / n())) %>%
                                          Mod %>% `^` (2))) %>%
        mutate(Hz = (row_number() - 1) / (n() * input$samp_per)) %>%            # Set Hz corresponding to each amplitude
        ungroup() %>%
        filter(dense_rank(Hz) - 1 <= floor(n_distinct(Hz) / 2),
               Hz < input$xmax) %>%
        select(-CTI) %>%
        pivot_longer(locations, z, values_to = q) %>%
        group_by(!!!syms(...))
    }
    
    amps <- function(w, ...){amplitude(cmbd_w, 1, "Location", w, ...)}
    
    fft_x <- 1 / input$duration
    xaxis_r <- RoundTo(input$xmax, fft_x)

    
    # Set Up Graphing
    
    t_srs_g <- function(x){
      cmbd_w %>%
        gather(Location, !!dep_var_abbr, -c(CTI, !!!grouping_cnsts)) %>%
        right_join(gather(conf_int(group_by(cmbd_w, CTI), locations),
                          Location, Conf_Int, -CTI),
                   by = c("CTI", "Location")) %>%
        group_by(CTI, Location, Conf_Int, !!!head(grouping_cnsts, x)) %>%
        summarise(!!dep_var_abbr := mean(!!dep_var_abbr)) %>%                               # Keeps either `RT` or `Acc` column—depending on whether `dep_var_abbr` parameter == `RT` or `Acc`
        ggplot(aes(CTI, !!dep_var_abbr, group = Location, color = Location,
                   fill = Location, ymin = !!dep_var_abbr - Conf_Int,
                   ymax = !!dep_var_abbr + Conf_Int)) + 
        labs(title = paste(dep_var_abbr, "by Cue-Target Interval, Exp.", input$dset),
             x = "Cue-Target Interval (ms)")
    }
    fft_g <- function(x) { x +
        labs(title = paste0("FFT of Target ", input$dep_var, ", Exp. ", input$dset),
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
                              labels = sapply(locations, simplify = TRUE, function(x) gsub("_", " ", x),
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
                   ncol = RoundTo(length(pcpts) / 4, 1, ceiling),  scales = "free_x")                                   # `free` means the y_axis isn't fixed from participant to participant
    }
    cmbd_g <- function(x, y) {
      graph(y, x) + geom_line(size = 1.5) +
        theme(legend.key.size = unit(.55, "in")) +
        labs(subtitle = paste( "Data from", as.character(length(pcpts)),
                               "participants"))
    }
    
    if (input$display == "Time-Series + FFT by Individual"){
      # Produces left half of final graph
      ts_facets <- idvl_g(1, t_srs_g) +
        geom_line(alpha = I(2 / 10), color = "grey", show.legend = FALSE) +       # Graphs unsmoothed data in light gray
        stat_smooth(method = tolower(input$smooth_method), span = 0.2, se = FALSE,               # Smoothes data depending on `smooth_method` parameter
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
      fft_facets <- amps("Power0", c("participant", "Hz")) %>%
        summarise_all(mean) %>%
        gather(Flash_and_or_field, Power, -Hz, -samp_shuff,
               -c(!!!grouping_cnsts)) %>%
        ggplot(aes(Hz, Power, color = Flash_and_or_field)) %>%
        idvl_g(fft_g) +
        geom_line() +
        scale_x_continuous(
          name = "Frequency (Hz)", limits = c(0, input$xmax), breaks = seq(
            0, input$xmax, ifelse(input$xmax > 10 | input$xmax != xaxis_r,
                            xaxis_r %>% `/`(seq(fft_x, xaxis_r, fft_x)) %>%
                              Closest(5) %>% `/` (xaxis_r) %>% max %>% `^` (-1),
                            fft_x))) +
        labs(caption = paste("Data from", as.character(length(pcpts)),
                             "participants")) +
        geom_text(data = as.data.frame(plot_label), inherit.aes = FALSE, size = 2.5,# Sets location for label overlayed onto graph
                  aes(label = lab, x = Inf, y = Inf), vjust = 1.15, hjust = 1.05)
      output$mygraph <- renderPlot({ grid.arrange(ts_facets, fft_facets, ncol = 2)
      }, height = 800, width = 1350)                # Combines time series and FFT graphs into one plot
    }  else if (input$display == "Time-Series Across Participants") {
      output$mygraph <- renderPlot({
        (cmbd_g(0, t_srs_g) +
           theme(panel.grid = element_blank()) +
           geom_ribbon(alpha = 0.15, aes(color = NULL))) %>%
          move_layers("GeomRibbon", position = "bottom")
      }, height = 800, width = 1350)
    } else { # Graph combined FFT ------------------------------------------------
      
      # Produces `shuff` # of null hypothesis permutations
      amps_shuff <- map_dfr(1:input$shuff, function(x){
        set.seed(x)
        cmbd_w %>%
          group_by(participant) %>%                                               # The three lines randomize row order of CTI's for each participant; `group_by` is for  shuffling only...
          sample_n(n()) %>%                                                       # ... within rows with the same particpant number; `mutate_at` below allows shuffling only for CTI...
          mutate_at(vars(CTI), list(~seq(min(CTIs), max(CTIs), samp_per)))        # ... column; specifically, keeps rows intact except resets CTI's in descending order, even though...
      }) %>%                                                                      # ... previous line was just randomizing, thereby randomizing CTI vs. performance
        amplitude(shuff, "Location", "Power",                                     # Tabulates amplitude for each participant's data for `shuff` number of shuffles
                  c("Hz", "samp_shuff", "Location")) %>%
        summarise_at(vars(Power), mean) %>%                                       # Finds average amplitude at each CTI across all participants per shuffle
        group_by(Hz, Location) %>%
        mutate_at(vars(Power), list(ntile = ~1.02 - ecdf(.)(.))) 
      
      amp_and_shf <- amps("Power0", c("Hz", "Location")) %>%
        summarise_at(vars(Power0), mean) %>%
          ungroup() %>%
          mutate(ntile = pmap_dbl(., function(Hz, Location, Power0, ...){           # For each row in `amps`...
            amps_shuff %>%                                                          # ... finds row in `amps_shuff` that...
              rename(Hz0 = Hz, Location0 = Location) %>%
              filter(Hz0 == Hz,                                                     # ... 1) matches in the Hz column...
                     Location0 == Location) %>%      
              ungroup() %>%
              top_n(1, -abs(Power - Power0)) %>%                # check why this pmap function is giving me slightly different results than what I'd written before
              select(ntile) %>%
              as.double()
          })) %>%
          group_by(Location) %>%
          mutate_at(vars(ntile), list(~p.adjust(., method = input$mult_correcs) / .)) %>%
          select(Hz, Location, ntile) %>%
          safe_right_join(amps_shuff, by = c("Hz", "Location"),
                          conflict = `*`) %>%
          group_by(Hz, Location) %>%
          filter(abs(ntile - input$α) == min(abs(ntile - input$α))) %>%
          ungroup() %>%
          select(-c(samp_shuff, ntile)) %>%
          combine(amps("Power", c("Location", "Hz")) %>%                        # Finds average amplitude at each CTI for real (not surrogate) data, then merges that data with surrogate
                    summarise(Power = mean(Power)) %>%
                      ungroup(),
                  names = c("Significance Cutoff", "Observed Data")) %>%
          right_join(conf_int(amps("Conf_Int", c("Hz", "Location")), "Conf_Int"),
                     by = c("Hz", "Location"))
      
        output$mygraph <- renderPlot({
          amp_and_shf2 <- amp_and_shf %>%
            ggplot(aes(Hz, Power, col = Location,
                       linetype = source, fill = Location,
                       ymin = Power - Conf_Int,
                       ymax = Power + Conf_Int)) %>%
            cmbd_g(fft_g) +
            scale_linetype_manual(values = c("solid", "dashed")) +
            scale_x_continuous(name = "Frequency (Hz)",
                               limits = c(0, input$xmax),
                               breaks = seq(0, input$xmax, ifelse(fft_x > .5,
                                                            round(fft_x, 2), 1))) +
            labs(linetype = "",
                 caption = paste("Significance threshold at p < ",
                                 as.character(input$α))) +
            geom_ribbon(data = filter(amp_and_shf, source == "Observed Data"),
                        alpha = 0.15, aes(color = NULL)) +
            geom_point(size = 3,
                       data = amp_and_shf %>%
                         spread(source, Power) %>%
                         filter(`Observed Data` > `Significance Cutoff`) %>%
                         select(-c(`Significance Cutoff`, Conf_Int)) %>%
                         gather(source, Power, -Hz, -Location), 
                       aes(ymin = NULL, ymax = NULL))
          move_layers(amp_and_shf2, "GeomRibbon", position = "bottom")
        }, height = 800, width = 1350)
      
    }
    
  })

}
shinyApp(ui, server)
