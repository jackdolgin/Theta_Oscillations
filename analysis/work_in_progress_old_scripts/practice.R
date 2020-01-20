# if (!require(devtools)) install.packages('devtools')
# if (!require(smisc)) devtools::install_github("stevenworthington/smisc")
# smisc::ipak(c("utils", "tidyr", "dplyr", "ggplot2", "DescTools", "bspec", "pracma", "gridExtra", "data.table", "tables", "zoo", "parallel", "scales", "lazyeval", "stats", "gdata", "viridis", "gginnards", "shiny", "shinyWidgets"))

library(utils); library(tidyr); library(dplyr); library(ggplot2); library(DescTools); library(bspec); library(pracma); library(gridExtra); library(data.table); library(tables); library(zoo); library(parallel); library(scales); library(lazyeval); library(stats); library(gdata); library(viridis); library(gginnards); library(shiny); library(shinyWidgets)

ui <- fluidPage(
  title = "Theta Oscillations Analysis",
  plotOutput("mygraph",  height = 800, width = 1350),
  hr(),       
  fluidRow(class = "text-center",
           column(4, h3( "Data Set, Task, and Graph Choice"), offset = 3)), br(), br(),
  fluidRow(class = "text-center",
           column(2, radioGroupButtons("dset", choices = c("Pilot", "Experimental"), selected = "Experimental", status = "primary")),
           column(2, radioGroupButtons("ext_objects", choices = c("2-object Task", "3-object Task"), status = "primary")),
           column(2, switchInput("wm_exp", "Working Memory Experiment", labelWidth = 300)),
           column(4, radioGroupButtons("display", choices = c("Time-Series Across Participants", "FFT Across Participants", "Time-Series + FFT by Individual"), selected = "FFT Across Participants", status = "primary"))
  ),
  br(), br(), br(), br(), 
  fluidRow(class = "text-center", 
           column(4, h3("Quantification and Statistical Analyses"), offset = 3)), br(), br(),
  fluidRow(column(4, br(),
                  switchInput("iso_sides", "Separate Hemifields", labelWidth = 150), br(),
                  switchInput("sbtr", "Analyze Invalid - Valid", labelWidth = 150), helpText("Subtract the dependent variable values at each CTI (valid - invalid) before performing analyses rather than analyzing valid and invalid trials independently"), br(),
                  numericInput("samp_per", "Sampling Period", min = round(1 / 60, 4), max = round(30 / 60, 4), value = round(1 / 60, 4), step = round(1 / 60, 4)), helpText(paste0("Spacing between CTI intevals (in seconds); the data was originally sampled at ", round(1 / 60, 4), ", but one could re-sample at a different rate, which would just clump neighboring CTI's together (whereas the below field groups neighbors but doesn't combine them, maintaining the total number of bins)")), br(),
                  numericInput("clumps", "Neighbors to average at each CTI", min = 0, max = 14, value = 2, step = 2), helpText("`0` means this function does nothing, `2` means each CTI is the average of that CTI and its neighboring CTI's on each sides, etc...")
  ),
  column(3,
         radioGroupButtons("dep_var", "Dependent Variable", c("Accuracy", "Response Time"), selected = "Accuracy", status = "primary"), br(),
         numericInput("pval", "P-value", max = .99, value = .05), helpText("The p-value to use for drawing the significance cutoff on the graphs"), br(),
         numericInput("shuff", "Surrogate Shuffles for Null Hypothesis", min = 1, max = 10000, value = 50, step = 1), helpText("NOTE: increasing this number slows down the run time")
  ),
  column(5, br(),
         checkboxGroupButtons("trends", choices = c("Detrending", "Demeaning"), selected = c("Detrending", "Demeaning"), status = "primary"), br(),
         radioGroupButtons("smooth_method", "Smoothing Individuals\" Time-Series", c("GLM", "GAM", "Loess", "LM"), selected = "Loess", status = "primary"), br(),
         radioGroupButtons("win_func", "Windowing Function", c("Cosine", "Hamming", "Hann", "Kaiser", "Square", "Triangle", "Tukey", "Welch"), selected = "Tukey", status = "primary"), br(),
         numericInput("xaxisvals", "Max Hz Displayed", min = 1, value = 15), br(),
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
         numericInput("wm_floor", "Working Memory Accuracy Floor", min = .01, max = .99, value = .7), helpText("Filter out participants whose working memory task accuracy is < this cutoff")
  ),
  column(4, br(),
         sliderInput("pre_range", "Pre-Filtered Accuracy Cutoffs", min = 0, max = 1, value = c(.45, .85)), helpText("Remove participants whose unfiltered accuracy is outside this range"), br(),
         sliderInput("post_range", "Post-Filtered Accuracy Cutoffs", min = 0, max = 1, value = c(.45, .85)), helpText("Remove participants whose filtered data is outside of the selected range")
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
    pcpts <- if (input$dset == "Pilot") {
      blocksize <- 54
      if(input$ext_objects == "2-object Task") 301:324 else 401:427
    } else {
      blocksize <- 80
      if (input$wm_exp){
        c(701:726)
      } else {
        if(input$ext_objects == "2-object Task") c(501:522, 524:534) else c(601:631)}
    }
    
    dep_var_abbr <- as.name(ifelse(input$dep_var == "Accuracy", "Acc", "RT"))
    
    dem_df <- fread(file.path("data", "Demographics.csv")) %>%
      mutate_at(vars(SubjID), as.numeric)
    
    # For each participant, function reads in data, filters it, and transforms it to prepare for interpolation and FFT'ing
    pcpts_combine <- function(pcpt){
      fread(file.path("data", pcpt, paste0(pcpt, ".csv"))) %>%# Reads in participant data (the `select` part is because PsychoPy created an empty column at the end of the data frame for the first few participants, which meant those dataframes had different dimensions than subsequent data frames, which `do.call` doesn't like)
        filter(Trial > 0) %>%                                                     # Filters out practice trials
        mutate(CatchAcc = mean(ifelse(Opacity != 0, NA, Acc), na.rm = TRUE)) %>%   # Creates column indicating mean accuracy for catch trials
        left_join(dem_df, by = c("participant" = "SubjID")) %>%
        filter(CatchAcc >= input$catch_floor,                                   # Filters out participants whose catch accuracy is below desired threshold
               # grepl(ifelse(input$attn_filter, "task" , ""), Q9),
               
               Opacity != 0) %>%
        group_by(CorrSide) %>%
        mutate(Side_Acc = mean(Acc, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(Side_Diff = max(Side_Acc) - min(Side_Acc)) %>%
        filter(Side_Diff <= input$side_bias) %>%
        mutate(Acc_prefilter = mean(Acc, na.rm = TRUE),                           # Creates column indicating mean accuracy before we've filtered for `block_range`, unlike `Acc_postfilter`
               CTI = RoundTo(RoundTo(lilsquareStartTime - flash_circleEndTime,
                                     1 / 60), input$samp_per)) %>% 
        filter(between(Acc_prefilter, input$pre_range[1], input$pre_range[2],   # Filters out participants whose non-catch, pre-block-filtering accuracy is outside of desired range
                       incbounds = TRUE),                       
               between(CTI,  min(input$CTI_range), max(input$CTI_range))) %>%
        mutate(block = RoundTo(Trial, blocksize, ceiling) / blocksize,                          # Creates column indicating trial's block
               RT = ifelse(Acc == 1 & ButtonPressTime - lilsquareStartTime > .1,  #                indicating RT after target appeared on screen, only for correct trials with an RT < 100 ms
                           ButtonPressTime - lilsquareStartTime, NA),
               Stim_Sides = as.character(
                 ifelse(CorrSide == FlashSide, "Valid", "Invalid")),              #                indicating whether cue was valid or invalid
               CorrSide = case_when(CorrSide == 1 ~ "Right",                      #                indicating which side the target appeared on
                                    CorrSide == -1 ~ "Left", TRUE ~ "Bottom"),
               Stim_Sides = case_when(input$iso_sides ~ paste(CorrSide, Stim_Sides, sep = "_"), # Overwrites `Stim_Sides` column if `sep_vis_fidels` parameter == `TRUE` to include which side of screen target was on,
                                      TRUE ~ Stim_Sides)) %>%             # as well as whether it was valid with cue; if `iso_sides` parameter == `No`, leaves `StimSides` unchanged
        group_by(block) %>%
        mutate(block_acc = mean(Acc)) %>%                                         # Creates column indicating block's mean accuracy
        ungroup() %>%
        mutate(rown = row_number(),
               miniblock = RoundTo(rown, 16, ceiling) / 16) %>%                   # Creates column indicating trial's mini-block (every 16 trials the opacity was readjusted)
        group_by(miniblock) %>%
        mutate(miniblock_avg = mean(Acc)) %>%
        ungroup() %>%
        mutate_at(vars(Acc, RT),
                  list(~ifelse(!between(block_acc, input$block_range[1],                # Changes `Acc` and `RT` column values to NA if trial's
                                        input$block_range[2]) |                         # block accuracy outside of `block_range`
                                 !between(miniblock_avg, input$miniblock_range[1],      # or miniblock was not in the desired range
                                          input$miniblock_range[2]), NA, .))) %>%
        mutate(Trials_filtered_out = sum(is.na(Acc)) / n(),
               Acc_postfilter = mean(Acc, na.rm = TRUE)) %>%                      # Creates column indicating mean accuracy for non-catch trials; note this is after before we've filtered for `block_range`, unlike `Acc_prefilter`
        filter(between(Acc_postfilter, input$post_range[1], input$post_range[2])) %>%                # Filters out participants whose non-catch, post-block-filtering accuracy is outside of desired range
        group_by(CTI, Stim_Sides, !!!grouping_cnsts) %>%
        summarise_at(vars(Acc, RT), list(~mean(., na.rm = TRUE))) %>%              # Overwrites `Acc` and `RT` columns according to mean of each combination of `CTI` and `Stim_Sides`
        arrange(CTI) %>%
        group_by(Stim_Sides) %>%
        mutate_at(vars(Acc, RT), list(~na.approx(., na.rm = FALSE, rule = 2))) %>%
        mutate_at(vars(Acc, RT), list(~rollapply(., input$clumps + 1, mean, partial = TRUE)))
    }
    
    cmbd <- do.call(rbind, lapply(pcpts, pcpts_combine)) %>%             # Calls `pcpts_combine` function for argumenet `pcpts`; then combines each participant's dataframe into one
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
    conf_int <- function(x, ...){ x %>%
        group_by(.dots = lazy_dots(...)) %>%
        summarise_at(vars(locations), list(~qnorm(.975) * std_err(.)))
    }
    
    # Transforms from Time to Frequency Domain
    amplitude <- function(x, y){
      pre_pad <- length(pcpts) * (length(unique(cmbd_w$CTI))) * y
      x %>%
        group_by(participant) %>%
        mutate_at(vars(locations), list(~case_when("Detrending" %in% input$trends ~ . - polyval(polyfit(CTI, ., 2), CTI), TRUE ~ .))) %>%
        mutate_at(vars(locations), list(~ case_when("Demeaning" %in% input$trends ~ . - mean(.), TRUE ~ .))) %>%
        mutate_at(vars(locations), list(~ . * win / Norm(win))) %>%
        ungroup() %>%
        add_row(participant = rep(pcpts, y * (ceiling((input$duration - diff(range(cmbd_w$CTI)))/input$samp_per)))) %>%
        head(-length(pcpts) * y) %>%
        mutate_at(vars(locations), list(~coalesce(., 0))) %>%
        mutate(samp_shuff = ifelse(row_number() <= pre_pad, 
                                   RoundTo(row_number(), pre_pad / y,
                                           ceiling) / (pre_pad / y), 
                                   RoundTo(row_number() - pre_pad, (n() - pre_pad) / (y),
                                           ceiling) / ((n() - pre_pad) / y))) %>%
        group_by(participant, samp_shuff) %>%
        mutate_at(vars(locations), list(~Mod(sqrt(2/n()) * fft(.))^2)) %>%
        mutate(Hz = round((row_number() - 1) / (n() * input$samp_per), 2)) %>%
        ungroup() %>%
        select(-CTI)
    }
    
    amps <- amplitude(cmbd_w, 1)
    fft_x <- round(1 / (length(unique(amps$Hz)) * input$samp_per),1)
    xaxis_r <- RoundTo(input$xaxisvals, fft_x)
    
    
    # Set Up Graphing
    
    t_srs_g <- function(x){
      cmbd_w %>%
        gather(Location, !!dep_var_abbr, -c(CTI, !!!grouping_cnsts)) %>%
        right_join(gather(conf_int(cmbd_w, CTI), Location, Conf_Int, -CTI),
                   by = c("CTI", "Location")) %>%
        group_by(CTI, Location, Conf_Int, !!!head(grouping_cnsts, x)) %>%
        summarise(!!dep_var_abbr := mean(!!dep_var_abbr)) %>%                               # Keeps either `RT` or `Acc` column—depending on whether `dep_var_abbr` parameter == `RT` or `Acc`
        ggplot(aes(CTI, !!dep_var_abbr, group = Location, color = Location,
                   fill = Location, ymin = !!dep_var_abbr - Conf_Int,
                   ymax = !!dep_var_abbr + Conf_Int)) + 
        labs(title = paste0(input$dep_var, " by Cue-Target Interval ", input$ext_objects),
             x = "Cue-Target Interval (ms)")
    }
    fft_g <- function(x) { x +
        labs(title = paste0("FFT of Target ", input$dep_var, ", ", input$ext_objects),
             col = "Target Location", y = "Spectral Power") +
        theme(panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_blank())
    }
    
    viridis_cols <- .7 + RoundTo(.0001 * RoundTo(length(locations),
                                                 4, floor), .2, ceiling)
    graph <- function(y, x) {
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
    
    idvl_g <- function(y, x){
      graph(y, x) +
        theme(plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5)) +
        facet_wrap(~factor(participant, levels = pcpts),
                   ncol = RoundTo(length(pcpts) / 4, 1, ceiling),  scales = "free_x")                                   # `free` means the y_axis isn't fixed from participant to participant
    }
    cmbd_g <- function(y, x) {
      graph(y, x) + geom_line(size = 1.5) +
        theme(legend.key.size = unit(.55, "in")) +
        labs(subtitle = paste( "Data from", as.character(length(pcpts)),
                               "participants"))
    }
    
    if (input$display == "Time-Series + FFT by Individual"){
      # Produces left half of final graph
      ts_facets <- idvl_g(t_srs_g, 1) +
        geom_line(alpha = I(2 / 10), color = "grey", show.legend = FALSE) +       # Graphs unsmoothed data in light gray
        stat_smooth(method = tolower(input$smooth_method), span = 0.2, se = FALSE,               # Smoothes data depending on `smooth_method` parameter
                    size = .5, show.legend = FALSE)
      
      # Produces label for each right-side graph
      plot_label <- amps %>%
        group_by(!!!grouping_cnsts) %>%
        summarise() %>%
        ungroup() %>%
        drop_na() %>%
        mutate_at(vars(!!!tail(grouping_cnsts, -1)),
                  funs(paste(quo_name(quo(.)), "=", percent(.)))) %>%
        unite(lab, !!!tail(grouping_cnsts, -1), sep = "\n", remove = FALSE)
      
      # Produces right half of final graph
      fft_facets <- idvl_g(fft_g, amps %>%
                             group_by(participant, Hz) %>%
                             summarise_all(mean) %>%
                             gather(Flash_and_or_field, Power, -Hz, -samp_shuff, -c(!!!grouping_cnsts)) %>%
                             ggplot(aes(Hz, Power, color = Flash_and_or_field))) +
        geom_line() +
        scale_x_continuous(name = "Frequency (Hz)", limits = c(0, xaxis_r),
                           breaks = seq(0, input$xaxisvals, ifelse(input$xaxisvals > 10 | input$xaxisvals != xaxis_r, 1/max(Closest(xaxis_r/ seq(fft_x, xaxis_r, fft_x), 5) /xaxis_r), fft_x))) +
        labs(caption = paste("Data from", as.character(length(pcpts)),
                             "participants")) +
        geom_text(data = as.data.frame(plot_label), inherit.aes = FALSE, size = 2.5,# Sets location for label overlayed onto graph
                  aes(label = lab, x = Inf, y = Inf), vjust = 1.15, hjust = 1.05)
      output$mygraph <- renderPlot({ grid.arrange(ts_facets, fft_facets, ncol = 2)
      }, height = 800, width = 1350)                # Combines time series and FFT graphs into one plot
    }  else if (input$display == "Time-Series Across Participants") {
      output$mygraph <- renderPlot({
        move_layers(cmbd_g(t_srs_g, 0) +
                      theme(panel.grid = element_blank()) +
                      geom_ribbon(alpha = 0.15, aes(color = NULL)), "GeomRibbon", position = "bottom")
      }, height = 800, width = 1350)
    } else { # Graph combined FFT ------------------------------------------------
      
      set.seed(123)
      
      # Produces `shuff` # of null hypothesis permutations
      shuffle <- function(x){
        cmbd_w %>%
          group_by(participant) %>%
          sample_n(length(CTIs), weight = CTI) %>%
          mutate_at(vars(CTI), list(~round(seq(min(CTIs), max(CTIs), input$samp_per), 2))) %>%
          mutate(samp_shuff = x)
      }
      
      # Produces and save graph
      amps_shuff <- do.call(rbind, lapply(1:input$shuff, shuffle)) %>%
        amplitude(input$shuff) %>%
        group_by(Hz, samp_shuff) %>%
        summarise_at(vars(locations), mean) %>%
        group_by(Hz) %>%
        summarise_at(vars(locations), list(~quantile(., probs = 1 - input$pval))) %>%
        combine(amps %>% group_by(Hz) %>% summarise_at(vars(locations), mean),
                names = (c("Significance Cutoff", "Observed Data"))) %>%
        gather(Location, Power, -c(Hz, source)) %>%
        right_join(gather(conf_int(amps, Hz), Location, Conf_Int, -Hz),
                   by = c("Hz", "Location"))
      output$mygraph <- renderPlot({
        (move_layers(cmbd_g(fft_g, ggplot(amps_shuff, aes(Hz, Power, col = Location, linetype = source, 
                                                          ymin = Power - Conf_Int, ymax = Power + Conf_Int, fill = Location))) +
                       scale_linetype_manual(values = c("solid", "dashed")) +
                       scale_x_continuous(name = "Frequency (Hz)", limits = c(0, input$xaxisvals),
                                          breaks = seq(0, input$xaxisvals, ifelse(fft_x > .5, round(fft_x, 2), 1))) +
                       labs(linetype = "",
                            caption = paste("Significance threshold at p < ",
                                            as.character(input$pval))) +
                       geom_ribbon(data = filter(amps_shuff, source == "Observed Data"), alpha = 0.15, aes(color = NULL)) +
                       geom_point(size = 3, data = amps_shuff %>% spread(source, Power) %>%
                                    filter(`Observed Data` > `Significance Cutoff`) %>%
                                    select(-c(`Significance Cutoff`, Conf_Int)) %>%
                                    gather(source, Power, -Hz, -Location), 
                                  aes(ymin = NULL, ymax = NULL)), "GeomRibbon", position = "bottom"))
      }, height = 800, width = 1350)
      
    }
    
  })
  
}
shinyApp(ui, server)
