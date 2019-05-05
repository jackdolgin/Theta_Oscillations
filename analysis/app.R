# if (!require(devtools)) install.packages('devtools')
# if (!require(smisc)) devtools::install_github("stevenworthington/smisc")
# smisc::ipak(c("utils", "tidyr", "dplyr", "ggplot2", "DescTools", "bspec", "pracma", "gridExtra", "data.table", "tables", "zoo", "parallel", "scales", "purrr", "lazyeval", "stats", "gdata", "viridis", "gginnards", "shiny"))

library(utils); library(tidyr); library(dplyr); library(ggplot2); library(DescTools); library(bspec); library(pracma); library(gridExtra); library(data.table); library(tables); library(zoo); library(parallel); library(scales); library(purrr); library(lazyeval); library(stats); library(gdata); library(viridis); library(gginnards); library(shiny)

ui <- fluidPage(
  title = "Theta Oscillations Analysis",
  plotOutput('mygraph',  height = 800, width = 1350),
  hr(),
      fluidRow(
        column(4,
               selectizeInput('ext_objects', 'Task by external objects', choices = c("2", "3")),  helpText("`2` corresponds to the two-object task, `3` to the three-object task"), br(),
               selectizeInput('dep_var', 'Dependent Variable', choices = c("Acc", "RT")), br(),
               selectizeInput('iso_sides', 'Group by Side', choices = c('No', 'Yes')), helpText("Group by not only valid and invalid but also by the side of the screen for each trial (i.e. going from `Valid` and `Invalid` to `Right Valid`, `Left Valid`, `Right Invalid`, `Left Invalid`)"), br(),
               selectizeInput('sbtr', 'Subtract Invalid from Valid', choices = c('No', 'Yes')), helpText("Subtract the dependent variable values at each CTI (valid - invalid) before performing analyses rather than analyzing valid and invalid trials independently"), br(),
               numericInput('pval', 'P-value', max = .99, value = .05), helpText('The p-value to use for drawing the significance cutoff on the graphs'), br(),
               numericInput('shuff', 'Surrogate Shuffles for Null Hypothesis', min = 1, max = 10000, value = 50, step = 1), helpText("NOTE: increasing this number slows down the run time"), br()),
        column(4,
               selectizeInput('display', 'Graph', choices = c("FFT Across Participants", "Time-Series Across Participants", "FFT + Time-Series By Individual")), br(),
               numericInput('latestart', 'Late Start', min = 0, max = .8, value = 0), helpText("Filters out trials within the first `latestart` seconds of CTI bins"), br(),
               numericInput('earlyend', 'Early End', min = 0, max = .8, value = 0), helpText('Filters out trials within the last `earlyend` seconds of CTI bins'), br(),
               numericInput('clumps', 'Points to average at each CTI', min = 1, max = 15, value = 1, step = 2), helpText("`1` means this function does nothing, `3` means each CTI is the average of that CTI and its neighboring CTI's on each sides, etc..."), br(),
               numericInput('samp_per', 'Sampling Period', min = round(1 / 60, 2), max = round(30 / 60, 2), value = round(1 / 60, 2), step = round(1 / 60, 2)), helpText("Spacing between CTI intevals (in seconds); the data was originally sampled at .02, but one could re-sample at a different rate, which would just clump neighboring CTI's together (whereas the above slider groups neighbors but doesn't combine them, keeping the same total number of bins)"), br(),
               numericInput('zeropads', 'Zero Pads', min = 0, max = 2000, value = 0, step = 1), helpText("The number of zero pads to use; can be set from 0 upwards"), br(),
               numericInput('xaxisvals', 'Greatest x-value Included in Graph', min = 1, value = 10), br(),
               selectizeInput('win_func', 'Windowing Function', choices = c("tukey", "square", "hann", "welch", "triangle", "hamming", "cosine", "kaiser")), br(),
               selectizeInput('smooth_method', 'Smoothing Methods for Time-Series Data of Individuals', choices = c("loess", "lm", "glm", "gam"))
 ),
        column(4,
               selectizeInput('attn_filter', 'Attention Filter', choices = c("Off", "On")), helpText("Either `Off` or `On`, which filters out any participants who indicated in the post-task questionnaire that they either dozed off at one point or were not fully focused on at least two of the eight blocks (the other options were that they were fully alert on all blocks or fully alert on all but one block)"), br(),
               sliderInput('pre_range', 'Pre-Filtered Accuracy Cutoffs', min = 0, max = 1, value = c(.45, .75)), helpText("Filter out participants whose unfiltered/scrutinized data is outside of the selected range"), br(),
               sliderInput('post_range', 'Post-Filtered Accuracy Cutoffs', min = 0, max = 1, value = c(.45, .75)), helpText("Filter out participants whose filtered data is outside of the selected range"), br(),
               numericInput('catch_floor', 'Catch-Trial Accuracy Floor', min = 0, max = 1, value = .85), helpText("Filter out participants whose filtered data is outside of the selected range"), br(),
               numericInput('block_floor', 'Block Accuracy Floor', min = 0, max = 1, value = .65), helpText("Interpolate over trials if the average hit rate in that block, every 48 trials, was below this value"), br(),
               sliderInput('miniblock_range', 'Mini-Block Accuracy Cutoffs', min = 0, max = 1, value = c(.45, .85)), helpText("Interpolate over trials if the average hit rate in that mini-block, every 16 trials which is how often the task difficulty was adjusted to titrate to 65%, is outside this range"), br(),
               numericInput('side_bias', 'Side Bias Ceiling', min = .01, max = .99, value = .3), helpText("Filter out participants whose hit rate at one visual field - another visual field is > `side_bias`"), br()
)))


server <- function(input, output, session) {
  
  grouping_cnsts <- quos(participant, Trials_filtered_out, Acc_prefilter,       # Columns that are frequently used for grouping, variable means
                         Acc_postfilter, CatchAcc)                              # don't have to type them out every time we use them for grouping
  
  observe({
    pcpts <- if(input$ext_objects == 2) 301:324 else 401:427

    dep_var <- as.name(input$dep_var)
    
    dem_df <- fread(file.path("data", "Demographics.csv")) %>%
      mutate_at(vars(SubjID), as.numeric)
    
    # For each participant, function reads in data, filters it, and transforms it to prepare for interpolation and FFT'ing
    pcpts_combine <- function(pcpt){
      fread(file.path("data", pcpt, paste0(pcpt, ".csv"))) %>%# Reads in participant data (the `select` part is because PsychoPy created an empty column at the end of the data frame for the first few participants, which meant those dataframes had different dimensions than subsequent data frames, which `do.call` doesn't like)
        filter(Trial > 0) %>%                                                     # Filters out practice trials
        mutate(CatchAcc = mean(ifelse(Opacity > 0, NA, Acc), na.rm = TRUE)) %>%   # Creates column indicating mean accuracy for catch trials
        group_by(CorrSide) %>%
        mutate(Side_Acc = mean(Acc, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(Side_Diff = max(Side_Acc) - min(Side_Acc)) %>%
        left_join(dem_df, by = c("participant" = "SubjID")) %>%
        filter(CatchAcc >= input$catch_floor,                                   # Filters out participants whose catch accuracy is below desired threshold
               grepl(ifelse(input$attn_filter == "On", "fully alert" , ""), Q9),
               Opacity > 0,
               Side_Diff <= input$side_bias) %>%
        mutate(Acc_prefilter = mean(Acc, na.rm = TRUE),                           # Creates column indicating mean accuracy before we've filtered for `block_floor`, unlike `Acc_postfilter`
               CTI = RoundTo(RoundTo(lilsquareStartTime - flash_circleEndTime,
                                     1 / 60), input$samp_per)) %>% 
        filter(between(Acc_prefilter, input$pre_range[1], input$pre_range[2], incbounds = TRUE),                       # Filters out participants whose non-catch, pre-block-filtering accuracy is outside of desired range
               between(CTI, min(CTI) + input$latestart, max(CTI) - input$earlyend)) %>%
        mutate(block = RoundTo(Trial, 54, ceiling) / 54,                          # Creates column indicating trial's block
               RT = ifelse(Acc == 1 & ButtonPressTime - lilsquareStartTime > .1,  #                indicating RT after target appeared on screen, only for correct trials with an RT < 100 ms
                           ButtonPressTime - lilsquareStartTime, NA),
               Stim_Sides = as.character(
                 ifelse(CorrSide == FlashSide, "Valid", "Invalid")),              #                indicating whether cue was valid or invalid
               CorrSide = case_when(CorrSide == 1 ~ "Right",                      #                indicating which side the target appeared on
                                    CorrSide == -1 ~ "Left", TRUE ~ "Bottom"),
               Stim_Sides = case_when(input$iso_sides == "No" ~ Stim_Sides,        # Overwrites 'Stim_Sides' column if `sep_vis_fidels` parameter == `Yes` to include which side of screen target was on,
                                      TRUE ~ paste(CorrSide, Stim_Sides, sep = "_"))) %>%             # as well as whether it was valid with cue; if `iso_sides` parameter == `No`, leaves `StimSides` unchanged
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
                  list(~ifelse(block_acc <= input$block_floor | !between(miniblock_avg,  # Changes `Acc` and `RT` column values to NA if trial's block accuracy < `block_floor`
                                                                  input$miniblock_range[1], input$miniblock_range[2]), NA, .))) %>%    # or miniblock was not in the desired range
        mutate(Trials_filtered_out = sum(is.na(Acc)) / n(),
               Acc_postfilter = mean(Acc, na.rm = TRUE)) %>%                      # Creates column indicating mean accuracy for non-catch trials; note this is after before we've filtered for `block_floor`, unlike `Acc_prefilter`
        filter(between(Acc_postfilter, input$post_range[1], input$post_range[2])) %>%                # Filters out participants whose non-catch, post-block-filtering accuracy is outside of desired range
        group_by(CTI, Stim_Sides, !!!grouping_cnsts) %>%
        summarise_at(vars(Acc, RT), list(~mean(., na.rm = TRUE))) %>%              # Overwrites `Acc` and `RT` columns according to mean of each combination of `CTI` and `Stim_Sides`
        arrange(CTI) %>%
        group_by(Stim_Sides) %>%
        mutate_at(vars(Acc, RT), list(~na.approx(., na.rm = FALSE, rule = 2))) %>%
        mutate_at(vars(Acc, RT), list(~rollapply(., input$clumps, mean, partial = TRUE)))
    }
    
    cmbd <- do.call(rbind, lapply(pcpts, pcpts_combine)) %>%             # Calls `pcpts_combine` function for argumenet `pcpts`; then combines each participant's dataframe into one
      arrange(Acc_prefilter, participant, Stim_Sides, CTI)

    CTIs <- unique(cmbd$CTI)
    if (input$win_func == "tukey"){ win <- tukeywindow(length(CTIs), .5)} else { # Creates window, which if `tukey` will add the parameter `r` == `.5` —so 'only' half the data length will be non-flat
      win <- match.fun(paste0(input$win_func, "window"))(length(CTIs))}
    locations <- unique(cmbd$Stim_Sides)                                   # Creates vector of column names representing sides locations of target in reference to cue (and also potentially side of screen)
    pcpts <- unique(cmbd$participant)                                       # Creates vector of remaining participant numbers after `pcpts_combine` filtering
    cmbd_w <- cmbd %>%
      pivot_wide(CTI:CatchAcc, Stim_Sides, !!dep_var) %>%
      arrange(Acc_prefilter, participant, CTI) %>%
      group_by(participant) %>%
      mutate_at(vars(locations), list(~na.approx(., na.rm = FALSE, rule = 2)))
    
    # Analyzes Invalid - Valid instead of them separetely, if sbtr == `Yes`
    if (input$sbtr == "Yes"){
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
    amplitude <- function(x, z){
      pre_pad <- length(pcpts) * (length(unique(cmbd_w$CTI))) * z
      x %>%
        group_by(participant) %>%
        mutate_at(vars(locations), list(~detrend(.) * win)) %>%
        ungroup() %>%
        add_row(participant = rep(pcpts, times = z * (input$zeropads + 1))) %>%
        head(-length(pcpts) * z) %>%
        mutate_at(vars(locations), list(~coalesce(., 0))) %>%
        mutate(samp_shuff = ifelse(row_number() <= pre_pad, 
                                   RoundTo(row_number(), pre_pad / z,
                                           ceiling) / (pre_pad / z), 
                                   RoundTo(row_number() - pre_pad, (n() - pre_pad) / (z),
                                           ceiling) / ((n() - pre_pad) /  z))) %>%
        group_by(participant, samp_shuff) %>%
        mutate_at(vars(locations), list(~Mod(fft(.)))) %>%          # Detrends, multiplies by window, applies FFT, and then takes magnitude
        mutate(Hz = (row_number() - 1) / (n() * input$samp_per)) %>%
        ungroup() %>%
        select(-CTI)
    }
    
    amps <- amplitude(cmbd_w, 1)
    
    # Set Up Graphing
    
    t_srs_g <- function(x){
      cmbd_w %>%
        gather(Location, !!dep_var, -c(CTI, !!!grouping_cnsts)) %>%
        right_join(gather(conf_int(cmbd_w, CTI), Location, Conf_Int, -CTI),
                   by = c("CTI", "Location")) %>%
        group_by(CTI, Location, Conf_Int, !!!head(grouping_cnsts, x)) %>%
        summarise(!!dep_var := mean(!!dep_var)) %>%                               # Keeps either `RT` or `Acc` column—depending on whether `dep_var` parameter == `RT` or `Acc`
        ggplot(aes(CTI, !!dep_var, group = Location, color = Location,
                   fill = Location, ymin = !!dep_var - Conf_Int,
                   ymax = !!dep_var + Conf_Int)) + 
        labs(title = paste0(dep_var," by Cue-Target Interval, ", input$ext_objects, "-object task"),
             x = "Cue-Target Interval (ms)")
    }
    fft_x <- 1 / (length(unique(amps$Hz)) * input$samp_per)
    fft_g <- function(x) { x +
        labs(title = paste0("FFT of Target ", dep_var, ", ", input$ext_objects, "-object task"),
             col = "Target Location") +
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
                              labels = sapply(locations, function(x) gsub("_", " ", x),
                                              USE.NAMES = FALSE, simplify = TRUE)) +
        scale_fill_viridis_d(option = "C",
                             end = viridis_cols) +
        guides(colour = guide_legend(reverse = TRUE), fill = FALSE)
    }
    
    idvl_g <- function(y, x){
      graph(y, x) +
        theme(plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5)) +
        facet_wrap(~factor(participant, levels = pcpts),
                   ncol = RoundTo(length(pcpts) / 4, 1, ceiling),  scales = 'free_x')                                   # `free` means the y_axis isn't fixed from participant to participant
    }
    cmbd_g <- function(y, x) {
      graph(y, x) + geom_line(size = 1.5) +
        theme(legend.key.size = unit(.55, "in")) +
        labs(subtitle = paste( "Data from", as.character(length(pcpts)),
                               "participants"))
    }
    
    if (input$display == "FFT + Time-Series By Individual"){
      # Produces left half of final graph
      ts_facets <- idvl_g(t_srs_g, 1) +
        geom_line(alpha = I(2 / 10), color = "grey", show.legend = FALSE) +       # Graphs unsmoothed data in light gray
        stat_smooth(method = input$smooth_method, span = 0.2, se = FALSE,               # Smoothes data depending on `smooth_method` parameter
                    size = .5, show.legend = FALSE)
      
      # Produces label for each right-side graph
      plot_label <- amps %>%
        group_by(!!!grouping_cnsts) %>%
        summarise() %>%
        ungroup() %>%
        mutate_at(vars(!!!tail(grouping_cnsts, -1)),
                  list(~paste(quo_name(quo(.)), "=", percent(.)))) %>%
        unite(lab, !!!tail(grouping_cnsts, -1), sep = "\n", remove = FALSE)

      # Produces right half of final graph
      fft_facets <- idvl_g(fft_g, amps %>%
                             group_by(participant, Hz) %>%
                             summarise_all(mean) %>%
                             gather(Flash_and_or_field, Magnitude, -Hz, -samp_shuff, -c(!!!grouping_cnsts)) %>%
                             ggplot(aes(Hz, Magnitude, color = Flash_and_or_field))) +
        geom_line() +
        scale_x_continuous(name = "Frequency (Hz)", limits = c(0, input$xaxisvals)) +
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
      fft_x <- 1 / (length(unique(amps$Hz)) * samp_per)
      
      # Produces 'shuff' # of null hypothesis permutations
      shuffle <- function(x){
        cmbd_w %>%
          group_by(participant) %>%
          sample_n(length(CTIs), weight = CTI) %>%
          mutate_at(vars(CTI), funs(seq(min(CTIs), max(CTIs), samp_per))) %>%
          mutate(samp_shuff = x)
      }
      
      # Produces and save graph
      amps_shuff <- do.call(rbind, lapply(1:shuff, shuffle)) %>%
        amplitude(shuff) %>%
        group_by(Hz, samp_shuff) %>%
        summarise_at(vars(locations), mean) %>%
        group_by(Hz) %>%
        summarise_at(vars(locations), list(~quantile(., probs = 1 - pval))) %>%
        combine(amps %>% group_by(Hz) %>% summarise_at(vars(locations), mean),
                names = (c("Significance Cutoff", "Observed Data"))) %>%
        gather(Location, Magnitude, -c(Hz, source)) %>%
        right_join(gather(conf_int(amps, Hz), Location, Conf_Int, -Hz),
                   by = c("Hz", "Location"))
      output$mygraph <- renderPlot({
      (move_layers(cmbd_g(fft_g, ggplot(amps_shuff, aes(Hz, Magnitude, col = Location, linetype = source, 
                                                        ymin = Magnitude - Conf_Int, ymax = Magnitude + Conf_Int, fill = Location))) +
                     scale_linetype_manual(values = c("solid", "dashed")) +
                     scale_x_continuous(name = "Frequency (Hz)", limits = c(0, xaxisvals),
                                        breaks = seq(0, xaxisvals, ifelse(fft_x > .5, round(fft_x, 2), 1))) +
                     labs(linetype = "",
                          caption = paste0(input$zeropads, " zero pads added","\nSignificance threshold at p < ",
                                           as.character(pval))) +
                     geom_ribbon(data = filter(amps_shuff, source == "Observed Data"), alpha = 0.15, aes(color = NULL)) +
                     geom_point(size = 3, data = amps_shuff %>% spread(source, Magnitude) %>%
                                  filter(`Observed Data` > `Significance Cutoff`) %>%
                                  select(-c(`Significance Cutoff`, Conf_Int)) %>%
                                  gather(source, Magnitude, -Hz, -Location), 
                                aes(ymin = NULL, ymax = NULL)), "GeomRibbon", position = "bottom"))
      }, height = 800, width = 1350)
      
    }
    
  })

}
shinyApp(ui, server)
