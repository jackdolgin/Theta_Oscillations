# file.info("pre_registration_2_objects_no_wm.R")$ctime
# [1] "2020-01-20 23:41:35 EST"

if (!require(devtools)) install.packages("devtools")
if (!require(scienceverse)) devtools::install_github("scienceverse/scienceverse")
library("scienceverse")

source("FFT_analysis.R")

theta_2_objects_ext <- study(
  "Rhythms in Internal and External Attentional Focus: Replicating Oscillations in Two-Object, External Attention Task",
  author = c("Jack Dolgin", "Ian C. Fiebelkorn", "Tobias Egner")) %>%
  
  add_hypothesis("Accuracy during invalid trials will oscillate significantly only at 4 Hz., though accounting for multiple comparisons we will test frequencies 3-8 Hz. We will keep running participants until there are 30 usable, complete data sets.",
                 id = "hypothesis") %>%
  
  add_analysis(
    func = "main_function",
    params = list(
      averaging = "Aggregate",
      visualize = "Aggregate",
      batch = "4a",
      wm_exp = FALSE,
      iso_sides = FALSE,
      sbtr = FALSE,
      samp_per = 3/60,
      clumps = 0,
      dep_var = "Accuracy",
      α = .05,
      shuff = 50000,
      mult_correcs = "BH",
      trends = c("Detrending", "Demeaning"),
      smooth_method = "Loess",
      show_CI = TRUE,
      display = c("time_series", "FFT"),
      win_func = "Tukey",
      xrange = c(3, 8),
      duration = 1,
      attn_filter = FALSE,
      catch_floor = .85,
      side_bias = .5,
      wm_floor = .7,
      invalid_floor = .02,
      pre_range = c(.25, .85),
      post_range = c(.25, .85),
      filtered_cap = .5,
      block_range = c(.20, .80),
      blocks_desired = 1:8,
      miniblock_range = c(.20, .80),
      CTI_range = c(.3, 1.09)                                 
    ),
    id = "analyses") %>%
  
  study_analyze()

study_save(theta_2_objects_ext, "theta_2_objects_ext_preregistration.json")
