# > file.info("pre_registration_2_and_3_objects.R")$ctime
# [1] "2019-06-23 15:29:54 EDT

if (!require(devtools)) install.packages("devtools")
if (!require(scienceverse)) devtools::install_github("scienceverse/scienceverse")
library("scienceverse")

theta_1_2_study <- study(
  "Rhythms in Internal and External Attentional Focus, Exps 1 + 2",
  author = c("Jack Dolgin", "Ian C. Fiebelkorn", "Tobias Egner")) %>%
  
  add_hypothesis("In the 2-object task, detection rates at invalidly-cued locations will peak at 4 Hz., and this peak will reach statistical significance.",
                 id = "h1") %>%
  
  add_hypothesis("In the 2-object task, response times at invalidly-cued locations will peak at 4 Hz., and this peak will reach statistical significance.",
                 id = "h2") %>%
  
  add_hypothesis("In the 3-object task, detection rates at invalidly-cued locations will peak at a lower frequency—2 or 3 Hz.—than in the 2-object task, and this peak will reach statistical significance.",
  id = "h3") %>%
  
  add_hypothesis("In the 3-object task, response times at invalidly-cued locations will peak at a lower frequency—2 or 3 Hz.—than in the 2-object task, and this peak will reach statistical significance.",
                 id = "h4") %>%
  
  add_analysis(
    func = "main_function",
    params = list(
      dset = "Experimental",
      display = "FFT Across Participants",
      ext_objects = 2,
      iso_sides = FALSE,
      sbtr = FALSE,
      samp_per = round(1 / 60, 4),
      clumps = 2,
      dep_var = "Accuracy",
      pval = .025,
      shuff = 50,
      trends = c("Detrending", "Demeaning"),
      smooth_method = "Loess",
      win_func = "Tukey",
      xaxisvals = 15,
      duration = 1,
      attn_filter = FALSE,
      catch_floor = .85,
      side_bias = .2,
      pre_range = c(.45, .75),
      post_range = c(.45, .75),
      block_range = c(.40, .80),
      miniblock_range = c(.40, .80),
      CTI_range = c(.3, 1.29)
    ),
    id = "main_analysis_h1") %>%
  
  add_analysis(
    func = "main_function",
    params = list(
      dset = "Experimental",
      display = "FFT Across Participants",
      ext_objects = 2,
      iso_sides = FALSE,
      sbtr = FALSE,
      samp_per = round(1 / 60, 4),
      clumps = 2,
      dep_var = "Response Time",
      pval = .025,
      shuff = 50,
      trends = c("Detrending", "Demeaning"),
      smooth_method = "Loess",
      win_func = "Tukey",
      xaxisvals = 15,
      duration = 1,
      attn_filter = FALSE,
      catch_floor = .85,
      side_bias = .2,
      pre_range = c(.45, .75),
      post_range = c(.45, .75),
      block_range = c(.40, .80),
      miniblock_range = c(.40, .80),
      CTI_range = c(.3, 1.29)
    ),
    id = "main_analysis_h2") %>%
  
  add_analysis(
    func = "main_function",
    params = list(
      dset = "Experimental",
      display = "FFT Across Participants",
      ext_objects = 3,
      iso_sides = FALSE,
      sbtr = FALSE,
      samp_per = round(1 / 60, 4),
      clumps = 2,
      dep_var = "Accuracy",
      pval = .025,
      shuff = 50,
      trends = c("Detrending", "Demeaning"),
      smooth_method = "Loess",
      win_func = "Tukey",
      xaxisvals = 15,
      duration = 1,
      attn_filter = FALSE,
      catch_floor = .85,
      side_bias = .2,
      pre_range = c(.45, .75),
      post_range = c(.45, .75),
      block_range = c(.40, .80),
      miniblock_range = c(.40, .80),
      CTI_range = c(.3, 1.29)
    ),
    id = "main_analysis_h3") %>%
  
  add_analysis(
    func = "main_function",
    params = list(
      dset = "Experimental",
      display = "FFT Across Participants",
      ext_objects = 3,
      iso_sides = FALSE,
      sbtr = FALSE,
      samp_per = round(1 / 60, 4),
      clumps = 2,
      dep_var = "Response Time",
      pval = .025,
      shuff = 50,
      trends = c("Detrending", "Demeaning"),
      smooth_method = "Loess",
      win_func = "Tukey",
      xaxisvals = 15,
      duration = 1,
      attn_filter = FALSE,
      catch_floor = .85,
      side_bias = .2,
      pre_range = c(.45, .75),
      post_range = c(.45, .75),
      block_range = c(.40, .80),
      miniblock_range = c(.40, .80),
      CTI_range = c(.3, 1.29)
    ),
    id = "main_analysis_h4") %>%
  
  study_analyze()
  
study_save(theta_1_2_study, "pre_data_theta_1_2_study.json")
