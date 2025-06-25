context("Test plot-prism")

# prep data ----
mae <- gDRutils::get_synthetic_data("combo_matrix")
se_sa <- mae[[gDRutils::get_supported_experiments("sa")]]
dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se_sa,
                                               assay_name = "Metrics")
dt_average <- gDRutils::convert_se_assay_to_dt(se = se_sa,
                                               assay_name = "Averaged")
dt_metrics_capped <-
  gDRutils::cap_assay_infinities(conc_assay_dt = dt_average,
                                 assay_dt = dt_metrics,
                                 experiment_name = gDRutils::get_supported_experiments("sa"),
                                 capping_fold = 5)

se_combo <- mae[[gDRutils::get_supported_experiments("combo")]]
dt_metrics_combo <- gDRutils::convert_se_assay_to_dt(se = se_combo,
                                                     assay_name = "Metrics")
dt_scores <- gDRutils::convert_se_assay_to_dt(se = se_combo,
                                              assay_name = "scores")

d_name <- "drug_002"
d_name2 <- "drug_026"
dt_response_met <- 
  prep_dt_response_metric_sa(dt_metrics, d_name,
                             metric = c("xc50", "x_mean", "x_max"))
dt_response_met_capped <- 
  prep_dt_response_metric_sa(dt_metrics_capped, d_name,
                             metric = c("xc50", "x_mean", "x_max"))
dt_response_dose <- 
  prep_dt_response_dose_sa(dt_average, d_name)

dt_response_score <- 
  prep_dt_response_scores(dt_scores, d_name, d_name2, 
                          metric = c("hsa_score", "bliss_score"))
dt_response_diff <-
  prep_dt_response_metric_diff(dt_metrics_combo, d_name, d_name2,
                               metric = c("xc50", "x_mean", "x_max"))

#_assoc
obj_depmap_feat <- 
  prep_dt_depmap_feat(feat_data_path = system.file("testdata", package = "gDRplots"),
                      meta_data_path = system.file("testdata/Model.csv", package = "gDRplots"))
obj_depmap_feat_2 <- 
  prep_dt_depmap_feat(feat_data_path = system.file("testdata", package = "gDRplots"),
                      meta_data_path = system.file("testdata/Model.csv", package = "gDRplots"),
                      feature_set = "OmicsSomaticMutationsMatrixHotspot")

ls_feat <- setdiff(names(obj_depmap_feat[["dt_depmap"]]), c("ModelID", "CCLEName"))
ls_feat_2 <- setdiff(names(obj_depmap_feat_2[["dt_depmap"]]), c("ModelID", "CCLEName"))
dist_q <- withr::with_seed(42, rnorm(400, mean = 0.3, sd = 0.05))

dt_assoc_sa <- data.table::data.table(
  feature = ls_feat,
  response = rep("RV_gDR_xc50", NROW(ls_feat)),
  rho = withr::with_seed(42, rnorm(NROW(ls_feat), mean = 0, sd = 0.05)),
  q_value = withr::with_seed(42, sample(dist_q[dist_q < 0.3], NROW(ls_feat), replace = TRUE))
)

obj_assoc_sa <- list(dt_assoc = dt_assoc_sa,
                     condition_info = unique(dt_response_met[["rId"]]),
                     selected_metric = "RV_gDR_xc50",
                     selected_feat_meta_col = obj_depmap_feat[["selected_feat_meta_col"]])

obj_depmap_meta <- 
  prep_dt_depmap_meta(meta_data_path = system.file("testdata/Model.csv", package = "gDRplots"))

ls_meta <- setdiff(names(obj_depmap_meta[["dt_depmap"]]), c("ModelID", "CCLEName"))
dist_q <- withr::with_seed(42, rnorm(400, mean = 1, sd = 0.1))

dt_assoc_combo <- data.table::data.table(
  feature = ls_meta,
  response = rep("hsa_score", NROW(ls_meta)),
  rho = withr::with_seed(42, rnorm(NROW(ls_meta), mean = 0, sd = 0.05)),
  q_value = withr::with_seed(42, sample(dist_q[dist_q < 0.75], NROW(ls_meta), replace = TRUE))
)

obj_assoc_combo <- list(dt_assoc = dt_assoc_combo,
                        condition_info = unique(dt_response_score[["rId"]]),
                        selected_metric = "hsa_score",
                        selected_feat_meta_col = obj_depmap_meta[["selected_feat_meta_col"]])

# tests ----
test_that("plot_volcano_assoc works as expected", {
  plt_1 <- plot_volcano_assoc(dt_assoc = obj_assoc_sa[["dt_assoc"]],
                              selected_feat_meta_col = obj_assoc_sa[["selected_feat_meta_col"]],
                              selected_metric = obj_assoc_sa[["selected_metric"]]) # default
  
  expect_is(plt_1, "gg")
  expect_length(plt_1[["layers"]], 2)
  expect_equal(plt_1[["labels"]][["x"]], "rho") # predef for x axis
  expect_equal(plt_1[["labels"]][["y"]], "neglog_q_value") # predef for y axis
  expect_equal(plt_1[["labels"]][["title"]], 
               paste0("RV_gDR_xc50__", obj_assoc_sa[["selected_feat_meta_col"]])) # <metric>__<feat>
  
  # scenario: check number of label to be shown
  no_lbl <- 3
  q_alpha <- 0.25
  plt_2 <- plot_volcano_assoc(dt_assoc = obj_assoc_sa[["dt_assoc"]],
                              selected_metric = obj_assoc_sa[["selected_metric"]],
                              selected_feat_meta_col = obj_assoc_sa[["selected_feat_meta_col"]],
                              condition_info = obj_assoc_sa[["condition_info"]],
                              alpha = q_alpha,
                              named_p_top = no_lbl)
  expect_is(plt_2, "gg")
  expect_equal(plt_2[["labels"]][["subtitle"]], obj_assoc_sa[["condition_info"]])
  expect_length(unique(unlist(ggplot2::ggplot_build(plt_2)$data[[1]]$colour)), 
                NROW(unique(obj_assoc_sa[["dt_assoc"]]$q_value <= q_alpha))) # stat significant
  expect_length( 
    ggplot2::ggplot_build(plt_2)$data[[1]]$label[ggplot2::ggplot_build(plt_2)$data[[1]]$label != ""],
    no_lbl) # lbl for top feat
  
  plt_3 <- plot_volcano_assoc(dt_assoc = obj_assoc_combo[["dt_assoc"]],
                              selected_feat_meta_col = obj_assoc_combo[["selected_feat_meta_col"]],
                              selected_metric = obj_assoc_combo[["selected_metric"]])
  expect_is(plt_3, "gg")
  expect_length(plt_3[["layers"]], 2)
  expect_equal(plt_3[["labels"]][["x"]], "rho") # predef for x axis
  expect_equal(plt_3[["labels"]][["y"]], "neglog_q_value") # predef for y axis
  expect_equal(plt_3[["labels"]][["title"]], 
               paste0("hsa_score__", obj_assoc_combo[["selected_feat_meta_col"]])) # <metric>__<meta>
  
  # scenario: check high alpha
  q_alpha_2 <- 0.71
  plt_4 <- plot_volcano_assoc(dt_assoc = obj_assoc_combo[["dt_assoc"]],
                              selected_feat_meta_col = obj_assoc_combo[["selected_feat_meta_col"]],
                              selected_metric = obj_assoc_combo[["selected_metric"]],
                              alpha = q_alpha_2)
  expect_is(plt_4, "gg")
  plt_4_data <- data.table::as.data.table(ggplot2::ggplot_build(plt_4)$data[[1]])
  expect_length(unique(unlist(plt_4_data$colour)), 
                NROW(unique(obj_assoc_combo[["dt_assoc"]]$q_value <= q_alpha_2))) # stat significant
  expect_equal(sort(plt_4_data$label[which(unlist(plt_4_data$colour) == "black")]),
               sort(obj_assoc_combo[["dt_assoc"]][q_value <= q_alpha_2, ]$feature))
  
  # scenario: NAs in depmap
  dt_assoc_na <- data.table::copy(obj_assoc_sa[["dt_assoc"]])
  plt_5 <- plot_volcano_assoc(dt_assoc = dt_assoc_na[0, ],
                              selected_feat_meta_col = obj_assoc_sa[["selected_feat_meta_col"]],
                              selected_metric = obj_assoc_sa[["selected_metric"]])
  expect_is(plt_5, "gg")
  expect_length(plt_5[["layers"]], 0) # empty plot
  expect_equal(plt_5[["labels"]][["x"]], "rho") # predef for x axis
  expect_equal(plt_5[["labels"]][["y"]], "neglog_q_value") # predef for y axis
  expect_true(grepl(": all NAs", plt_5[["labels"]][["title"]]))
  
  dt_assoc_na[["q_value"]] <- NA
  plt_6 <- plot_volcano_assoc(dt_assoc = dt_assoc_na,
                              selected_feat_meta_col = obj_assoc_sa[["selected_feat_meta_col"]],
                              selected_metric = obj_assoc_sa[["selected_metric"]])
  expect_is(plt_6, "gg")
  expect_length(plt_6[["layers"]], 0) # empty plot
  expect_equal(plt_6[["labels"]][["x"]], "rho") # predef for x axis
  expect_equal(plt_6[["labels"]][["y"]], "neglog_q_value") # predef for y axis
  expect_true(grepl(": all NAs", plt_6[["labels"]][["title"]]))
  
  # scenario: check max_N
  dt_assoc_big <- data.table::data.table(
    feature = sprintf("XZ_A%02dTT", 1:50),
    response = rep("RV_gDR_xc50", 50),
    rho = withr::with_seed(42, rnorm(n = 50, mean = 0, sd = 0.035)),
    q_value = withr::with_seed(42, rnorm(n = 50, mean = 0.15, sd = 0.05))
  )
  max_non_stat_sig <- 20
  plt_7 <- plot_volcano_assoc(dt_assoc = dt_assoc_big,
                              selected_metric = "RV_gDR_xc50",
                              selected_feat_meta_col = "XZ_fatures",
                              condition_info = NULL,
                              max_N = max_non_stat_sig)
  expect_is(plt_7, "gg")
  expect_length(plt_7[["layers"]], 2)
  expect_equal(sum(ggplot2::ggplot_build(plt_7)$data[[1]][["label"]] != ""), 10) # default named_p_top
  expect_equal(NROW(ggplot2::ggplot_build(plt_7)$data[[1]]),
               NROW(dt_assoc_big[q_value < 0.05]) + max_non_stat_sig) # default alpha
  
  # scenario: dt_assoc has more than 1 values in `response` and `selected_metric` is not one of them
  dt_assoc_mix <- data.table::data.table(
    feature = sprintf("XZ_A%02dTT", 1:50),
    response = rep(c("RV_gDR_hsa_score", "RV_gDR_bliss_score"), length.out = 50),
    rho = withr::with_seed(42, rnorm(n = 50, mean = 0, sd = 0.035)),
    q_value = withr::with_seed(42, rnorm(n = 50, mean = 0.15, sd = 0.05))
  )
  expect_warning({
    plt_8 <- plot_volcano_assoc(dt_assoc = dt_assoc_mix,
                                selected_metric = "RV_gDR_xc50",
                                selected_feat_meta_col = "XZ_fatures",
                                condition_info = NULL)
  }, "Association data is not consistent - there is more than one value in the `response` column.")
  expect_is(plt_8, "gg")
  expect_length(plt_8[["layers"]], 2)
  expect_equal(NROW(ggplot2::ggplot_build(plt_8)$data[[1]]), NROW(dt_assoc_mix))
  
  # scenario: dt_assoc has more than 1 values in `response` and `selected_metric` is one of them
  sel_met <- "RV_gDR_hsa_score"
  expect_warning({
    plt_9 <- plot_volcano_assoc(dt_assoc = dt_assoc_mix,
                                selected_metric = sel_met,
                                selected_feat_meta_col = "XZ_fatures",
                                condition_info = NULL)
  }, "Association data was filtered based on `selected_metric`")
  expect_is(plt_9, "gg")
  expect_length(plt_9[["layers"]], 2)
  expect_equal(NROW(ggplot2::ggplot_build(plt_9)$data[[1]]), NROW(dt_assoc_mix[response == sel_met, ]))
  
  # scenario: dt_assoc does not have column `response`
  plt_10 <- plot_volcano_assoc(dt_assoc = obj_assoc_sa[["dt_assoc"]][, -c("response")],
                               selected_feat_meta_col = obj_assoc_sa[["selected_feat_meta_col"]],
                               selected_metric = obj_assoc_sa[["selected_metric"]]) # default
  expect_is(plt_10, "gg")
  expect_equal(NROW(plt_10[["layers"]]), NROW(plt_1[["layers"]]))
  expect_equal(plt_10[["labels"]][["x"]], plt_1[["labels"]][["x"]]) # predef for x axis
  expect_equal(plt_10[["labels"]][["y"]], plt_1[["labels"]][["y"]]) # predef for y axis
  expect_equal(plt_10[["labels"]][["title"]], plt_1[["labels"]][["title"]]) # <metric>__<feat>
  
  # testing assertions
  expect_error(plot_volcano_assoc(dt_assoc = obj_assoc_sa,
                                  selected_feat_meta_col = obj_assoc_sa[["selected_feat_meta_col"]]),
               "Assertion on 'dt_assoc' failed: Must be a data.table")
  expect_error(plot_volcano_assoc(dt_assoc = obj_assoc_sa[["dt_assoc"]],
                                  selected_feat_meta_col = obj_assoc_sa),
               "Assertion on 'selected_feat_meta_col' failed: Must be of type 'string'")
  expect_error(plot_volcano_assoc(dt_assoc = obj_assoc_sa[["dt_assoc"]],
                                  selected_feat_meta_col = obj_assoc_sa[["selected_feat_meta_col"]],
                                  selected_metric = 1),
               "Assertion on 'selected_metric' failed: Must be of type 'string'")
  expect_error(plot_volcano_assoc(dt_assoc = obj_assoc_sa[["dt_assoc"]],
                                  selected_feat_meta_col = obj_assoc_sa[["selected_feat_meta_col"]],
                                  selected_metric = obj_assoc_sa[["selected_metric"]],
                                  condition_info = 123),
               "Assertion on 'condition_info' failed: Must be of type 'string'")
  expect_error(plot_volcano_assoc(dt_assoc = obj_assoc_sa[["dt_assoc"]],
                                  selected_feat_meta_col = obj_assoc_sa[["selected_feat_meta_col"]],
                                  selected_metric = obj_assoc_sa[["selected_metric"]],
                                  alpha = "0.1"),
               "Assertion on 'alpha' failed: Must be of type 'number'")
  expect_error(plot_volcano_assoc(dt_assoc = obj_assoc_sa[["dt_assoc"]],
                                  selected_feat_meta_col = obj_assoc_sa[["selected_feat_meta_col"]],
                                  selected_metric = obj_assoc_sa[["selected_metric"]],
                                  named_p_top = "5"),
               "Assertion on 'named_p_top' failed: Must be of type 'number'")
  expect_error(plot_volcano_assoc(dt_assoc = obj_assoc_sa[["dt_assoc"]],
                                  selected_feat_meta_col = obj_assoc_sa[["selected_feat_meta_col"]],
                                  selected_metric = obj_assoc_sa[["selected_metric"]],
                                  max_N = "only 10"),
               "Assertion on 'max_N' failed: Must be of type 'number'")
  expect_error(plot_volcano_assoc(dt_assoc = obj_assoc_sa[["dt_assoc"]],
                                  selected_feat_meta_col = obj_assoc_sa[["selected_feat_meta_col"]],
                                  selected_metric = obj_assoc_sa[["selected_metric"]],
                                  max_N = 2:6),
               "Assertion on 'max_N' failed: Must have length 1.")
  expect_error(plot_volcano_assoc(dt_assoc = obj_assoc_sa[["dt_assoc"]],
                                  selected_feat_meta_col = obj_assoc_sa[["selected_feat_meta_col"]],
                                  selected_metric = obj_assoc_sa[["selected_metric"]],
                                  max_N = 6),
               "Assertion on 'max_N' failed: Element 1 is not >= 10.")
})

test_that("plot_scatter_with_corr works as expected", {
  selected_feat <- ls_feat[1]
  selected_metric <- "RV_gDR_x_0.01"
  dt_response <- dt_response_dose[, c("rId", "cId", "CellLineName", selected_metric), with = FALSE]
  
  plt_1 <- plot_scatter_with_corr(dt_response = dt_response,
                                  dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                  selected_feat = selected_feat)
  expect_is(plt_1, "gg")
  expect_length(plt_1[["layers"]], 3)
  expect_equal(plt_1[["labels"]][["x"]], selected_feat)
  expect_equal(plt_1[["labels"]][["y"]], selected_metric)
  expect_equal(plt_1[["labels"]][["title"]], NULL)
  expect_equal(plt_1[["labels"]][["caption"]], unique(dt_response$rId))
  expect_true(all(vapply(c("corr", "slope", "intercept"), 
                         function(i) grepl(i, plt_1[["labels"]][["subtitle"]]), logical(1))))
  
  plt_2 <- plot_scatter_with_corr(dt_response = dt_response,
                                  dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                  selected_feat = selected_feat,
                                  selected_feat_meta_col = obj_depmap_feat[["selected_feat_meta_col"]])
  expect_is(plt_2, "gg")
  expect_length(plt_2[["layers"]], 3)
  expect_equal(plt_2[["labels"]][["title"]], obj_depmap_feat[["selected_feat_meta_col"]])
  
  selected_feat_2 <- "XZ_A5BN"
  selected_metric_2 <- "RV_gDR_bliss_score"
  dt_response_2 <- dt_response_score[, c("rId", "cId", "CellLineName", selected_metric_2), with = FALSE]
  
  plt_3 <- plot_scatter_with_corr(dt_response = dt_response_2,
                                  dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                  selected_feat = selected_feat_2,
                                  selected_feat_meta_col = obj_depmap_feat[["selected_feat_meta_col"]])
  expect_is(plt_3, "gg")
  expect_length(plt_3[["layers"]], 3)
  expect_equal(plt_3[["labels"]][["x"]], selected_feat_2)
  expect_equal(plt_3[["labels"]][["y"]], selected_metric_2)
  expect_equal(plt_3[["labels"]][["title"]], obj_depmap_feat[["selected_feat_meta_col"]])
  expect_equal(plt_3[["labels"]][["caption"]], unique(dt_response_2$rId))
  
  # NAs in response
  dt_response_na <- data.table::copy(dt_response)
  dt_response_na[[selected_metric]] <- NA
  plt_4 <- plot_scatter_with_corr(dt_response = dt_response_na,
                                  dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                  selected_feat = selected_feat,
                                  selected_feat_meta_col = obj_depmap_feat[["selected_feat_meta_col"]])
  expect_is(plt_4, "gg")
  expect_length(plt_4[["layers"]], 0) # empty plot
  expect_equal(plt_4[["labels"]][["x"]], selected_feat)
  expect_equal(plt_4[["labels"]][["y"]], selected_metric)
  expect_equal(plt_4[["labels"]][["title"]], 
               paste0(obj_depmap_feat[["selected_feat_meta_col"]], ": all NAs"))
  expect_equal(plt_4[["labels"]][["caption"]], unique(dt_response_na$rId))
  
  # NAs in depmap
  dt_depmap_na <- data.table::copy(obj_depmap_feat[["dt_depmap"]])
  dt_depmap_na[[selected_feat]] <- NA
  plt_5 <- plot_scatter_with_corr(dt_response = dt_response,
                                  dt_depmap = dt_depmap_na, 
                                  selected_feat = selected_feat,
                                  selected_feat_meta_col = obj_depmap_feat[["selected_feat_meta_col"]])
  expect_is(plt_5, "gg")
  expect_length(plt_5[["layers"]], 0) # empty plot
  expect_equal(plt_5[["labels"]][["x"]], selected_feat)
  expect_equal(plt_5[["labels"]][["y"]], selected_metric)
  expect_equal(plt_5[["labels"]][["title"]], 
               paste0(obj_depmap_feat[["selected_feat_meta_col"]], ": all NAs"))
  expect_equal(plt_5[["labels"]][["caption"]], unique(dt_response_na$rId))
  
  # testing assertions
  expect_error(plot_scatter_with_corr(dt_response = unlist(dt_response),
                                      dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                      selected_feat = selected_feat),
               "Assertion on 'dt_response' failed: Must be a data.table")
  expect_error(plot_scatter_with_corr(dt_response = dt_response,
                                      dt_depmap = obj_depmap_feat, 
                                      selected_feat = selected_feat),
               "Assertion on 'dt_depmap' failed: Must be a data.table")
  expect_error(plot_scatter_with_corr(dt_response = dt_response,
                                      dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                      selected_feat = 1),
               "Assertion on 'selected_feat' failed: Must be of type 'string'")
  expect_error(plot_scatter_with_corr(dt_response = dt_response,
                                      dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                      selected_feat = "not_existen_feat"),
               "Assertion on 'names\\(dt_depmap\\)' failed: Names must include the elements")
  expect_error(plot_scatter_with_corr(dt_response = dt_response,
                                      dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                      selected_feat = selected_feat,
                                      selected_feat_meta_col = 1),
               "Assertion on 'selected_feat_meta_col' failed: Must be of type 'string'")
}) 

test_that("plot_scatter_with_corr_panel works as expected", {
  selected_feats <- ls_feat[seq_len(3)]
  selected_metric <- "RV_gDR_x_10"
  dt_response <- dt_response_dose[, c("rId", "cId", "CellLineName", selected_metric), with = FALSE]
  
  plt_1 <- plot_scatter_with_corr_panel(dt_response = dt_response,
                                        dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                        selected_feats = selected_feats)
  expect_is(plt_1, "gg")
  expect_length(plt_1[["layers"]], 3)
  expect_equal(plt_1[["labels"]][["x"]], "")
  expect_equal(plt_1[["labels"]][["y"]], selected_metric)
  expect_equal(plt_1[["labels"]][["title"]], NULL)
  expect_equal(plt_1[["labels"]][["caption"]], unique(dt_response$rId))
  
  plt_2 <- plot_scatter_with_corr_panel(dt_response = dt_response,
                                        dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                        selected_feats = selected_feats[1],
                                        selected_feat_meta_col = obj_depmap_feat[["selected_feat_meta_col"]])
  expect_is(plt_2, "gg")
  expect_length(plt_2[["layers"]], 3)
  expect_equal(plt_2[["labels"]][["title"]], obj_depmap_feat[["selected_feat_meta_col"]])
  expect_equal(plt_2[["labels"]][["caption"]], unique(dt_response$rId))
  
  # selected feat is not present in dt_depmap
  new_selected_feats <- c(selected_feats, "XZ_non_avial")
  plt_3 <- plot_scatter_with_corr_panel(dt_response = dt_response,
                                        dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                        selected_feats = new_selected_feats)
  expect_is(plt_3, "gg")
  expect_length(plt_3[["layers"]], 3)
  expect_equal(plt_3[["labels"]][["y"]], selected_metric)
  expect_equal(
    NROW(data.table::data.table(ggplot2::ggplot_build(plt_3)[["data"]][[1]])[alpha == 0, .N, by = PANEL]), # nolint 
    NROW(new_selected_feats[!new_selected_feats %in% names(obj_depmap_feat[["dt_depmap"]])])) 
  
  # only NAs in selected_feats 
  plt_4 <- 
    plot_scatter_with_corr_panel(dt_response = dt_response,
                                 dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                 selected_feats = rep(NA, 2),
                                 selected_feat_meta_col = obj_depmap_feat[["selected_feat_meta_col"]])
  expect_is(plt_4, "gg")
  expect_length(plt_4[["layers"]], 0) # empty plot
  expect_equal(plt_4[["labels"]][["x"]], "")
  expect_equal(plt_4[["labels"]][["y"]], selected_metric)
  expect_equal(plt_4[["labels"]][["title"]], 
               paste0(obj_depmap_feat[["selected_feat_meta_col"]], ": all NAs"))
  
  # some NAs in selected_feats 
  selected_feats_with_NAs <- c(NA, selected_feats, NA) 
  plt_5 <- plot_scatter_with_corr_panel(dt_response = dt_response,
                                        dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                        selected_feats = selected_feats_with_NAs)
  expect_is(plt_5, "gg")
  expect_length(plt_5[["layers"]], 3)
  expect_equal(plt_5[["labels"]][["y"]], selected_metric)
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_5)[["data"]][[1]]$PANEL)),
               NROW(selected_feats_with_NAs))
  expect_equal(
    NROW(data.table::data.table(ggplot2::ggplot_build(plt_5)[["data"]][[1]])[alpha == 0, .N, by = PANEL]), # nolint 
    NROW(selected_feats_with_NAs[!selected_feats_with_NAs %in% names(obj_depmap_feat[["dt_depmap"]])])) 
  
  # NAs in response
  dt_response_na <- data.table::copy(dt_response)
  dt_response_na[[selected_metric]] <- NA
  plt_6 <- plot_scatter_with_corr_panel(dt_response = dt_response_na,
                                        dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                        selected_feats = selected_feats)
  expect_is(plt_6, "gg")
  expect_equal(plt_6[["labels"]][["x"]], "")
  expect_equal(plt_6[["labels"]][["y"]], selected_metric)
  expect_equal(plt_6[["labels"]][["title"]], NULL)
  expect_equal(plt_6[["labels"]][["caption"]], unique(dt_response$rId))
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_6)[["data"]][[1]]$PANEL)), 
               NROW(selected_feats))
  expect_equal(NROW(ggplot2::ggplot_build(plt_6)[["data"]]), NROW(selected_feats))
  
  # NAs in depmap
  dt_depmap_na <- data.table::copy(obj_depmap_feat[["dt_depmap"]])
  dt_depmap_na[[selected_feats[2]]] <- NA
  plt_7 <- plot_scatter_with_corr_panel(dt_response = dt_response,
                                        dt_depmap = dt_depmap_na, 
                                        selected_feats = selected_feats)
  expect_is(plt_7, "gg")
  expect_equal(plt_7[["labels"]][["x"]], "")
  expect_equal(plt_7[["labels"]][["y"]], selected_metric)
  expect_equal(plt_7[["labels"]][["title"]], NULL)
  expect_equal(plt_7[["labels"]][["caption"]], unique(dt_response$rId))
  expect_equal(
    NROW(data.table::data.table(ggplot2::ggplot_build(plt_7)[["data"]][[1]])[alpha == 0, .N, by = PANEL]), # nolint 
    sum(vapply(selected_feats, function(nm) all(is.na(dt_depmap_na[[nm]])), logical(1)))) 
  
  # testing assertions
  expect_error(plot_scatter_with_corr_panel(dt_response = unlist(dt_response),
                                            dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                            selected_feats = selected_feats),
               "Assertion on 'dt_response' failed: Must be a data.table")
  expect_error(plot_scatter_with_corr_panel(dt_response = dt_response,
                                            dt_depmap = obj_depmap_feat, 
                                            selected_feats = selected_feats),
               "Assertion on 'dt_depmap' failed: Must be a data.table")
  expect_error(plot_scatter_with_corr_panel(dt_response = dt_response,
                                            dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                            selected_feats = 1),
               "Assertion on 'selected_feats' failed: Must be of type 'character'")
  expect_error(plot_scatter_with_corr_panel(dt_response = dt_response,
                                            dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                            selected_feats = selected_feats,
                                            selected_feat_meta_col = 1),
               "Assertion on 'selected_feat_meta_col' failed: Must be of type 'string'")
})

test_that("plot_boxplot_num works as expected", {
  selected_feat <- ls_feat_2[1]
  selected_metric <- "RV_gDR_x_0.01"
  dt_response <- dt_response_dose[, c("rId", "cId", "CellLineName", selected_metric), with = FALSE]
  
  plt_1 <- plot_boxplot_num(dt_response = dt_response,
                            dt_depmap = obj_depmap_feat_2[["dt_depmap"]],
                            selected_feat = selected_feat) # default
  expect_is(plt_1, "gg")
  expect_length(plt_1[["layers"]], 3)
  expect_equal(plt_1[["labels"]][["x"]], selected_feat)
  expect_equal(plt_1[["labels"]][["y"]], selected_metric)
  expect_equal(plt_1[["labels"]][["title"]], NULL)
  expect_equal(plt_1[["labels"]][["caption"]], unique(dt_response$rId))
  expect_equal(plt_1[["labels"]][["caption"]], unique(dt_response$rId))
  expect_equal(# check the uniqueness of points
    NROW(ggplot2::ggplot_build(plt_1)$data[[3]]), 
    NROW(obj_depmap_feat_2[["dt_depmap"]][!is.na(get(selected_feat)) & CCLEName %in% dt_response[["CellLineName"]], ]))
  
  # scenario: only one level
  selected_feat_2 <- ls_feat_2[5]
  selected_metric_2 <- "RV_gDR_bliss_score"
  dt_response_2 <- dt_response_score[, c("rId", "cId", "CellLineName", selected_metric_2), with = FALSE]
  
  plt_2 <- plot_boxplot_num(dt_response = dt_response_2,
                            dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                            selected_feat = selected_feat_2,
                            selected_feat_meta_col = obj_depmap_feat_2[["selected_feat_meta_col"]])
  expect_is(plt_2, "gg")
  expect_length(plt_2[["layers"]], 3)
  expect_equal(plt_2[["labels"]][["x"]], selected_feat_2)
  expect_equal(plt_2[["labels"]][["y"]], selected_metric_2)
  expect_equal(plt_2[["labels"]][["title"]], obj_depmap_feat_2[["selected_feat_meta_col"]])
  expect_equal(plt_2[["labels"]][["caption"]], unique(dt_response_2$rId))
  
  res_count_3 <- obj_depmap_feat_2[["dt_depmap"]][CCLEName %in% dt_response_2$CellLineName]
  res_count_3 <-
    res_count_3[!is.na(get(selected_feat)), .N, by = selected_feat][, lbl := sprintf("%s (%s)", get(selected_feat), N)]
  plt_3 <- plot_boxplot_num(dt_response = dt_response_2,
                            dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                            selected_feat = selected_feat,
                            selected_feat_meta_col = obj_depmap_feat_2[["selected_feat_meta_col"]])
  expect_is(plt_3, "gg")
  expect_length(plt_3[["layers"]], 3)
  expect_equal(plt_3[["labels"]][["x"]], selected_feat)
  expect_equal(sort(ggplot2::layer_scales(plt_3)$x$get_labels()), sort(res_count_3$lbl))
  
  # NAs in response
  dt_response_na <- data.table::copy(dt_response)
  dt_response_na[[selected_metric]] <- NA
  plt_4 <- plot_boxplot_num(dt_response = dt_response_na,
                            dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                            selected_feat = selected_feat,
                            selected_feat_meta_col = obj_depmap_feat_2[["selected_feat_meta_col"]])
  expect_is(plt_4, "gg")
  expect_length(plt_4[["layers"]], 0) # empty plot
  expect_equal(plt_4[["labels"]][["x"]], selected_feat)
  expect_equal(plt_4[["labels"]][["y"]], selected_metric)
  expect_equal(plt_4[["labels"]][["title"]], 
               paste0(obj_depmap_feat_2[["selected_feat_meta_col"]], ": all NAs"))
  expect_equal(plt_4[["labels"]][["caption"]], unique(dt_response_na$rId))
  
  # NAs in depmap
  dt_depmap_na <- data.table::copy(obj_depmap_feat_2[["dt_depmap"]])
  dt_depmap_na[[selected_feat]] <- NA
  plt_5 <- plot_boxplot_num(dt_response = dt_response,
                            dt_depmap = dt_depmap_na, 
                            selected_feat = selected_feat,
                            selected_feat_meta_col = obj_depmap_feat_2[["selected_feat_meta_col"]])
  expect_is(plt_5, "gg")
  expect_length(plt_5[["layers"]], 0) # empty plot
  expect_equal(plt_5[["labels"]][["x"]], selected_feat)
  expect_equal(plt_5[["labels"]][["y"]], selected_metric)
  expect_equal(plt_5[["labels"]][["title"]], 
               paste0(obj_depmap_feat_2[["selected_feat_meta_col"]], ": all NAs"))
  expect_equal(plt_5[["labels"]][["caption"]], unique(dt_response_na$rId))
  
  # scenario: xc50 capped
  selected_metric <- "RV_gDR_log10_xc50"
  dt_response <- dt_response_met[, c("rId", "cId", "CellLineName", selected_metric), with = FALSE]
  dt_response_capped <- 
    dt_response_met_capped[, c("rId", "cId", "CellLineName", selected_metric), with = FALSE]
  
  res_count_6 <- obj_depmap_feat_2[["dt_depmap"]][CCLEName %in% dt_response$CellLineName]
  res_count_6 <- 
    res_count_6[!is.na(get(selected_feat)), .N, by = selected_feat][, lbl := sprintf("%s (%s)", get(selected_feat), N)]
  
  plt_6_raw <- plot_boxplot_num(dt_response = dt_response,
                                dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                                selected_feat = selected_feat)
  plt_6_cap <- plot_boxplot_num(dt_response = dt_response_capped,
                                dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                                selected_feat = selected_feat)
  expect_is(plt_6_raw, "gg")
  expect_is(plt_6_cap, "gg")
  expect_length(plt_6_raw[["layers"]], 3)
  expect_length(plt_6_cap[["layers"]], 3)
  expect_equal(plt_6_raw[["labels"]][["x"]], plt_6_cap[["labels"]][["x"]])
  expect_equal(plt_6_raw[["labels"]][["y"]], plt_6_cap[["labels"]][["y"]])
  expect_true(any(is.infinite(ggplot2::ggplot_build(plt_6_raw)$data[[3]]$y)))
  expect_false(any(is.infinite(ggplot2::ggplot_build(plt_6_cap)$data[[3]]$y)))
  expect_equal(sort(ggplot2::layer_scales(plt_6_raw)$x$get_labels()), sort(res_count_6$lbl))
  expect_equal(sort(ggplot2::layer_scales(plt_6_cap)$x$get_labels()), sort(res_count_6$lbl))
  
  # testing assertions
  expect_error(plot_boxplot_num(dt_response = unlist(dt_response),
                                dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                                selected_feat = selected_feat),
               "Assertion on 'dt_response' failed: Must be a data.table")
  expect_error(plot_boxplot_num(dt_response = dt_response,
                                dt_depmap = obj_depmap_feat_2, 
                                selected_feat = selected_feat),
               "Assertion on 'dt_depmap' failed: Must be a data.table")
  expect_error(plot_boxplot_num(dt_response = dt_response,
                                dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                                selected_feat = 1),
               "Assertion on 'selected_feat' failed: Must be of type 'string'")
  expect_error(plot_boxplot_num(dt_response = dt_response,
                                dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                                selected_feat = "not_existen_feat"),
               "Assertion on 'names\\(dt_depmap\\)' failed: Names must include the elements")
  expect_error(plot_boxplot_num(dt_response = dt_response,
                                dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                                selected_feat = selected_feat,
                                selected_feat_meta_col = 1),
               "Assertion on 'selected_feat_meta_col' failed: Must be of type 'string'")
})

test_that("plot_boxplot_num_panel works as expected", {
  selected_feats <- c("NU_X1QW", "NU_X3OP", "NU_X5BN")
  selected_feats <- c("NU_X1QW (456)", "NU_X3OP (523)", "NU_X5BN")
  selected_feats <- ls_feat_2[3:5]
  selected_metric <- "RV_gDR_x_10"
  dt_response <- dt_response_dose[, c("rId", "cId", "CellLineName", selected_metric), with = FALSE]
  
  plt_1 <- 
    plot_boxplot_num_panel(dt_response = dt_response,
                           dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                           selected_feats = selected_feats) # default
  expect_is(plt_1, "gg")
  expect_length(plt_1[["layers"]], 3)
  expect_equal(plt_1[["labels"]][["x"]], "")
  expect_equal(plt_1[["labels"]][["y"]], selected_metric)
  expect_equal(plt_1[["labels"]][["title"]], NULL)
  expect_equal(plt_1[["labels"]][["caption"]], unique(dt_response$rId))
  
  # selected feat is not present in dt_depmap
  new_selected_feats <- c(selected_feats, "NU_non_avial")
  plt_2 <- 
    plot_boxplot_num_panel(dt_response = dt_response,
                           dt_depmap = obj_depmap_feat_2[["dt_depmap"]],
                           selected_feats = new_selected_feats,
                           selected_feat_meta_col = obj_depmap_feat_2[["selected_feat_meta_col"]])
  expect_is(plt_2, "gg")
  expect_length(plt_2[["layers"]], 3)
  expect_equal(plt_2[["labels"]][["y"]], selected_metric)
  expect_equal(NROW(ggplot2::ggplot_build(plt_2)[["data"]][[1]]), NROW(new_selected_feats))
  expect_equal(# check the uniqueness of points
    unique(data.table::data.table(ggplot2::ggplot_build(plt_2)[["data"]][[3]])[!is.na(x) & !is.na(y), .N, by = "PANEL"]$N), #nolint
    sum(stats::complete.cases(obj_depmap_feat_2[["dt_depmap"]][CCLEName %in% dt_response[["CellLineName"]], ])))
  
  # only NAs in selected_feats 
  plt_3 <- 
    plot_boxplot_num_panel(dt_response = dt_response,
                           dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                           selected_feats = rep(NA, 2),
                           selected_feat_meta_col = obj_depmap_feat_2[["selected_feat_meta_col"]])
  expect_is(plt_3, "gg")
  expect_length(plt_3[["layers"]], 0) # empty plot
  expect_equal(plt_3[["labels"]][["x"]], "")
  expect_equal(plt_3[["labels"]][["y"]], selected_metric)
  expect_equal(plt_3[["labels"]][["title"]], 
               paste0(obj_depmap_feat_2[["selected_feat_meta_col"]], ": all NAs"))
  
  # some NAs in selected_feats 
  selected_feats_with_NAs <- c(NA, selected_feats, NA) 
  plt_4 <- 
    plot_boxplot_num_panel(dt_response = dt_response,
                           dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                           selected_feats = selected_feats_with_NAs,
                           ncol = 2)
  expect_is(plt_4, "gg")
  expect_length(plt_4[["layers"]], 3)
  expect_equal(plt_4[["labels"]][["y"]], selected_metric)
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_4)[["data"]][[1]]$PANEL)),
               NROW(selected_feats_with_NAs))
  expect_equal(plt_4[["facet"]][["params"]][["ncol"]], 2)
  
  # NAs in response
  dt_response_na <- data.table::copy(dt_response)
  dt_response_na[[selected_metric]] <- NA
  plt_5 <-
    plot_boxplot_num_panel(dt_response = dt_response_na,
                           dt_depmap = obj_depmap_feat_2[["dt_depmap"]],
                           selected_feats = selected_feats)
  expect_is(plt_5, "gg")
  expect_length(plt_5[["layers"]], 0) # empty plot
  expect_equal(plt_5[["labels"]][["x"]], "")
  expect_equal(plt_5[["labels"]][["y"]], selected_metric)
  expect_true(grepl("all NAs", plt_5[["labels"]][["title"]]))
  expect_equal(NROW(ggplot2::ggplot_build(plt_5)[["data"]][[1]]), 0)
  
  # NAs in depmap
  dt_depmap_na <- data.table::copy(obj_depmap_feat_2[["dt_depmap"]])
  dt_depmap_na[[selected_feats[2]]] <- NA
  plt_6 <- 
    plot_boxplot_num_panel(dt_response = dt_response,
                           dt_depmap = dt_depmap_na, 
                           selected_feats = selected_feats,
                           ncol = 1)
  expect_is(plt_6, "gg")
  expect_equal(plt_6[["labels"]][["x"]], "")
  expect_equal(plt_6[["labels"]][["y"]], selected_metric)
  expect_equal(NROW(ggplot2::ggplot_build(plt_6)[["data"]][[1]]), NROW(selected_feats)) 
  expect_equal(plt_6[["facet"]][["params"]][["ncol"]], 1)
  
  # capped
  selected_feats_cap <- ls_feat_2[1:3]
  selected_metric <- "RV_gDR_log10_xc50"
  dt_response <- dt_response_met[, c("rId", "cId", "CellLineName", selected_metric), with = FALSE]
  dt_response_capped <- 
    dt_response_met_capped[, c("rId", "cId", "CellLineName", selected_metric), with = FALSE]
  
  plt_7_raw <- plot_boxplot_num_panel(dt_response = dt_response,
                                      dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                                      selected_feats = selected_feats_cap)
  plt_7_cap <- plot_boxplot_num_panel(dt_response = dt_response_capped,
                                      dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                                      selected_feats = selected_feats_cap)
  expect_is(plt_7_raw, "gg")
  expect_is(plt_7_cap, "gg")
  expect_length(plt_7_raw[["layers"]], 3)
  expect_length(plt_7_cap[["layers"]], 3)
  expect_equal(plt_7_raw[["labels"]][["x"]], plt_7_cap[["labels"]][["x"]])
  expect_equal(plt_7_raw[["labels"]][["y"]], plt_7_cap[["labels"]][["y"]])
  expect_true(any(is.infinite(ggplot2::ggplot_build(plt_7_raw)$data[[3]]$y)))
  expect_false(any(is.infinite(ggplot2::ggplot_build(plt_7_cap)$data[[3]]$y)))
  
  # all selected_feats are not avialave in data
  plt_8 <- 
    plot_boxplot_num_panel(dt_response = dt_response,
                           dt_depmap = dt_depmap_na, 
                           selected_feats = LETTERS[1:5])
  expect_is(plt_8, "gg")
  expect_length(plt_8[["layers"]], 0) # empty plot
  expect_true(grepl("all NAs", plt_8[["labels"]][["title"]]))
  
  # there are more category - one feat is not avialable
  dt_depmap_multi <- data.table::copy(obj_depmap_feat_2[["dt_depmap"]])
  dt_depmap_multi[[selected_feats[2]]][1:10] <- 2
  
  plt_9 <- 
    plot_boxplot_num_panel(dt_response = dt_response,
                           dt_depmap = dt_depmap_multi, 
                           selected_feats = new_selected_feats)
  expect_is(plt_9, "gg")
  expect_length(plt_9[["layers"]], 3)
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_9)[["data"]][[2]]$x)),
               sum(!is.na(unique(as.vector(as.matrix(
                 dt_depmap_multi[, .SD, .SDcols = intersect(new_selected_feats, names(dt_depmap_multi))])))))
  )
  
  # testing assertions
  expect_error(plot_boxplot_num_panel(dt_response = unlist(dt_response),
                                      dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                                      selected_feats = selected_feats),
               "Assertion on 'dt_response' failed: Must be a data.table")
  expect_error(plot_boxplot_num_panel(dt_response = dt_response,
                                      dt_depmap = obj_depmap_feat_2, 
                                      selected_feats = selected_feats),
               "Assertion on 'dt_depmap' failed: Must be a data.table")
  expect_error(plot_boxplot_num_panel(dt_response = dt_response,
                                      dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                                      selected_feats = 1),
               "Assertion on 'selected_feats' failed: Must be of type 'character'")
  expect_error(plot_boxplot_num_panel(dt_response = dt_response,
                                      dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                                      selected_feats = selected_feats,
                                      selected_feat_meta_col = 1),
               "Assertion on 'selected_feat_meta_col' failed: Must be of type 'string'")
  expect_error(plot_boxplot_num_panel(dt_response = dt_response,
                                      dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                                      selected_feats = selected_feats,
                                      ncol = "1"),
               "Assertion on 'ncol' failed: Must be of type 'single integerish value'")
  expect_error(plot_boxplot_num_panel(dt_response = dt_response,
                                      dt_depmap = obj_depmap_feat_2[["dt_depmap"]], 
                                      selected_feats = selected_feats,
                                      ncol = 2.5),
               "Assertion on 'ncol' failed: Must be of type 'single integerish value'")
  
})

test_that("plot_boxplot_meta works as expected", {
  id_col <- c("ModelID", "CCLEName")
  
  selected_metric_1 <- "RV_gDR_log10_xc50"
  dt_response_1 <- dt_response_met[, c("rId", "cId", "CellLineName", selected_metric_1), with = FALSE]
  dt_response_capped_1 <-
    dt_response_met_capped[, c("rId", "cId", "CellLineName", selected_metric_1), with = FALSE]
  
  selected_metric_2 <- "RV_gDR_x_10"
  dt_response_2 <- dt_response_dose[, c("rId", "cId", "CellLineName", selected_metric_2), with = FALSE]
  
  # scenairo: default
  plt_1 <- 
    plot_boxplot_meta(dt_response = dt_response_1,
                      dt_depmap = obj_depmap_meta[["dt_depmap"]], 
                      selected_feat_meta_col = obj_depmap_meta[["selected_feat_meta_col"]]) # default
  expect_is(plt_1, "gg")
  expect_equal(plt_1[["labels"]][["y"]], selected_metric_1)
  expect_equal(plt_1[["labels"]][["title"]], obj_depmap_meta[["selected_feat_meta_col"]])
  expect_length(plt_1[["layers"]], 4)
  common_cellline_1 <- merge(dt_response_1[!is.infinite(get(selected_metric_1))], 
                             obj_depmap_meta[["dt_depmap"]], 
                             by.x = "CellLineName", by.y = "CCLEName")[["CellLineName"]]
  res_1 <- obj_depmap_meta[["dt_depmap"]][CCLEName %in% common_cellline_1, .SD, .SDcols = -id_col]
  expect_equal(sort(ggplot2::layer_scales(plt_1)$x$range$range),
               sort(names(res_1)[colSums(res_1) > 0]))
  expect_equal(
    sort(ggplot2::ggplot_build(plt_1)$data[[3]]$y),
    sort(dt_response_1[!is.infinite(get(selected_metric_1))][[selected_metric_1]]))
  expect_equal(ggplot2::ggplot_build(plt_1)$data[[4]]$label, 
               c(colSums(res_1)[colSums(res_1) > 1], use.names = FALSE))
  
  # scenario: Inf value
  plt_1_inf <- 
    plot_boxplot_meta(dt_response = dt_response_1,
                      dt_depmap = obj_depmap_meta[["dt_depmap"]], 
                      selected_feat_meta_col = obj_depmap_meta[["selected_feat_meta_col"]], 
                      with_inf = TRUE)
  expect_is(plt_1_inf, "gg")
  common_cellline_1_inf <- merge(dt_response_1, 
                                 obj_depmap_meta[["dt_depmap"]], 
                                 by.x = "CellLineName", by.y = "CCLEName")[["CellLineName"]]
  res_1_inf <- obj_depmap_meta[["dt_depmap"]][CCLEName %in% common_cellline_1_inf, .SD, .SDcols = -id_col]
  expect_equal(sort(ggplot2::layer_scales(plt_1_inf)$x$range$range),
               sort(names(res_1)[colSums(res_1) > 0]))
  expect_equal(sort(ggplot2::ggplot_build(plt_1_inf)$data[[3]]$y),
               sort(dt_response_1[[selected_metric_1]]))
  expect_equal(ggplot2::ggplot_build(plt_1_inf)$data[[4]]$label, 
               c(colSums(res_1_inf)[colSums(res_1_inf) > 1], use.names = FALSE))
  
  # scenario: capped values
  plt_1_cap <- 
    plot_boxplot_meta(dt_response = dt_response_capped_1,
                      dt_depmap = obj_depmap_meta[["dt_depmap"]], 
                      selected_feat_meta_col = obj_depmap_meta[["selected_feat_meta_col"]]) # default
  expect_is(plt_1_cap, "gg")
  expect_equal(plt_1_cap[["labels"]][["y"]], selected_metric_1)
  expect_equal(plt_1_cap[["labels"]][["title"]], obj_depmap_meta[["selected_feat_meta_col"]])
  expect_length(plt_1_cap[["layers"]], 4)
  expect_equal(sort(ggplot2::layer_scales(plt_1_cap)$x$range$range),
               sort(names(res_1)[colSums(res_1) > 0]))
  expect_equal(sort(ggplot2::ggplot_build(plt_1_cap)$data[[3]]$y),
               sort(dt_response_capped_1[[selected_metric_1]]))
  expect_equal(ggplot2::ggplot_build(plt_1_cap)$data[[4]]$label, 
               c(colSums(res_1_inf)[colSums(res_1_inf) > 1], use.names = FALSE))
  
  # scenario: x-label with max 8 character
  obj_depmap_meta_lng <- 
    prep_dt_depmap_meta(meta_data_path = system.file("testdata/Model.csv", package = "gDRplots"),
                        metadata_col = "TreatmentStatus")
  plt_2 <- 
    plot_boxplot_meta(dt_response = dt_response_2,
                      dt_depmap = obj_depmap_meta_lng[["dt_depmap"]], 
                      selected_feat_meta_col = obj_depmap_meta_lng[["selected_feat_meta_col"]],
                      max_x_lbl_length = 8)
  expect_is(plt_2, "gg")
  expect_length(plt_2[["layers"]], 5)
  common_cellline_3 <- merge(dt_response_2, 
                             obj_depmap_meta_lng[["dt_depmap"]], 
                             by.x = "CellLineName", by.y = "CCLEName")[["CellLineName"]]
  res_2 <- obj_depmap_meta_lng[["dt_depmap"]][CCLEName %in% common_cellline_1_inf, .SD, .SDcols = -id_col]
  expect_equal(sort(ggplot2::layer_scales(plt_2)$x$range$range),
               sort(names(res_2)[colSums(res_2) > 0]))
  
  ls_x_lbl <- unique(ggplot2::layer_scales(plt_2)$x$labels)
  short_lbl <- paste0(
    substr(names(res_2)[colSums(res_2) > 0][nchar(names(res_2)[colSums(res_2) > 0]) > 8], 1, 8 - 3),
    "...")
  expect_true(
    all(names(res_2)[colSums(res_2) > 0][nchar(names(res_2)[colSums(res_2) > 0]) <= 8] %in% ls_x_lbl))
  expect_true(all(short_lbl %in% ls_x_lbl))
  expect_equal(ggplot2::ggplot_build(plt_2)[["data"]][[5]]$yintercept, 1)
  
  # scenario: plot without one-item-group
  plt_3_with_1 <- 
    plot_boxplot_meta(dt_response = dt_response_2,
                      dt_depmap = obj_depmap_meta_lng[["dt_depmap"]], 
                      selected_feat_meta_col = obj_depmap_meta_lng[["selected_feat_meta_col"]])
  
  plt_3_without_1 <- 
    plot_boxplot_meta(dt_response = dt_response_2,
                      dt_depmap = obj_depmap_meta_lng[["dt_depmap"]], 
                      selected_feat_meta_col = obj_depmap_meta_lng[["selected_feat_meta_col"]],
                      with_1_item_grp = FALSE)
  expect_is(plt_3_with_1, "gg")
  expect_is(plt_3_without_1, "gg")
  expect_length(plt_3_with_1[["layers"]], 5)
  expect_length(plt_3_without_1[["layers"]], 5)
  common_cellline_3 <- merge(dt_response_2, 
                             obj_depmap_meta_lng[["dt_depmap"]], 
                             by.x = "CellLineName", by.y = "CCLEName")[["CellLineName"]]
  res_3_with_1 <- 
    obj_depmap_meta_lng[["dt_depmap"]][CCLEName %in% common_cellline_3, .SD, .SDcols = -id_col]
  expect_equal(sort(ggplot2::layer_scales(plt_3_with_1)$x$range$range),
               sort(names(res_3_with_1)[colSums(res_3_with_1) > 0]))
  expect_equal(ggplot2::ggplot_build(plt_3_with_1)$data[[4]]$label, 
               c(colSums(res_3_with_1)[colSums(res_3_with_1) > 0], use.names = FALSE))
  expect_true(1 %in% c(colSums(res_3_with_1)[colSums(res_3_with_1) > 0], use.names = FALSE))
  res_3_without_1 <- 
    obj_depmap_meta_lng[["dt_depmap"]][CCLEName %in% common_cellline_3, .SD, .SDcols = -id_col]
  expect_equal(sort(ggplot2::layer_scales(plt_3_without_1)$x$range$range),
               sort(names(res_3_without_1)[colSums(res_3_without_1) > 1]))
  expect_equal(ggplot2::ggplot_build(plt_3_without_1)$data[[4]]$label, 
               c(colSums(res_3_without_1)[colSums(res_3_without_1) > 1], use.names = FALSE))
  expect_false(1 %in% c(colSums(res_3_without_1)[colSums(res_3_without_1) > 1], use.names = FALSE))
  
  # scenario: plot with combo metric
  selected_metric_3 <- "RV_gDR_bliss_score"
  dt_response_3 <- dt_response_score[, c("rId", "cId", "CellLineName", selected_metric_3), with = FALSE]
  
  plt_4 <- 
    plot_boxplot_meta(dt_response = dt_response_3,
                      dt_depmap = obj_depmap_meta_lng[["dt_depmap"]], 
                      selected_feat_meta_col = obj_depmap_meta_lng[["selected_feat_meta_col"]]) # default
  expect_is(plt_4, "gg")
  expect_length(plt_4[["layers"]], 4) # max(dt_response_3[["RV_gDR_bliss_score"]]) < 0.5 # nolint
  expect_equal(plt_4[["labels"]][["y"]], selected_metric_3)
  expect_equal(plt_4[["labels"]][["title"]], obj_depmap_meta_lng[["selected_feat_meta_col"]])
  expect_equal(plt_4[["labels"]][["caption"]], unique(dt_response_3$rId))
  
  # scenario: data with numeric levels
  obj_depmap_meta_num <- 
    prep_dt_depmap_meta(meta_data_path = system.file("testdata/Model.csv", package = "gDRplots"),
                        metadata_col = "Age")
  
  plt_5 <- 
    plot_boxplot_meta(dt_response = dt_response_2,
                      dt_depmap = obj_depmap_meta_num[["dt_depmap"]],
                      selected_feat_meta_col = obj_depmap_meta_num[["selected_feat_meta_col"]])
  expect_is(plt_5, "gg")
  expect_equal(plt_5[["labels"]][["y"]], selected_metric_2)
  common_cellline_5 <- merge(dt_response_2, 
                             obj_depmap_meta_num[["dt_depmap"]], 
                             by.x = "CellLineName", by.y = "CCLEName")[["CellLineName"]]
  res_5 <- obj_depmap_meta_num[["dt_depmap"]][CCLEName %in% common_cellline_5, .SD, .SDcols = -id_col]
  expect_equal(ggplot2::layer_scales(plt_5)$x$range$range,
               sort(names(res_5)[colSums(res_5) > 0]))
  
  # scenario: data with logical levels
  obj_depmap_meta_log <- 
    prep_dt_depmap_meta(meta_data_path = system.file("testdata/Model.csv", package = "gDRplots"),
                        metadata_col = "SourceDetail")
  
  plt_6 <- 
    plot_boxplot_meta(dt_response = dt_response_3,
                      dt_depmap = obj_depmap_meta_log[["dt_depmap"]],
                      selected_feat_meta_col = obj_depmap_meta_log[["selected_feat_meta_col"]])
  expect_is(plt_6, "gg")
  expect_equal(plt_6[["labels"]][["y"]], selected_metric_3)
  common_cellline_6 <- merge(dt_response_3, 
                             obj_depmap_meta_log[["dt_depmap"]], 
                             by.x = "CellLineName", by.y = "CCLEName")[["CellLineName"]]
  res_6 <- obj_depmap_meta_log[["dt_depmap"]][CCLEName %in% common_cellline_6, .SD, .SDcols = -id_col]
  expect_equal(ggplot2::layer_scales(plt_6)$x$range$range,
               sort(names(res_6)[colSums(res_6) > 0]))
  
  # scenario: lack of the intersection (empty plot)
  dt_depmap_meta_empty <- data.table::copy(obj_depmap_meta[["dt_depmap"]])
  meta_name <- setdiff(names(dt_depmap_meta_empty), c("CCLEName", "ModelID"))
  dt_depmap_meta_empty <- dt_depmap_meta_empty[, (meta_name) := 0]
  
  plt_7 <- 
    plot_boxplot_meta(dt_response = dt_response_1,
                      dt_depmap = dt_depmap_meta_empty, 
                      obj_depmap_meta_num[["selected_feat_meta_col"]])
  expect_is(plt_7, "gg")
  expect_equal(plt_7[["labels"]][["y"]], selected_metric_1)
  expect_true(grepl(obj_depmap_meta_num[["selected_feat_meta_col"]], 
                    plt_7[["labels"]][["title"]]))
  expect_true(grepl("all NAs", plt_7[["labels"]][["title"]]))
  expect_equal(NROW(ggplot2::ggplot_build(plt_7)$data[[1]]), 0)
  
  # scenario: lack of one-to-one relationship
  dt_depmap_meta_multi <- data.table::copy(obj_depmap_meta[["dt_depmap"]])
  meta_name <- setdiff(names(dt_depmap_meta_multi), c("CCLEName", "ModelID"))
  dt_depmap_meta_multi <- dt_depmap_meta_multi[, (meta_name[seq_len(2)]) := 1]
  
  expect_warning({
    plt_8 <- 
      plot_boxplot_meta(dt_response = dt_response_1,
                        dt_depmap = dt_depmap_meta_multi, 
                        selected_feat_meta_col = obj_depmap_meta[["selected_feat_meta_col"]])
    expect_is(plt_8, "gg")
    expect_equal(plt_8[["labels"]][["y"]], selected_metric_1)
    common_cellline_8 <- merge(dt_response_1[!is.infinite(get(selected_metric_1)), ], 
                               dt_depmap_meta_multi, 
                               by.x = "CellLineName", by.y = "CCLEName")[["CellLineName"]]
    res_8 <- dt_depmap_meta_multi[CCLEName %in% common_cellline_8, .SD, .SDcols = -id_col]
    expect_true(NROW(ggplot2::ggplot_build(plt_8)$data[[3]]) == sum(colSums(res_8)))
    expect_true(NROW(ggplot2::ggplot_build(plt_8)$data[[3]]) > sum(colSums(res_1)))
    expect_equal(ggplot2::layer_scales(plt_8)$x$range$range,
                 sort(names(res_8)[colSums(res_8) > 0]))
  })
  expect_warning(
    plot_boxplot_meta(dt_response = dt_response_1,
                      dt_depmap = dt_depmap_meta_multi, 
                      selected_feat_meta_col = obj_depmap_meta[["selected_feat_meta_col"]]),
    "The data does not appear to be categorical"
  )
  
  # testing assertions
  expect_error(plot_boxplot_meta(dt_response = unlist(dt_response_1),
                                 dt_depmap = obj_depmap_meta[["dt_depmap"]],
                                 selected_feat_meta_col = selected_meta),
               "Assertion on 'dt_response' failed: Must be a data.table")
  expect_error(plot_boxplot_meta(dt_response = dt_response_1,
                                 dt_depmap = unlist(obj_depmap_meta), 
                                 selected_feat_meta_col = selected_meta),
               "Assertion on 'dt_depmap' failed: Must be a data.table")
  expect_error(plot_boxplot_meta(dt_response = dt_response_1,
                                 dt_depmap = obj_depmap_meta[["dt_depmap"]], 
                                 selected_feat_meta_col = 1),
               "Assertion on 'selected_feat_meta_col' failed: Must be of type 'string'")
  expect_error(plot_boxplot_meta(dt_response = dt_response_1,
                                 dt_depmap = obj_depmap_meta[["dt_depmap"]], 
                                 selected_feat_meta_col = obj_depmap_meta[["selected_feat_meta_col"]],
                                 with_1_item_grp = "str"),
               "Assertion on 'with_1_item_grp' failed: Must be of type 'logical flag'")
  expect_error(plot_boxplot_meta(dt_response = dt_response_1,
                                 dt_depmap = obj_depmap_meta[["dt_depmap"]], 
                                 selected_feat_meta_col = obj_depmap_meta[["selected_feat_meta_col"]],
                                 max_x_lbl_length = "ten"),
               "Assertion on 'max_x_lbl_length' failed: Must be of type 'number'")
  expect_error(plot_boxplot_meta(dt_response = dt_response_1,
                                 dt_depmap = obj_depmap_meta[["dt_depmap"]], 
                                 selected_feat_meta_col = obj_depmap_meta[["selected_feat_meta_col"]],
                                 max_x_lbl_length = 1:5),
               "Assertion on 'max_x_lbl_length' failed: Must have length 1")
  expect_error(plot_boxplot_meta(dt_response = dt_response_1,
                                 dt_depmap = obj_depmap_meta[["dt_depmap"]], 
                                 selected_feat_meta_col = obj_depmap_meta[["selected_feat_meta_col"]],
                                 with_inf = "yes"),
               "Assertion on 'with_inf' failed: Must be of type 'logical flag'")
}) 


test_that("plot_volcano_assoc_panel works as expected", {
  plt_1 <- plot_volcano_assoc_panel(dt_response = dt_response_met,
                                    dt_depmap = obj_depmap_feat[["dt_depmap"]],
                                    selected_metric = "RV_gDR_x_max",  
                                    selected_feat_meta_col = obj_depmap_feat[["selected_feat_meta_col"]])
  expect_is(plt_1, "gg")
  expect_true(any(grepl("PANEL", names(ggplot2::ggplot_build(plt_1)[["data"]][[1]]))))
  
  plt_2 <- plot_volcano_assoc_panel(dt_response = dt_response_score,
                                    dt_depmap = obj_depmap_feat_2[["dt_depmap"]],
                                    selected_metric = "RV_gDR_bliss_score",  
                                    selected_feat_meta_col = obj_depmap_feat_2[["selected_feat_meta_col"]])
  expect_is(plt_2, "gg")
  expect_true(any(grepl("PANEL", names(ggplot2::ggplot_build(plt_2)[["data"]][[1]]))))
  
  plt_3 <- plot_volcano_assoc_panel(dt_response = dt_response_diff,
                                    dt_depmap = obj_depmap_meta[["dt_depmap"]],
                                    selected_metric = "RV_gDR_x_max_cotrt_diff_0.1_col_fittings",  
                                    selected_feat_meta_col = obj_depmap_meta[["selected_feat_meta_col"]])
  expect_is(plt_3, "gg")
  expect_true(any(grepl("PANEL", names(ggplot2::ggplot_build(plt_3)[["data"]][[1]]))))
  
  expect_error(plot_volcano_assoc_panel(dt_response = unlist(dt_response_met),
                                        dt_depmap = obj_depmap_feat[["dt_depmap"]],
                                        selected_metric = "RV_gDR_x_max",  
                                        selected_feat_meta_col = obj_depmap_feat[["selected_feat_meta_col"]]),
               "Assertion on 'dt_response' failed: Must be a data.table")
  expect_error(plot_volcano_assoc_panel(dt_response = dt_response_met,
                                        dt_depmap = obj_depmap_feat,
                                        selected_metric = "RV_gDR_x_max",  
                                        selected_feat_meta_col = obj_depmap_feat[["selected_feat_meta_col"]]),
               "Assertion on 'dt_depmap' failed: Must be a data.table")
  expect_error(plot_volcano_assoc_panel(dt_response = dt_response_met,
                                        dt_depmap = obj_depmap_feat[["dt_depmap"]],
                                        selected_metric = 1,  
                                        selected_feat_meta_col = obj_depmap_feat[["selected_feat_meta_col"]]),
               "Assertion on 'selected_metric' failed: Must be of type 'string'")
  expect_error(plot_volcano_assoc_panel(dt_response = dt_response_met,
                                        dt_depmap = obj_depmap_feat[["dt_depmap"]],
                                        selected_metric = "RV_gDR_x_max",  
                                        selected_feat_meta_col = NA),
               "Assertion on 'selected_feat_meta_col' failed: May not be NA.")
  expect_error(plot_volcano_assoc_panel(dt_response = dt_response_met,
                                        dt_depmap = obj_depmap_feat[["dt_depmap"]],
                                        selected_metric = "not_known_met",  
                                        selected_feat_meta_col = obj_depmap_feat[["selected_feat_meta_col"]]),
               "Assertion on 'names\\(dt_response\\)' failed: Names must include the elements")
  expect_error(plot_volcano_assoc_panel(dt_response = dt_response_met,
                                        dt_depmap = obj_depmap_feat[["dt_depmap"]][, 4:6],
                                        selected_metric = "RV_gDR_x_max",  
                                        selected_feat_meta_col = obj_depmap_feat[["selected_feat_meta_col"]]),
               "Assertion on 'names\\(dt_depmap\\)' failed: Names must include the elements")
})

test_that(".get_data_type works as expected", {
  tab_cat <- data.table::data.table(
    "A" = c(0, 0, 0, 1),
    "B" = c(0, 1, 1, 0),
    "C" = c(1, 0, 0, 0)
  )
  
  tab_cat_2 <- data.table::data.table(
    "A" = c(0, 0, 0, 1),
    "B" = c(0, 1, -1, 0),
    "C" = c(-1, 0, 0, 0)
  )
  
  tab_cat_3 <- data.table::data.table(
    "A" = c(0, 0, 0, 1),
    "B" = c(0, 1, 2, 0),
    "C" = c(2, 0, 0, 0)
  )
  
  tab_cat_na <- data.table::copy(tab_cat)
  tab_cat_na[1:2, c(1, 2)] <- NA
  
  tab_cat_id <- cbind(
    data.table::data.table(id = sprintf("ID_%s", seq_len(NROW(tab_cat)))),
    tab_cat
  )
  
  tab_num <- data.table::data.table(
    "A" = 1:5,
    "B" = 11:15,
    "C" = 101:105
  )
  
  tab_num_na <- data.table::copy(tab_num)
  tab_num_na[1:2, c(1, 2)] <- NA
  
  tab_num_id <- cbind(
    data.table::data.table(
      id = sprintf("ID_%s", seq_len(NROW(tab_num))),
      grp = LETTERS[seq_len(NROW(tab_num))]
    ),
    tab_num
  )
  
  tab_not_cat <- data.table::copy(tab_cat) * 2.12
  
  tab_num_as_cat <- data.table::copy(tab_cat)
  tab_num_as_cat$C <- 1
  
  tab_num_as_cat_2 <- data.table::copy(tab_cat_2)
  tab_num_as_cat_2[1, ][["B"]] <- -1
  
  tab_mix <- data.table::data.table(
    "A" = LETTERS[1:5],
    "B" = 11:15,
    "C" = 101:105
  )
  
  tab_fact <- data.table::data.table(
    "A" = LETTERS[1:5],
    "B" = factor(c(0, 1, 1, 0, 0)),
    "C" = c(1, 0, 0, 1, 0)
  )
  
  expect_equal(.get_data_type(dt_ = tab_cat), "categorical")
  expect_equal(.get_data_type(dt_ = tab_cat_2), "categorical")
  expect_equal(.get_data_type(dt_ = tab_cat_3), "categorical")
  expect_equal(.get_data_type(dt_ = tab_cat_id, desc_col = "id"), "categorical")
  expect_equal(.get_data_type(dt_ = tab_cat_na), "categorical")
  expect_equal(.get_data_type(dt_ = tab_not_cat), "numeric")
  expect_equal(.get_data_type(dt_ = tab_num_as_cat), "num_as_cat")
  expect_equal(.get_data_type(dt_ = tab_num_as_cat_2), "num_as_cat")
  expect_equal(.get_data_type(dt_ = tab_num), "numeric")
  expect_equal(.get_data_type(dt_ = tab_num_na), "numeric")
  expect_equal(.get_data_type(dt_ = tab_num_id, desc_col = c("id", "grp")), "numeric")
  expect_equal(.get_data_type(dt_ = tab_mix), "unknown")
  expect_equal(.get_data_type(dt_ = tab_mix, desc_col = "A"), "numeric")
  expect_equal(.get_data_type(dt_ = tab_fact), "unknown")
  expect_equal(.get_data_type(dt_ = tab_fact, desc_col = "A"), "unknown")
  
  expect_error(.get_data_type(dt_ = NULL),
               "Assertion on 'dt_' failed: Must be a data.table")
  expect_error(.get_data_type(dt_ = tab_cat, desc_col = 1),
               "Assertion on 'desc_col' failed: Must be of type 'character'")
  expect_error(.get_data_type(dt_ = tab_cat, desc_col = "str"),
               "Assertion on 'desc_col' failed: Must be a subset")
})

test_that(".get_n_top_asssoc works as expected", {
  tab_assoc <- data.table::data.table(
    feature = sprintf("ID_%02d", 1:12),
    rho = c(0.216, 0.082, 0.079, 0.067, 0.059, 0.024, 
            0.008, 0.002, -0.097, -0.166, -0.172, -0.245),
    q_value = seq(0.68, 0.02, length.out = 12)
  )
  
  expect_equal(.get_n_top_asssoc(tab_assoc), 
               tab_assoc[order(q_value)][["feature"]][1:4]) # default
  expect_equal(.get_n_top_asssoc(tab_assoc, 14), 
               tab_assoc[order(q_value)][["feature"]])
  expect_equal(.get_n_top_asssoc(tab_assoc[1:3, ]),
               tab_assoc[1:3, ][order(q_value)][["feature"]])
  
  tab_assoc_rho <- data.table::copy(tab_assoc)
  tab_assoc_rho$q_value <- 1
  expect_equal(.get_n_top_asssoc(tab_assoc_rho), 
               tab_assoc_rho[order(-abs(get("rho")))][["feature"]][1:4])
  
  tab_assoc_alpha <- data.table::copy(tab_assoc)
  tab_assoc_alpha$q_value <- 1
  tab_assoc_alpha$rho <- 0.2
  expect_equal(.get_n_top_asssoc(tab_assoc_alpha), tab_assoc_alpha[["feature"]][1:4])
  
  expect_error(.get_n_top_asssoc(dt_ = NULL),
               "Assertion on 'dt_assoc' failed: Must be a data.table")
  expect_error(.get_n_top_asssoc(dt_ = tab_assoc[, 1:2]),
               "failed: Names must include the elements")
  expect_error(.get_n_top_asssoc(dt_ = tab_assoc, n_top = "1"),
               "Assertion on 'n_top' failed: Must be of type 'number'")
  expect_error(.get_n_top_asssoc(dt_ = tab_assoc, n_top = 2:5),
               "Assertion on 'n_top' failed: Must have length 1")
})
