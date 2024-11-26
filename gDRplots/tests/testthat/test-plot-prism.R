context("Test plot-prism")

# data ----
mae <- gDRutils::get_synthetic_data("combo_matrix")
se_sa <- mae[[gDRutils::get_supported_experiments("sa")]]
dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se_sa,
                                               assay_name = "Metrics")
dt_average <- gDRutils::convert_se_assay_to_dt(se = se_sa,
                                               assay_name = "Averaged")

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
dt_response_dose <- 
  prep_dt_response_dose_sa(dt_average, d_name)

dt_response_score <- 
  prep_dt_response_scores(dt_scores, d_name, d_name2, 
                          metric = c("hsa_score", "bliss_score"))
dt_response_diff <-
  prep_dt_response_metric_diff(dt_metrics_combo, d_name, d_name2,
                               metric = c("xc50", "x_mean", "x_max"))


# fake depmap data
cell_lines <- gDRtestData::create_synthetic_cell_lines()[["CellLineName"]]
dt_model <- data.table::data.table(
  ModelID = sprintf("ACH-%06d", seq_along(cell_lines)),
  CCLEName = cell_lines
)
drugs <- gDRtestData::create_synthetic_drugs()[["DrugName"]]

#_feature_sets
dt_depmap_feat <- data.table::data.table(
  CCLEName = cell_lines,
  "XZ_A1QW" = withr::with_seed(42, rnorm(n = NROW(cell_lines), mean = -0.05, sd = 0.11)),
  "XZ_A2GH" = withr::with_seed(42, rnorm(n = NROW(cell_lines), mean = -0.03, sd = 0.13)),
  "XZ_A3OP" = withr::with_seed(42, rnorm(n = NROW(cell_lines), mean = 0.045, sd = 0.10)),
  "XZ_A4RT" = withr::with_seed(42, rnorm(n = NROW(cell_lines), mean = 0.05, sd = 0.10)),
  "XZ_A5BN" = withr::with_seed(42, rnorm(n = NROW(cell_lines), mean = 0.11, sd = 0.13))
)
dt_depmap_feat[CCLEName %in% c("cellline_FD", "cellline_NE"), 
               setdiff(names(dt_depmap_feat), "CCLEName") := NA]
dt_depmap_feat <- merge(dt_model, dt_depmap_feat, by = "CCLEName")
data.table::setkey(dt_depmap_feat, NULL)

obj_depmap_feat <- list(
  dt_depmap = dt_depmap_feat,
  selected_feat_meta_col = "XZ_fatures"
)

#_meta
dt_depmap_meta_lng <- data.table::data.table(
  CCLEName = cell_lines,
  meta_xx = withr::with_seed(42, sample(x = sprintf("meta_%s", c("AA", "BB", "CC")), 
                                        size = NROW(cell_lines), replace = TRUE))
)
dt_depmap_meta_lng[CCLEName %in% c("cellline_OO", "cellline_AA"), ][["meta_xx"]] <- NA
dt_depmap_meta_lng[CCLEName == "cellline_FD", ][["meta_xx"]] <- "meta_DD"
dt_depmap_meta_lng[CCLEName == "cellline_NE", ][["meta_xx"]] <- "longer_than_other_meta_EE"

dt_depmap_meta <- data.table::dcast(data = dt_depmap_meta_lng, 
                                    formula = CCLEName ~ meta_xx, 
                                    fun.aggregate = length)
data.table::setkey(dt_depmap_meta, NULL)
dt_depmap_meta[, "NA" := NULL] # as in kaleidoscope
dt_depmap_meta <- merge(dt_model, dt_depmap_meta, by = "CCLEName")

obj_depmap_meta <- list(
  dt_depmap = dt_depmap_meta,
  selected_feat_meta_col = "meta_xx"
)

#_assoc
ls_feat <- setdiff(names(obj_depmap_feat[["dt_depmap"]]), c("ModelID", "CCLEName"))
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
                     selected_feat_meta_col = "XZ_fatures")

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
                        selected_feat_meta_col = "meta_xx")

# tests ----
test_that("plot_volcano_assoc works as expected", {
  plt_1 <- plot_volcano_assoc(dt_assoc = obj_assoc_sa[["dt_assoc"]],
                              selected_feat_meta_col = obj_assoc_sa[["selected_feat_meta_col"]],
                              selected_metric = obj_assoc_sa[["selected_metric"]])
  
  expect_is(plt_1, "gg")
  expect_length(plt_1[["layers"]], 2)
  expect_equal(plt_1[["labels"]][["x"]], "rho") # predef for x axis
  expect_equal(plt_1[["labels"]][["y"]], "neglog_q_value") # predef for y axis
  expect_equal(plt_1[["labels"]][["title"]], "RV_gDR_xc50__XZ_fatures") # <metric>__<feat>
  
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
  expect_equal(plt_3[["labels"]][["title"]], "hsa_score__meta_xx") # <metric>__<meta>
  
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
  
  # NAs in depmap
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
  
  # check max_N
  dt_assoc_big <- data.table::data.table(
    feature = sprintf("XZ_A%02dTT", 1:50),
    response = rep("RV_gDR_xc50", 50),
    rho = withr::with_seed(42, rnorm(n = 50, mean = 0, sd = 0.035)),
    q_value = withr::with_seed(42, rnorm(n = 50, mean = 0.15, sd = 0.05))
  )
  max_non_stat_sig <- 20
  plt_7 <- plot_volcano_assoc(dt_assoc = dt_assoc_big,
                              selected_metric = obj_assoc_sa[["selected_metric"]],
                              selected_feat_meta_col = obj_assoc_sa[["selected_feat_meta_col"]],
                              condition_info = obj_assoc_sa[["condition_info"]],
                              max_N = max_non_stat_sig)
  expect_is(plt_7, "gg")
  expect_length(plt_7[["layers"]], 2)
  expect_equal(sum(ggplot2::ggplot_build(plt_7)$data[[1]][["label"]] != ""), 10) # default named_p_top
  expect_equal(NROW(ggplot2::ggplot_build(plt_7)$data[[1]]),
               NROW(dt_assoc_big[q_value < 0.05]) + max_non_stat_sig) # default alpha
  
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
  selected_feat <- "XZ_A3OP"
  selected_metric <- "RV_gDR_x_0.01"
  dt_response <- dt_response_dose[, c("rId", "cId", "CellLineName", selected_metric), with = FALSE]
  
  plt_1 <- 
    plot_scatter_with_corr(dt_response = dt_response,
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
  
  plt_2 <- 
    plot_scatter_with_corr(dt_response = dt_response,
                           dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                           selected_feat = selected_feat,
                           selected_feat_meta_col = obj_depmap_feat[["selected_feat_meta_col"]])
  expect_is(plt_2, "gg")
  expect_length(plt_2[["layers"]], 3)
  expect_equal(plt_2[["labels"]][["title"]], obj_depmap_feat[["selected_feat_meta_col"]])
  
  selected_feat_2 <- "XZ_A5BN"
  selected_metric_2 <- "RV_gDR_bliss_score"
  dt_response_2 <- dt_response_score[, c("rId", "cId", "CellLineName", selected_metric_2), with = FALSE]
  
  plt_3 <-
    plot_scatter_with_corr(dt_response = dt_response_2,
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
  plt_4 <- 
    plot_scatter_with_corr(dt_response = dt_response_na,
                           dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                           selected_feat = selected_feat,
                           selected_feat_meta_col = obj_depmap_feat[["selected_feat_meta_col"]])
  expect_is(plt_4, "gg")
  expect_length(plt_4[["layers"]], 0) # empty plot
  expect_equal(plt_4[["labels"]][["x"]], selected_feat)
  expect_equal(plt_4[["labels"]][["y"]], selected_metric)
  expect_equal(plt_4[["labels"]][["title"]], 
               paste(obj_depmap_feat[["selected_feat_meta_col"]], ": all NAs"))
  expect_equal(plt_4[["labels"]][["caption"]], unique(dt_response_na$rId))
  
  # NAs in depmap
  dt_depmap_na <- data.table::copy(obj_depmap_feat[["dt_depmap"]])
  dt_depmap_na[[selected_feat]] <- NA
  plt_5 <- 
    plot_scatter_with_corr(dt_response = dt_response,
                           dt_depmap = dt_depmap_na, 
                           selected_feat = selected_feat,
                           selected_feat_meta_col = obj_depmap_feat[["selected_feat_meta_col"]])
  expect_is(plt_5, "gg")
  expect_length(plt_5[["layers"]], 0) # empty plot
  expect_equal(plt_5[["labels"]][["x"]], selected_feat)
  expect_equal(plt_5[["labels"]][["y"]], selected_metric)
  expect_equal(plt_5[["labels"]][["title"]], 
               paste(obj_depmap_feat[["selected_feat_meta_col"]], ": all NAs"))
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
  selected_feats <- c("XZ_A1QW", "XZ_A3OP", "XZ_A5BN")
  selected_metric <- "RV_gDR_x_10"
  dt_response <- dt_response_dose[, c("rId", "cId", "CellLineName", selected_metric), with = FALSE]
  
  plt_1 <- 
    plot_scatter_with_corr_panel(dt_response = dt_response,
                                 dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                 selected_feats = selected_feats)
  expect_is(plt_1, "gg")
  expect_length(plt_1[["layers"]], 3)
  expect_equal(plt_1[["labels"]][["x"]], "")
  expect_equal(plt_1[["labels"]][["y"]], selected_metric)
  expect_equal(plt_1[["labels"]][["title"]], NULL)
  expect_equal(plt_1[["labels"]][["caption"]], unique(dt_response$rId))
  
  plt_2 <- 
    plot_scatter_with_corr_panel(dt_response = dt_response,
                                 dt_depmap = obj_depmap_feat[["dt_depmap"]], 
                                 selected_feats = selected_feats[1],
                                 selected_feat_meta_col = obj_depmap_feat[["selected_feat_meta_col"]])
  expect_is(plt_2, "gg")
  expect_length(plt_2[["layers"]], 3)
  expect_equal(plt_2[["labels"]][["title"]], "XZ_fatures")
  expect_equal(plt_2[["labels"]][["caption"]], unique(dt_response$rId))
  
  # selected feat is not present in dt_depmap
  new_selected_feats <- c(selected_feats, "XZ_non_avial")
  plt_3 <- 
    plot_scatter_with_corr_panel(dt_response = dt_response,
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
  expect_equal(plt_4[["labels"]][["title"]], "XZ_fatures : all NAs")
  
  # some NAs in selected_feats 
  selected_feats_with_NAs <- c(NA, selected_feats, NA) 
  plt_5 <- 
    plot_scatter_with_corr_panel(dt_response = dt_response,
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
  plt_6 <- 
    plot_scatter_with_corr_panel(dt_response = dt_response_na,
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
  plt_7 <- 
    plot_scatter_with_corr_panel(dt_response = dt_response,
                                 dt_depmap = dt_depmap_na, 
                                 selected_feats = selected_feats)
  expect_is(plt_7, "gg")
  expect_equal(plt_7[["labels"]][["x"]], "")
  expect_equal(plt_7[["labels"]][["y"]], selected_metric)
  expect_equal(plt_7[["labels"]][["title"]], NULL)
  expect_equal(plt_6[["labels"]][["caption"]], unique(dt_response$rId))
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

test_that("plot_boxplot_meta works as expected", {
  selected_meta <- "meta_xx"
  selected_metric <- "RV_gDR_xc50"
  dt_response <- dt_response_met[, c("rId", "cId", "CellLineName", selected_metric), with = FALSE]
  
  grp_stat <- dt_depmap_meta_lng[CCLEName %in% dt_response$CellLineName, .N, by = meta_xx]
  
  plt_1 <- plot_boxplot_meta(dt_response = dt_response,
                             dt_depmap = dt_depmap_meta, 
                             selected_feat_meta_col = selected_meta)
  expect_is(plt_1, "gg")
  expect_equal(plt_1[["labels"]][["y"]], selected_metric)
  expect_equal(plt_1[["labels"]][["title"]], selected_meta)
  expect_length(plt_1[["layers"]], 4)
  expect_length(ggplot2::ggplot_build(plt_1)$data[[2]]$xid,
                NROW(grp_stat[!is.na(meta_xx)]))
  expect_equal(sort(ggplot2::layer_scales(plt_1)$x$range$range),
               sort(grp_stat[!is.na(meta_xx)]$meta_xx))
  
  plt_2 <- plot_boxplot_meta(dt_response = dt_response,
                             dt_depmap = dt_depmap_meta, 
                             selected_feat_meta_col = selected_meta,
                             with_1_item_grp = FALSE)
  expect_is(plt_2, "gg")
  expect_length(plt_2[["layers"]], 4)
  expect_length(ggplot2::ggplot_build(plt_2)$data[[2]]$xid,
                NROW(grp_stat[!is.na(meta_xx) & N > 1]))
  expect_equal(sort(ggplot2::layer_scales(plt_2)$x$range$range), 
               sort(grp_stat[!is.na(meta_xx) & N > 1]$meta_xx))
  
  plt_3 <- plot_boxplot_meta(dt_response = dt_response,
                             dt_depmap = dt_depmap_meta, 
                             selected_feat_meta_col = selected_meta,
                             max_x_lbl_length = 8)
  expect_is(plt_3, "gg")
  expect_length(plt_3[["layers"]], 4)
  expect_length(ggplot2::ggplot_build(plt_3)$data[[2]]$xid,
                NROW(grp_stat[!is.na(meta_xx)]))
  ls_x_lbl <- ggplot2::layer_scales(plt_3)$x$labels
  short_lbl <- paste0(substr(grp_stat[!is.na(meta_xx) & nchar(meta_xx) > 8]$meta_xx, 1, 8 - 3), "...")
  expect_true(all(grp_stat[!is.na(meta_xx) & nchar(meta_xx) < 8]$meta_xx %in% ls_x_lbl))
  expect_true(short_lbl %in% ls_x_lbl)
  expect_equal(sort(ggplot2::layer_scales(plt_3)$x$range$range), 
               sort(grp_stat[!is.na(meta_xx)]$meta_xx))
  
  # combo plot
  selected_metric_2 <- "RV_gDR_bliss_score"
  dt_response_2 <- dt_response_score[, c("rId", "cId", "CellLineName", selected_metric_2), with = FALSE]
  
  plt_4 <- plot_boxplot_meta(dt_response = dt_response_2,
                             dt_depmap = dt_depmap_meta, 
                             selected_feat_meta_col = selected_meta)
  expect_is(plt_4, "gg")
  expect_length(plt_4[["layers"]], 3) # max(dt_response_2[["RV_gDR_bliss_score"]]) < 0.5 # nolint
  expect_equal(plt_4[["labels"]][["y"]], selected_metric_2)
  expect_equal(plt_4[["labels"]][["title"]], selected_meta)
  expect_equal(plt_4[["labels"]][["caption"]], unique(dt_response_2$rId))
  
  # numeric meta
  dt_depmap_meta_numeric <- data.table::copy(dt_depmap_meta)
  meta_name <- setdiff(names(dt_depmap_meta_numeric), c("CCLEName", "ModelID"))
  names(dt_depmap_meta_numeric) <- c("CCLEName", "ModelID", seq_along(meta_name))
  dt_ <- dt_depmap_meta_numeric[CCLEName %in% dt_response[["CellLineName"]], ]
  lbl_ <- vapply(names(dt_[, as.character(seq_along(meta_name)), with = FALSE]), 
                 function(nm) !all(dt_[[nm]] == 0), logical(1))
  
  plt_5 <- plot_boxplot_meta(dt_response = dt_response,
                             dt_depmap = dt_depmap_meta_numeric, 
                             selected_feat_meta_col = selected_meta)
  expect_is(plt_5, "gg")
  expect_equal(plt_5[["labels"]][["y"]], selected_metric)
  expect_equal(ggplot2::layer_scales(plt_5)$x$range$range, names(lbl_)[lbl_])
  
  # testing assertions
  expect_error(plot_boxplot_meta(dt_response = unlist(dt_response),
                                 dt_depmap = dt_depmap_meta,
                                 selected_feat_meta_col = selected_meta),
               "Assertion on 'dt_response' failed: Must be a data.table")
  expect_error(plot_boxplot_meta(dt_response = dt_response,
                                 dt_depmap = unlist(dt_depmap_meta), 
                                 selected_feat_meta_col = selected_meta),
               "Assertion on 'dt_depmap' failed: Must be a data.table")
  expect_error(plot_boxplot_meta(dt_response = dt_response,
                                 dt_depmap = dt_depmap_meta,
                                 selected_feat_meta_col = 1),
               "Assertion on 'selected_feat_meta_col' failed: Must be of type 'string'")
  expect_error(plot_boxplot_meta(dt_response = dt_response,
                                 dt_depmap = dt_depmap_meta, 
                                 selected_feat_meta_col = selected_meta,
                                 with_1_item_grp = "str"),
               "Assertion on 'with_1_item_grp' failed: Must be of type 'logical flag'")
  expect_error(plot_boxplot_meta(dt_response = dt_response,
                                 dt_depmap = dt_depmap_meta, 
                                 selected_feat_meta_col = selected_meta,
                                 max_x_lbl_length = "ten"),
               "Assertion on 'max_x_lbl_length' failed: Must be of type 'number'")
  expect_error(plot_boxplot_meta(dt_response = dt_response,
                                 dt_depmap = dt_depmap_meta, 
                                 selected_feat_meta_col = selected_meta,
                                 max_x_lbl_length = 1:5),
               "Assertion on 'max_x_lbl_length' failed: Must have length 1")
}) 

test_that("plot_volcano_assoc_panel works as expected", {
  # TODO in GDR-2710
})

test_that(".get_data_type works as expected", {
  tab_cat <- data.table::data.table(
    "A" = c(0, 0, 0, 1),
    "B" = c(0, 1, 1, 0),
    "C" = c(1, 0, 0, 0)
  )
  
  tab_cat_na <- data.table::copy(tab_cat)
  tab_cat_na[1:2, c(1, 2)] <- NA
  
  tab_num <- data.table::data.table(
    "A" = 1:5,
    "B" = 11:15,
    "C" = 101:105
  )
  
  tab_num_na <- data.table::copy(tab_num)
  tab_num_na[1:2, c(1, 2)] <- NA
  
  tab_unkn <- data.table::copy(tab_cat) * 2
  
  tab_mix <- data.table::data.table(
    "A" = LETTERS[1:5],
    "B" = 11:15,
    "C" = 101:105
  )
  
  expect_equal(.get_data_type(tab_cat), "categorical")
  expect_equal(.get_data_type(tab_cat_na), "categorical")
  expect_equal(.get_data_type(tab_num), "numeric")
  expect_equal(.get_data_type(tab_num_na), "numeric")
  expect_equal(.get_data_type(tab_unkn), "unknown")
  expect_equal(.get_data_type(tab_mix), "unknown")
  expect_equal(.get_data_type(tab_mix, desc_col = "A"), "numeric")
  
  expect_error(.get_data_type(dt_ = NULL),
               "Assertion on 'dt_' failed: Must be a data.table")
  expect_error(.get_data_type(dt_ = tab_cat, desc_col = 1),
               "Assertion on 'desc_col' failed: Must be of type 'character'")
})
