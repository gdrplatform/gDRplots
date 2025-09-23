context("Test boxplot_metric")

test_that("plot_boxplot_metric_sa works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  
  plt_1 <- plot_boxplot_metric_sa(dt_metrics, 
                                  group_var = "CellLineName") # default
  expect_is(plt_1, "gg")
  expect_length(plt_1[["layers"]], 3)
  expect_true(grepl("GR", plt_1[["labels"]][["y"]]))
  expect_true(grepl("log10", plt_1[["labels"]][["y"]])) # xc50 in log10 scale
  expect_true(grepl("drug", plt_1[["labels"]][["title"]]))
  
  plt_2 <- plot_boxplot_metric_sa(dt_metrics, 
                                  group_var = "DrugName",
                                  normalization_type = "RV",
                                  metric = "x_max",
                                  colors_vec = "darkgreen")
  expect_is(plt_2, "gg")
  expect_true(grepl("E Max", plt_2[["labels"]][["y"]]))
  expect_true(grepl("celllines", plt_2[["labels"]][["title"]]))
  expect_equal(unique(ggplot2::ggplot_build(plt_2)[["data"]][[2]][["fill"]]), "darkgreen")
  
  plt_3 <- plot_boxplot_metric_sa(dt_metrics,
                                  group_var = "CellLineName",
                                  metric = "x_inf",
                                  colors_vec = c("blue", "yellow"))
  expect_is(plt_3, "gg")
  expect_true(grepl("GR", plt_3[["labels"]][["y"]]))
  expect_true(grepl("Inf", plt_3[["labels"]][["y"]]))
  expect_equal(unique(ggplot2::ggplot_build(plt_3)[["data"]][[2]][["fill"]]), "blue")
  
  expect_message(
    plt_4 <- plot_boxplot_metric_sa(dt_metrics, 
                                    group_var = "DrugName",
                                    normalization_type = "RV",
                                    grouped_flag = TRUE,
                                    colored_pts_flag = TRUE),
    "Please, choose only one coloring option"
  )
  expect_is(plt_4, "gg")
  expect_length(plt_4[["layers"]], 4) # grouped_flag
  expect_length(plt_4[["guides"]]$guides, 2) # col_var + point_var
  expect_true(all(c("colour", "fill") %in% names(plt_4[["guides"]]$guides)))

  ls_moa_col <- c("deeppink", "darkcyan", "orange", "darkblue", "gold")
  plt_5 <- plot_boxplot_metric_sa_by_drugs(dt_metrics,
                                           with_inf = TRUE,
                                           grouped_flag = TRUE,
                                           colors_vec = ls_moa_col) 
  expect_is(plt_5, "gg")
  expect_length(plt_5[["layers"]], 4)
  expect_equal(sort(ggplot2::ggplot_build(plt_5)[["data"]][[4]][["y"]]),
               sort(log10(dt_metrics[normalization_type == "GR", ][["xc50"]])))
  expect_equal(sort(ggplot2::get_panel_scales(plt_5)$x$get_labels()),
               sort(unique(dt_metrics[["DrugName"]])))
  expect_equal(unique(ggplot2::ggplot_build(plt_5)[["data"]][[2]][["fill"]]), ls_moa_col)
  
  plt_6 <- plot_boxplot_metric_sa(dt_metrics,
                                  group_var = "CellLineName",
                                  metric = "x_AOC",
                                  colored_pts_flag = TRUE)
  expect_is(plt_6, "gg")
  expect_length(plt_6[["layers"]], 3)
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_6)[["data"]][[3]][["colour"]])),
               NROW(unique(dt_metrics[["DrugName"]])))
  
  mae <- gDRutils::get_synthetic_data("medium")
  se_2 <- mae[[gDRutils::get_supported_experiments("sa")]]
  dt_metrics_2 <- gDRutils::convert_se_assay_to_dt(se_2, "Metrics")
  
  plt_7 <- plot_boxplot_metric_sa(dt_metrics_2,
                                  group_var = "CellLineName",
                                  metric = "x_mean",
                                  normalization_type = "RV",
                                  colored_pts_flag = TRUE)
  expect_is(plt_7, "gg")
  expect_length(plt_7[["layers"]], 3)
  expect_equal(sort(ggplot2::ggplot_build(plt_7)[["data"]][[3]][["y"]]),
               sort(dt_metrics_2[normalization_type == "RV", ][["x_mean"]]))
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_7)[["data"]][[3]][["colour"]])), 
               1) # too many points to be colored
  
  plt_8 <- plot_boxplot_metric_sa(dt_metrics_2,
                                  group_var = "DrugName",
                                  metric = "x_mean",
                                  normalization_type = "RV",
                                  colored_pts_flag = TRUE)
  expect_is(plt_8, "gg")
  expect_length(plt_8[["layers"]], 3)
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_8)[["data"]][[3]][["colour"]])), 
               1) # too many points to be colored
  
  expect_error(plot_boxplot_metric_sa(dt_metrics = unlist(dt_metrics),
                                      group_var = "CellLineName"),
               "Assertion on 'dt_metrics' failed: Must be a data.table")
  expect_error(plot_boxplot_metric_sa(dt_metrics = dt_metrics,
                                      group_var = "unknown"),
               "Assertion on 'group_var' failed: Must be element of set")
  expect_error(plot_boxplot_metric_sa(dt_metrics = dt_metrics,
                                      group_var = "CellLineName",
                                      normalization_type = "XX"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(plot_boxplot_metric_sa(dt_metrics = dt_metrics,
                                      group_var = "CellLineName",
                                      metric = "xxx"),
               "Assertion on 'metric' failed: Must be element of set")
  expect_error(plot_boxplot_metric_sa(dt_metrics = dt_metrics[, -c("CellLineName")],
                                      group_var = "CellLineName"),
               "Assertion on 'names\\(dt_metrics\\)' failed: Names must include the elements")
  expect_error(plot_boxplot_metric_sa(dt_metrics = dt_metrics,
                                      group_var = "CellLineName",
                                      fit_source = 1),
               "Assertion on 'fit_source' failed: Must be of type 'string'")
  expect_error(plot_boxplot_metric_sa(dt_metrics = dt_metrics,
                                      group_var = "CellLineName",
                                      grouped_flag = "yes"),
               "Assertion on 'grouped_flag' failed: Must be of type 'logical flag'")
  expect_error(plot_boxplot_metric_sa(dt_metrics = dt_metrics,
                                      group_var = "CellLineName",
                                      colored_pts_flag = "ADD"),
               "Assertion on 'colored_pts_flag' failed: Must be of type 'logical flag'")
  expect_error(plot_boxplot_metric_sa(dt_metrics = dt_metrics,
                                      group_var = "CellLineName",
                                      colors_vec = 1:3),
               "Assertion on 'colors_vec' failed: Must be of type 'character'")
  expect_error(plot_boxplot_metric_sa(dt_metrics = dt_metrics,
                                      group_var = "CellLineName",
                                      with_inf = "yes"),
               "Assertion on 'with_inf' failed: Must be of type 'logical flag'")
})

test_that("plot_boxplot_metric_sa_by_CLs works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  
  plt_1 <- plot_boxplot_metric_sa_by_CLs(dt_metrics) # default
  
  expect_is(plt_1, "gg")
  expect_true(grepl("GR", plt_1[["labels"]][["y"]]))
  expect_true(grepl("log10", plt_1[["labels"]][["y"]])) # xc50 in log10 scale
  expect_true(grepl("drug", plt_1[["labels"]][["title"]]))
  expect_length(plt_1[["layers"]], 3)
  expect_equal(sort(ggplot2::get_panel_scales(plt_1)$x$get_labels()),
               sort(unique(dt_metrics[["CellLineName"]])))
  
  plt_2 <- plot_boxplot_metric_sa_by_CLs(dt_metrics,
                                         metric = "x_max",
                                         grouped_flag = TRUE)
  
  ls_lbl_x <- unique(dt_metrics[, c("CellLineName", "Tissue"), with = FALSE])
  data.table::setorderv(ls_lbl_x, "Tissue")
  expect_is(plt_2, "gg")
  expect_length(plt_2[["layers"]], 4)
  expect_equal(ggplot2::get_panel_scales(plt_2)$x$get_labels(), ls_lbl_x[["CellLineName"]])
  expect_length(unique(ggplot2::ggplot_build(plt_2)[["data"]][[2]][["fill"]]), 
                NROW(unique(ls_lbl_x[["Tissue"]])))
  expect_true(grepl(NROW(unique(dt_metrics[["DrugName"]])), plt_2[["labels"]][["title"]]))
  
  plt_3 <- plot_boxplot_metric_sa_by_CLs(dt_metrics,
                                         colors_vec = "darkred",
                                         with_inf = TRUE) 
  expect_is(plt_3, "gg")
  expect_length(plt_3[["layers"]], 3)
  expect_equal(sort(ggplot2::ggplot_build(plt_3)[["data"]][[3]][["y"]]),
               sort(log10(dt_metrics[normalization_type == "GR", ][["xc50"]])))
  expect_equal(sort(ggplot2::get_panel_scales(plt_3)$x$get_labels()),
               sort(unique(dt_metrics[["CellLineName"]])))
  expect_equal(unique(ggplot2::ggplot_build(plt_3)[["data"]][[2]][["fill"]]), "darkred")
  
  plt_4 <- plot_boxplot_metric_sa_by_CLs(dt_metrics,
                                         metric = "p_value", 
                                         colored_pts_flag = TRUE) 
  expect_is(plt_4, "gg")
  expect_length(plt_4[["layers"]], 3)
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_4)[["data"]][[3]][["colour"]])),
               NROW(unique(dt_metrics[["DrugName"]])))
})

test_that("plot_boxplot_metric_sa_by_drugs works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  
  plt_1 <- plot_boxplot_metric_sa_by_drugs(dt_metrics) # default
  
  expect_is(plt_1, "gg")
  expect_length(plt_1[["layers"]], 3)
  expect_equal(plt_1[["labels"]][["y"]], get_hm_title("xc50", "GR"))
  expect_true(grepl("cellline", plt_1[["labels"]][["title"]]))
  expect_length(plt_1[["layers"]], 3)
  expect_equal(sort(ggplot2::get_panel_scales(plt_1)$x$get_labels()),
               sort(unique(dt_metrics[["DrugName"]])))
  
  plt_2 <- plot_boxplot_metric_sa_by_drugs(dt_metrics,
                                           normalization_type = "RV",
                                           metric = "x_AOC",
                                           colors_vec = "gold")
  expect_is(plt_2, "gg")
  expect_length(plt_2[["layers"]], 3)
  expect_equal(plt_2[["labels"]][["y"]], get_hm_title("x_AOC", "RV"))
  expect_length(plt_2[["layers"]], 3)
  expect_equal(unique(ggplot2::ggplot_build(plt_2)[["data"]][[2]][["fill"]]), "gold")
  
  plt_3 <- plot_boxplot_metric_sa_by_drugs(dt_metrics,
                                           metric = "x_max",
                                           colors_vec = c("darkred", "yellow"),
                                           colored_pts_flag = TRUE)
  expect_is(plt_3, "gg")
  expect_length(plt_3[["layers"]], 3)
  expect_equal(unique(ggplot2::ggplot_build(plt_3)[["data"]][[2]][["fill"]]), "darkred")
  expect_true(grepl(NROW(unique(dt_metrics[["CellLineName"]])), plt_3[["labels"]][["title"]]))
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_3)[["data"]][[3]][["colour"]])),
               NROW(unique(dt_metrics[["CellLineName"]])))
  
  plt_4 <- plot_boxplot_metric_sa_by_drugs(dt_metrics,
                                           grouped_flag = TRUE,
                                           colors_vec = c("#0000FF", "#00FF00"),
                                           with_inf = TRUE)
  
  ls_lbl_x <- unique(dt_metrics[, c("DrugName", "drug_moa"), with = FALSE])
  data.table::setorderv(ls_lbl_x, "drug_moa")
  expect_is(plt_4, "gg")
  expect_length(plt_4[["layers"]], 4)
  expect_equal(ggplot2::get_panel_scales(plt_4)$x$get_labels(), ls_lbl_x[["DrugName"]])
  expect_length(unique(ggplot2::ggplot_build(plt_4)[["data"]][[2]][["fill"]]), 
                NROW(unique(ls_lbl_x[["drug_moa"]])))
  expect_true(all(c("#0000FF", "#00FF00") %in% unique(ggplot2::ggplot_build(plt_4)[["data"]][[2]][["fill"]])))
  expect_equal(sort(ggplot2::ggplot_build(plt_4)[["data"]][[4]][["y"]]),
               sort(log10(dt_metrics[normalization_type == "GR", ][["xc50"]])))
  expect_equal(sort(ggplot2::get_panel_scales(plt_4)$x$get_labels()),
               sort(unique(dt_metrics[["DrugName"]])))
})

test_that("plot_boxplot_metric_combo works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_scores <- gDRutils::convert_se_assay_to_dt(se = se,
                                                assay_name = "scores")
  ls_comb <- unique(paste(dt_scores[["DrugName"]], "x", dt_scores[["DrugName_2"]]))
  
  plt_1 <- plot_boxplot_metric_combo(dt_scores, 
                                     group_var = "CellLineName") # default
  expect_is(plt_1, "gg")
  expect_length(plt_1[["layers"]], 3)
  expect_equal(plt_1[["labels"]][["y"]], get_hm_title("hsa_score", "GR"))
  expect_true(grepl("drug", plt_1[["labels"]][["title"]]))
  expect_equal(sort(ggplot2::get_panel_scales(plt_1)$x$get_labels()),
               sort(unique(dt_scores[["CellLineName"]])))
  expect_true(grepl(NROW(ls_comb), plt_1[["labels"]][["title"]]))
  
  plt_2 <- plot_boxplot_metric_combo(dt_scores, 
                                     group_var = "DrugName",
                                     normalization_type = "RV",
                                     metric = "bliss_score",
                                     colored_pts_flag = TRUE,
                                     colors_vec = "darkgreen")
  expect_is(plt_2, "gg")
  expect_length(plt_2[["layers"]], 3)
  expect_equal(plt_2[["labels"]][["y"]], get_hm_title("bliss_score", "RV"))
  expect_true(grepl("celllines", plt_2[["labels"]][["title"]]))
  expect_equal(unique(ggplot2::ggplot_build(plt_2)[["data"]][[2]][["fill"]]), "darkgreen")
  expect_true(grepl(NROW(unique(dt_scores[["CellLineName"]])), plt_2[["labels"]][["title"]]))
  expect_equal(sort(ggplot2::get_panel_scales(plt_2)$x$get_labels()),
               sort(ls_comb))
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_2)[["data"]][[3]][["colour"]])),
               NROW(unique(dt_scores[["CellLineName"]])))
  
  plt_3 <- plot_boxplot_metric_combo(dt_scores,
                                     group_var = "CellLineName",
                                     metric = "bliss_score",
                                     colors_vec = c("blue", "yellow"))
  expect_is(plt_3, "gg")
  expect_equal(unique(ggplot2::ggplot_build(plt_3)[["data"]][[2]][["fill"]]), "blue")
  
  expect_message(
    plt_4 <- plot_boxplot_metric_combo(dt_scores, 
                                       group_var = "CellLineName",
                                       normalization_type = "RV",
                                       grouped_flag = TRUE,
                                       colored_pts_flag = TRUE),
    "Please, choose only one coloring option"
  )
  expect_is(plt_4, "gg")
  expect_length(plt_4[["layers"]], 4) # grouped_flag
  expect_length(plt_4[["guides"]]$guides, 2) # col_var + point_var
  expect_true(all(c("colour", "fill") %in% names(plt_4[["guides"]]$guides)))
  
  expect_message(
    plt_5 <- plot_boxplot_metric_combo(dt_scores, 
                                       group_var = "DrugName",
                                       normalization_type = "RV",
                                       grouped_flag = TRUE,
                                       colored_pts_flag = TRUE),
    "Coloring box by group is not available for this scenario."
  )
  expect_is(plt_5, "gg")
  expect_length(plt_5[["layers"]], 3) # grouped_flag ignored
  expect_null(plt_5[["labels"]][["fill"]])
  expect_false("fill" %in% names(plt_5[["guides"]]$guides))
  expect_true("colour" %in% names(plt_5[["guides"]]$guides))
  expect_equal(sort(ggplot2::get_panel_scales(plt_5)$x$get_labels()),
               sort(ls_comb))
  
  dt_scores_2 <- data.table::copy(dt_scores)
  ls_col_x <- c("CellLineName", "clid", "DrugName", "Gnumber")
  dt_scores_2[, (ls_col_x) := lapply(.SD, paste0, "X"), .SDcols = ls_col_x]
  dt_scores_2 <- rbind(dt_scores, dt_scores_2)
  
  plt_6 <- plot_boxplot_metric_combo(dt_scores_2, 
                                     group_var = "DrugName",
                                     normalization_type = "RV",
                                     colored_pts_flag = TRUE)
  expect_is(plt_6, "gg")
  expect_length(plt_6[["layers"]], 3)
  expect_equal(sort(ggplot2::ggplot_build(plt_6)[["data"]][[3]][["y"]]),
               sort(dt_scores_2[normalization_type == "RV", ][["hsa_score"]]))
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_6)[["data"]][[3]][["colour"]])), 
               1) # too many points to be colored
  
  plt_7 <- plot_boxplot_metric_combo(dt_scores_2, 
                                     group_var = "CellLineName",
                                     metric = "bliss_score",
                                     normalization_type = "RV",
                                     colored_pts_flag = TRUE)
  expect_is(plt_7, "gg")
  expect_length(plt_7[["layers"]], 3)
  expect_equal(sort(ggplot2::ggplot_build(plt_7)[["data"]][[3]][["y"]]),
               sort(dt_scores_2[normalization_type == "RV", ][["bliss_score"]]))
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_7)[["data"]][[3]][["colour"]])), 
               1) # too many points to be colored
  
  expect_error(plot_boxplot_metric_combo(dt_scores = unlist(dt_scores),
                                         group_var = "CellLineName"),
               "Assertion on 'dt_scores' failed: Must be a data.table")
  expect_error(plot_boxplot_metric_combo(dt_scores = dt_scores,
                                         group_var = "unknown"),
               "Assertion on 'group_var' failed: Must be element of set")
  expect_error(plot_boxplot_metric_combo(dt_scores = dt_scores,
                                         group_var = "CellLineName",
                                         normalization_type = "XX"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(plot_boxplot_metric_combo(dt_scores = dt_scores,
                                         group_var = "CellLineName",
                                         metric = "xxx"),
               "Assertion on 'metric' failed: Must be element of set")
  expect_error(plot_boxplot_metric_combo(dt_scores = dt_scores[, -c("CellLineName")],
                                         group_var = "CellLineName"),
               "Assertion on 'names\\(dt_scores\\)' failed: Names must include the elements")
  expect_error(plot_boxplot_metric_combo(dt_scores = dt_scores,
                                         group_var = "CellLineName",
                                         fit_source = 1),
               "Assertion on 'fit_source' failed: Must be of type 'string'")
  expect_error(plot_boxplot_metric_combo(dt_scores = dt_scores,
                                         group_var = "CellLineName",
                                         grouped_flag = "yes"),
               "Assertion on 'grouped_flag' failed: Must be of type 'logical flag'")
  expect_error(plot_boxplot_metric_combo(dt_scores = dt_scores,
                                         group_var = "CellLineName",
                                         colored_pts_flag = "TRUE"),
               "Assertion on 'colored_pts_flag' failed: Must be of type 'logical flag'")
  expect_error(plot_boxplot_metric_combo(dt_scores = dt_scores,
                                         group_var = "CellLineName",
                                         colors_vec = 1:3),
               "Assertion on 'colors_vec' failed: Must be of type 'character'")
})

test_that("plot_boxplot_metric_combo_by_CLs works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_scores <- gDRutils::convert_se_assay_to_dt(se = se,
                                                assay_name = "scores")
  ls_comb <- unique(paste(dt_scores[["DrugName"]], "x", dt_scores[["DrugName_2"]]))
  
  plt_1 <- plot_boxplot_metric_combo_by_CLs(dt_scores) # default
  
  expect_is(plt_1, "gg")
  expect_equal(plt_1[["labels"]][["y"]], get_hm_title("hsa_score", "GR"))
  expect_true(grepl("drug", plt_1[["labels"]][["title"]]))
  expect_length(plt_1[["layers"]], 3)
  expect_equal(sort(ggplot2::get_panel_scales(plt_1)$x$get_labels()),
               sort(unique(dt_scores[["CellLineName"]])))
  
  plt_2 <- plot_boxplot_metric_combo_by_CLs(dt_scores,
                                            normalization_type = "RV",
                                            grouped_flag = TRUE,
                                            colors_vec = "#FF0000")
  expect_is(plt_2, "gg")
  expect_equal(plt_2[["labels"]][["y"]], get_hm_title("hsa_score", "RV"))
  expect_length(plt_2[["layers"]], 4)
  expect_equal(unique(ggplot2::ggplot_build(plt_2)[["data"]][[2]][["fill"]]), "#FF0000")
  
  plt_3 <- plot_boxplot_metric_combo_by_CLs(dt_scores,
                                            metric = "hsa_score",
                                            normalization_type = "RV",
                                            grouped_flag = TRUE,
                                            colors_vec = c("deeppink", "darkcyan", "orange", "darkblue"))
  
  ls_lbl_x <- unique(dt_scores[, c("CellLineName", "Tissue"), with = FALSE])
  data.table::setorderv(ls_lbl_x, "Tissue")
  expect_is(plt_3, "gg")
  expect_equal(ggplot2::get_panel_scales(plt_3)$x$get_labels(), ls_lbl_x[["CellLineName"]])
  expect_length(unique(ggplot2::ggplot_build(plt_3)[["data"]][[2]][["fill"]]), 
                NROW(unique(ls_lbl_x[["Tissue"]])))
  expect_true(grepl(NROW(ls_comb), plt_3[["labels"]][["title"]]))
})

test_that("plot_boxplot_metric_combo_by_CLs works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_scores <- gDRutils::convert_se_assay_to_dt(se = se,
                                                assay_name = "scores")
  ls_comb <- unique(paste(dt_scores[["DrugName"]], "x", dt_scores[["DrugName_2"]]))
  
  plt_1 <- plot_boxplot_metric_combo_by_drugs(dt_scores) # default
  
  expect_is(plt_1, "gg")
  expect_equal(plt_1[["labels"]][["y"]], get_hm_title("hsa_score", "GR"))
  expect_true(grepl("celllines", plt_1[["labels"]][["title"]]))
  expect_length(plt_1[["layers"]], 3)
  expect_equal(sort(ggplot2::get_panel_scales(plt_1)$x$get_labels()), sort(ls_comb))
  
  plt_2 <- plot_boxplot_metric_combo_by_drugs(dt_scores,
                                              normalization_type = "RV",
                                              colors_vec = "gold")
  expect_is(plt_2, "gg")
  expect_equal(plt_2[["labels"]][["y"]], get_hm_title("hsa_score", "RV"))
  expect_length(plt_2[["layers"]], 3)
  expect_equal(unique(ggplot2::ggplot_build(plt_2)[["data"]][[2]][["fill"]]), "gold")
  
  plt_4 <- plot_boxplot_metric_combo_by_drugs(dt_scores,
                                              normalization_type = "RV",
                                              metric = "bliss_score",
                                              colors_vec = c("#0000FF", "#00FF00"))
  
  expect_is(plt_4, "gg")
  expect_equal(plt_4[["labels"]][["y"]], get_hm_title("bliss_score", "RV"))
  expect_equal(ggplot2::get_panel_scales(plt_4)$x$get_labels(), ls_comb)
  expect_equal(unique(ggplot2::ggplot_build(plt_4)[["data"]][[2]][["fill"]]), "#0000FF")
  expect_true(grepl(NROW(unique(dt_scores[["CellLineName"]])), plt_4[["labels"]][["title"]]))
})
