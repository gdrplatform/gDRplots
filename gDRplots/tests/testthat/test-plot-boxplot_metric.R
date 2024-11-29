context("Test boxplot_metric")

test_that("plot_boxplot_metric_sa_by_CLs works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  
  plt_1 <- plot_boxplot_metric_sa_by_CLs(dt_metrics) # default
  
  expect_is(plt_1, "gg")
  expect_true(grepl("GR", plt_1[["labels"]][["y"]]))
  expect_true(grepl("xc50", plt_1[["labels"]][["y"]]))
  expect_length(plt_1[["layers"]], 3)
  expect_equal(sort(ggplot2::layer_scales(plt_1)$x$get_labels()),
               sort(unique(dt_metrics[["CellLineName"]])))
  
  plt_2 <- plot_boxplot_metric_sa_by_CLs(dt_metrics,
                                         normalization_type = "RV",
                                         metric = "x_AOC",
                                         colors_vec = "darkgreen")
  expect_is(plt_2, "gg")
  expect_true(grepl("RV", plt_2[["labels"]][["y"]]))
  expect_true(grepl("x_AOC", plt_2[["labels"]][["y"]]))
  expect_length(plt_2[["layers"]], 3)
  expect_equal(unique(ggplot2::ggplot_build(plt_2)[["data"]][[2]][["fill"]]), "darkgreen")
  
  plt_3 <- plot_boxplot_metric_sa_by_CLs(dt_metrics,
                                         colors_vec = c("blue", "yellow"))
  expect_is(plt_3, "gg")
  expect_equal(unique(ggplot2::ggplot_build(plt_3)[["data"]][[2]][["fill"]]), "blue")
  
  plt_4 <- plot_boxplot_metric_sa_by_CLs(dt_metrics,
                                         metric = "x_max",
                                         grouped_flag = TRUE)
  
  ls_lbl_x <- unique(dt_metrics[, c("CellLineName", "Tissue"), with = FALSE])
  data.table::setorderv(ls_lbl_x, "Tissue")
  expect_is(plt_4, "gg")
  expect_equal(ggplot2::layer_scales(plt_4)$x$get_labels(), ls_lbl_x[["CellLineName"]])
  expect_length(unique(ggplot2::ggplot_build(plt_4)[["data"]][[2]][["fill"]]), 
                NROW(unique(ls_lbl_x[["Tissue"]])))
  expect_true(grepl(NROW(unique(dt_metrics[["DrugName"]])), plt_4[["labels"]][["title"]]))
  
  expect_error(plot_boxplot_metric_sa_by_CLs(dt_metrics = unlist(dt_metrics)),
               "Assertion on 'dt_metrics' failed: Must be a data.table")
  expect_error(plot_boxplot_metric_sa_by_CLs(dt_metrics = dt_metrics,
                                             normalization_type = "XX"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(plot_boxplot_metric_sa_by_CLs(dt_metrics = dt_metrics,
                                             metric = "xxx"),
               "Assertion on 'metric' failed: Must be element of set")
  expect_error(plot_boxplot_metric_sa_by_CLs(dt_metrics = dt_metrics,
                                             fit_source = 1),
               "Assertion on 'fit_source' failed: Must be of type 'string'")
  expect_error(plot_boxplot_metric_sa_by_CLs(dt_metrics = dt_metrics,
                                             grouped_flag = "yes"),
               "Assertion on 'grouped_flag' failed: Must be of type 'logical flag'")
  expect_error(plot_boxplot_metric_sa_by_CLs(dt_metrics = dt_metrics,
                                             colors_vec = 1:3),
               "Assertion on 'colors_vec' failed: Must be of type 'character'")
})

test_that("plot_boxplot_metric_sa_by_drugs works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  
  plt_1 <- plot_boxplot_metric_sa_by_drugs(dt_metrics) # default
  
  expect_is(plt_1, "gg")
  expect_true(grepl("GR", plt_1[["labels"]][["y"]]))
  expect_true(grepl("xc50", plt_1[["labels"]][["y"]]))
  expect_true(grepl("cellline", plt_1[["labels"]][["title"]]))
  expect_length(plt_1[["layers"]], 3)
  expect_equal(sort(ggplot2::layer_scales(plt_1)$x$get_labels()),
               sort(unique(dt_metrics[["DrugName"]])))
  
  plt_2 <- plot_boxplot_metric_sa_by_drugs(dt_metrics,
                                           normalization_type = "RV",
                                           metric = "x_AOC",
                                           colors_vec = "gold")
  expect_is(plt_2, "gg")
  expect_true(grepl("RV", plt_2[["labels"]][["y"]]))
  expect_true(grepl("x_AOC", plt_2[["labels"]][["y"]]))
  expect_length(plt_2[["layers"]], 3)
  expect_equal(unique(ggplot2::ggplot_build(plt_2)[["data"]][[2]][["fill"]]), "gold")
  
  plt_3 <- plot_boxplot_metric_sa_by_drugs(dt_metrics,
                                           colors_vec = c("darkred", "yellow"))
  expect_is(plt_3, "gg")
  expect_equal(unique(ggplot2::ggplot_build(plt_3)[["data"]][[2]][["fill"]]), "darkred")
  
  plt_4 <- plot_boxplot_metric_sa_by_drugs(dt_metrics,
                                           metric = "x_max",
                                           grouped_flag = TRUE,
                                           colors_vec = c("#0000FF", "#00FF00"))
  
  ls_lbl_x <- unique(dt_metrics[, c("DrugName", "drug_moa"), with = FALSE])
  data.table::setorderv(ls_lbl_x, "drug_moa")
  expect_is(plt_4, "gg")
  expect_equal(ggplot2::layer_scales(plt_4)$x$get_labels(), ls_lbl_x[["DrugName"]])
  expect_length(unique(ggplot2::ggplot_build(plt_4)[["data"]][[2]][["fill"]]), 
                NROW(unique(ls_lbl_x[["drug_moa"]])))
  expect_true(all(c("#0000FF", "#00FF00") %in% unique(ggplot2::ggplot_build(plt_4)[["data"]][[2]][["fill"]])))
})
