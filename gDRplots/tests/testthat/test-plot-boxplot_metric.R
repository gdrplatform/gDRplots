context("Test boxplot_metric")

test_that("plot_boxplot_sa_metric_by_CLs works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]

  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")

  plt_1 <- plot_boxplot_sa_metric_by_CLs(dt_metrics) # defaul
  
  expect_is(plt_1, "gg")
  expect_true(grepl("GR", plt_1[["labels"]][["y"]]))
  expect_true(grepl("xc50", plt_1[["labels"]][["y"]]))
  expect_length(plt_1[["layers"]], 3)
  expect_equal(sort(ggplot2::layer_scales(plt_1)$x$get_labels()),
               sort(unique(dt_metrics[["CellLineName"]])))
})
  