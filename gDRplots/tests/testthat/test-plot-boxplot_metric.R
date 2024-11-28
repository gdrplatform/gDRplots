context("Test boxplot_metric")

test_that("plot_boxplot_sa_metric_by_CLs works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  
  plt_1 <- plot_boxplot_sa_metric_by_CLs(dt_metrics) # default
  
  expect_is(plt_1, "gg")
  expect_true(grepl("GR", plt_1[["labels"]][["y"]]))
  expect_true(grepl("xc50", plt_1[["labels"]][["y"]]))
  expect_length(plt_1[["layers"]], 3)
  expect_equal(sort(ggplot2::layer_scales(plt_1)$x$get_labels()),
               sort(unique(dt_metrics[["CellLineName"]])))
  
  plt_2 <- plot_boxplot_sa_metric_by_CLs(dt_metrics,
                                         normalization_type = "RV",
                                         metric = "x_AOC",
                                         colors_vec = "darkgreen")
  expect_is(plt_2, "gg")
  expect_true(grepl("RV", plt_2[["labels"]][["y"]]))
  expect_true(grepl("x_AOC", plt_2[["labels"]][["y"]]))
  expect_length(plt_2[["layers"]], 3)
  expect_equal(unique(ggplot2::ggplot_build(plt_2)[["data"]][[2]][["fill"]]), "darkgreen")
  
  plt_3 <- plot_boxplot_sa_metric_by_CLs(dt_metrics,
                                         colors_vec = c("blue", "yellow"))
  expect_is(plt_3, "gg")
  expect_equal(unique(ggplot2::ggplot_build(plt_3)[["data"]][[2]][["fill"]]), "blue")
})
