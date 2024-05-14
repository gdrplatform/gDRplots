context("Test SA_plots")

test_that("prepareCurves works as expected", {
  mae <- gDRutils::get_synthetic_data("small")
  SE <- mae[[1]]
  iR <- rownames(SE)[1]
  grouping <- "cId"
  dt_metrics <- gDRutils::convert_se_assay_to_dt(SE[iR], "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(SE[iR], "Averaged")

  plt <- grob_SA(dt_metrics = dt_metrics,
                 dt_average = dt_average,
                 grouping = grouping)
  expect_is(plt, "gg")
  expect_true(grepl("GR", plt[["labels"]][["y"]]))
  expect_length(plt[["layers"]], 4)
  
  normalization_type <- "RV"
  plt <- grob_SA(dt_metrics = dt_metrics,
                 dt_average = dt_average,
                 grouping = grouping,
                 normalization_type = normalization_type,
                 plot_fit_flag = FALSE)
  expect_is(plt, "gg")
  expect_true(grepl(normalization_type, plt[["labels"]][["y"]]))
  expect_length(plt[["layers"]], 3)
})
