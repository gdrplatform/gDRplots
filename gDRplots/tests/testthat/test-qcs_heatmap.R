context("Test qcs_heatmap")

test_that("heatmap_QCS works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix_small")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  response_metrics <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "Metrics")
  
  cdata <- SummarizedExperiment::colData(se)
  rdata <- SummarizedExperiment::rowData(se)
  
  plt_1 <-  heatmap_QCS(tab_response = response_metrics)
  expect_is(plt_1, "pheatmap")
  expect_equal(plt_1$gtable$grobs[[2]]$label, cdata[["CellLineName"]])
  expect_equal(plt_1$gtable$grobs[[3]]$label, rdata[["DrugName"]])
  expect_true(all(vapply(plt_1$gtable$grobs[[1]]$children[[1]]$gp$fill, isValidColor, logical(1))))
})


test_that("get_hm_title works as expected", {
  dataset_name <- "Dataset AB123"
  metric <- "x_max"
  metric_growth <- "RV"
  
  expect_equal(get_hm_title(dataset_name, metric, metric_growth), "Dataset AB123 (E_max)")
})
