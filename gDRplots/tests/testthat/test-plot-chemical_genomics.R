context("Test plot-chemical_genomics")

# Load test data

metrics_data <- qs::qread(system.file("testdata/cgs_data.qs", package = "gDRplots"))

test_that("analyze_cgs works correctly", {
  # test with a single cell line and single metric
  results <- analyze_cgs(metrics_data, metrics = "xc50", cell_line = "CellLineName_1")
  expect_is(results, "list")
  expect_equal(names(results), "CellLineName_1")
  expect_is(results$CellLineName_1$fgsea$xc50, "data.table")
  expect_gt(nrow(results$CellLineName_1$fgsea$xc50), 0) # Check for results
  expect_is(results$CellLineName_1$metrics_diff, "data.table")
  expect_gt(nrow(results$CellLineName_1$metrics_diff), 0)
  expect_is(results$CellLineName_1$moa_list, "list")
  expect_gt(length(results$CellLineName_1$moa_list), 0)
  
  
  # test with multiple metrics
  results <- analyze_cgs(metrics_data, metrics = c("xc50", "x_max"), cell_line = "CellLineName_1")
  expect_equal(names(results$CellLineName_1$fgsea), c("xc50", "x_max"))
  expect_is(results$CellLineName_1$fgsea$xc50, "data.table")
  expect_is(results$CellLineName_1$fgsea$x_max, "data.table")
  
  # test normalization type
  results1 <- analyze_cgs(metrics_data, metrics = "xc50", cl_name = "CellLineName_1", normalization_type = "RV")
  results2 <- analyze_cgs(metrics_data, metrics = "xc50", cl_name = "CellLineName_1", normalization_type = "GR")
  
  expect_false(identical(results1$CellLineName_1$metrics_diff$xc50, results2$CellLineName_1$metrics_diff$xc50)) 
})

test_that("plot_cgs_ranking works correctly", {
  
  results <- analyze_cgs(metrics_data, metrics = c("xc50", "x_max"), cell_line = "CellLineName_1")
  plt <- plot_cgs_ranking(results, cell_line = "CellLineName_1", metric = "xc50")
  expect_is(plt, "ggplot")
  
  # rest with a different metric
  results <- analyze_cgs(metrics_data, metrics = c("xc50", "x_max"), cell_line = "CellLineName_1")
  plt <- plot_cgs_ranking(results, cell_line = "CellLineName_1", metric = "x_max")
  expect_is(plt, "ggplot")
  
  # test case where no significant GSEA results exist (using a dummy metric)
  results$CellLineName_1$fgsea$xc50$padj <- 1  # Set all padj values to 1 to simulate no significance
  plt <- plot_cgs_ranking(results, cl_name = "CellLineName_1", metric = "xc50") # Should still produce a plot
  expect_is(plt, "ggplot")
})
