context("Test combo_plots")

test_that("heatmap_combo_metrics works as expected", {
  cl_name <- "cellline_FD"
  drug1_name <- "drug_011"
  drug2_name <- "drug_026"
  
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[["combination"]]
  
  plts_1 <- heatmap_combo_metrics(se, drug1_name, drug2_name, cl_name)
  expect_is(plts_1, "gg")
  expect_is(plts_1, "ggarrange")

  metric_growth <- "RV"
  plts_2 <- heatmap_combo_metrics(se, drug1_name, drug2_name, cl_name, metric_growth, as_panel = FALSE)
  expect_is(plts_2, "list")
  expect_equal(names(plts_2), c(names(gDRutils::get_combo_excess_field_names()), "iso_compare"))
  expect_true(all(vapply(seq_along(plts_2), 
                         function(i) grepl(metric_growth, plts_2[[i]]$labels$title), logical(1))))
  
  
  expect_error(heatmap_combo_metrics(se = mae,
                                     drug1_name = drug1_name,
                                     drug2_name = drug2_name, 
                                     cl_name = cl_name),
               "Assertion on 'se' failed: Must inherit from class 'SummarizedExperiment'")
  expect_error(heatmap_combo_metrics(se = se,
                                     drug1_name = "unknown_drug",
                                     drug2_name = drug2_name, 
                                     cl_name = cl_name),
               "Assertion on 'drug1_name' failed: Must be element of set")
  expect_error(heatmap_combo_metrics(se = se,
                                     drug1_name = drug1_name,
                                     drug2_name = drug2_name, 
                                     cl_name = cl_name,
                                     metric_growth = "AB"),
               "Assertion on 'metric_growth' failed: Must be element of set")
})
