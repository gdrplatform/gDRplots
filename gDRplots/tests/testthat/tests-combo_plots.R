context("Test combo_plots")

test_that("prepareCurves works as expected", {
  cellline_name <- "cellline_FD"
  drug1_name <- "drug_011"
  drug2_name <- "drug_026"
  norm_type <- "RV"
  
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[["combination"]]
  
  plts_1 <- gDR_combo_plot(se, drug1_name, drug2_name, cellline_name)
  expect_is(plts_1, "list")
  expect_equal(names(plts_1), c(names(gDRutils::get_combo_excess_field_names()), "iso_compare"))
  expect_true(all(vapply(seq_along(plts_1), 
                         function(i) grepl("GR", plts_1[[i]]$labels$title), logical(1))))

  normalization_type <- "RV"
  plts_2 <- gDR_combo_plot(se, drug1_name, drug2_name, cellline_name, normalization_type)
  expect_is(plts_2, "list")
  expect_equal(names(plts_2), c(names(gDRutils::get_combo_excess_field_names()), "iso_compare"))
  expect_true(all(vapply(seq_along(plts_2), 
                         function(i) grepl(normalization_type, plts_2[[i]]$labels$title), logical(1))))
  
  
  expect_error(gDR_combo_plot(se = mae,
                              drug1_name = drug1_name,
                              drug2_name = drug2_name, 
                              cellline_name = cellline_name),
               "Assertion on 'se' failed: Must inherit from class 'SummarizedExperiment'")
  expect_error(gDR_combo_plot(se = se,
                              drug1_name = "unknown_drug",
                              drug2_name = drug2_name, 
                              cellline_name = cellline_name),
               "Assertion on 'drug1_name' failed: Must be element of set")
})
