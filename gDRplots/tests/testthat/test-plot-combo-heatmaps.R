context("Test combo_plots")

test_that("heatmap_combo_metrics works as expected", {
  cl_name <- "cellline_FD"
  drug1_name <- "drug_011"
  drug2_name <- "drug_026"
  
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
  dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")
  
  plts_1 <- heatmap_combo_metrics(dt_excess, dt_isobolograms, 
                                  drug1_name, drug2_name, cl_name)
  expect_is(plts_1, "gg")

  normalization_type <- "RV"
  plts_2 <- heatmap_combo_metrics(dt_excess, dt_isobolograms, 
                                  drug1_name, drug2_name, cl_name, 
                                  normalization_type, as_panel = FALSE)
  expect_is(plts_2, "list")
  expect_equal(names(plts_2), c(names(gDRutils::get_combo_excess_field_names()), "iso_compare"))
  expect_true(all(vapply(seq_along(plts_2), 
                         function(i) grepl(normalization_type, plts_2[[i]]$labels$title), logical(1))))
  
  
  expect_error(heatmap_combo_metrics(dt_excess = unlist(dt_excess), 
                                     dt_isobolograms = dt_isobolograms,
                                     drug1_name = drug1_name,
                                     drug2_name = drug2_name, 
                                     cl_name = cl_name),
               "Check on 'dt_excess' failed: Must be a data.table")
  expect_error(heatmap_combo_metrics(dt_excess = dt_excess, 
                                     dt_isobolograms = dt_isobolograms,
                                     drug1_name = "unknown_drug",
                                     drug2_name = drug2_name, 
                                     cl_name = cl_name),
               "Assertion on 'drug1_name' failed: Must be element of set")
  expect_error(heatmap_combo_metrics(dt_excess = dt_excess, 
                                     dt_isobolograms = dt_isobolograms,
                                     drug1_name = drug1_name,
                                     drug2_name = drug2_name, 
                                     cl_name = cl_name,
                                     normalization_type = "AB"),
               "Assertion on 'normalization_type' failed: Must be element of set")
})
