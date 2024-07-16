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

test_that("prep_hm_limits works as expected", {
  vec <- c(1.0089, 0.9806, 0.1174, -0.1657, -0.2826)
  # GR smooth
  expect_equal(prep_hm_limits(vec), c(-0.2826, 1.0089))
  # RV smooth
  expect_equal(prep_hm_limits(vec, normalization_type = "RV"), c(0, 1.0089))

  # hsa_excess
  vec <- c(-0.1016, -0.0647, 0.0021, 0.0328, 0.6824)
  expect_equal(prep_hm_limits(vec, metric = "hsa_excess"), c(-0.25, 0.6824))
  # bliss_excess
  vec <- c(-0.2651, -0.1289, 0.0051, 0.0202, 0.0394)
  expect_equal(prep_hm_limits(vec, metric = "bliss_excess"), c(-0.2651, 0.25))
  
  expect_error(prep_hm_limits(LETTERS[1:5]),
               "Assertion on 'num_vec' failed: Must be of type 'numeric'")
  expect_error(prep_hm_limits(vec, metric = "str"),
               "Assertion on 'metric' failed: Must be element of set")
  expect_error(prep_hm_limits(vec, normalization_type = "str"),
               "Assertion on 'normalization_type' failed: Must be element of set")
})

test_that("transform_log_conc works as expected", {
  vec <- c(0, 0.003, 0.01, 0.03)
  result <- log10(vec)
  result[1] <- result[2] + (result[2] - result[3])
  expect_equal(transform_log_conc(vec), result)
  
  expect_error(transform_log_conc(LETTERS[seq_len(5)]), 
               "Assertion on 'conc_vec' failed: Must be of type 'numeric'")   
  
  expect_error(transform_log_conc(-1:2),
               "Assertion on 'conc_vec' failed: Element 1 is not >= 0")
  expect_error(transform_log_conc(c(0, 0.001)),
               "There are not enough values to handle 0.")
})
