context("Test sa_plots")

test_that("plot_dose_response_combo works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
  
  cl_name <- "cellline_BC"
  drug1_name <- "drug_011"
  drug2_name <- "drug_021"
  cl_clid <- unique(dt_average[CellLineName == cl_name, ]$clid)
  
  plt_1 <- plot_dose_response_combo(dt_average = dt_average,
                                    drug1_name = drug1_name,
                                    drug2_name = drug2_name,
                                    cl_name = cl_name)
  expect_is(plt_1, "gg")
  expect_true(grepl("GR", plt_1[["labels"]][["y"]]))
  expect_length(plt_1[["layers"]], 3)
  expect_equal(plt_1[["labels"]][["title"]], sprintf("%s (%s)", cl_name, cl_clid))
  
  expect_error(plot_dose_response_combo(dt_average = unlist(dt_average),
                                        drug1_name = drug1_name,
                                        drug2_name = drug2_name,
                                        cl_name = cl_name),
               "Assertion on 'dt_average' failed: Must be a data.table")
  
  expect_error(plot_dose_response_combo(dt_average = dt_average,
                                        drug1_name = 1,
                                        drug2_name = drug2_name,
                                        cl_name = cl_name),
               "Assertion on 'drug1_name' failed: Must be element of set")
  
  expect_error(plot_dose_response_combo(dt_average = dt_average,
                                        drug1_name = drug1_name,
                                        drug2_name = "drug_XX",
                                        cl_name = cl_name),
               "Assertion on 'drug2_name' failed: Must be element of set")
  
  expect_error(plot_dose_response_combo(dt_average = dt_average,
                                        drug1_name = drug1_name,
                                        drug2_name = drug2_name,
                                        cl_name = 1),
               "Assertion on 'cl_name' failed: Must be of type 'string'")
  
  expect_error(plot_dose_response_combo(dt_average = dt_average,
                                        drug1_name = drug1_name,
                                        drug2_name = drug2_name,
                                        cl_name = "cellline_XX"),
               "Assertion on 'cl_name' failed: Must be element of set")
})

test_that("plot_dose_response_combo_qc_panel works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
  
  cl_name <- "cellline_IB"
  
  plt_1 <- plot_dose_response_combo_qc_panel(dt_average = dt_average,
                                             cl_name = cl_name)
  
  expect_is(plt_1, "gg")
  expect_true(grepl("GR", plt_1[["labels"]][["y"]]))
  expect_length(plt_1[["layers"]], 3)
  expect_true(grepl(cl_name, plt_1[["labels"]][["title"]]))
})