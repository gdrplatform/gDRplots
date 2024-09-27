context("Test helpers-prism")

conc <- gDRutils::get_env_identifiers("concentration")
drug_name <- gDRutils::get_env_identifiers("drug_name")
cellline_name <- gDRutils::get_env_identifiers("cellline_name")
meta_col <- c("rId", "cId", cellline_name)

test_that("prep_dt_response_metric_sa works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix_small")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Metrics")
  d_name <- "drug_004"
  
  res <- dt_metrics[get(drug_name) == d_name & normalization_type == "RV", ]
  
  dt_response <- prep_dt_response_metric_sa(dt_metrics, d_name) # default
  expect_is(dt_response, "data.table")
  expect_named(dt_response, c(meta_col, "RV_gDR_xc50"))
  expect_equal(NROW(dt_response), NROW(res))
  
  dt_response <- prep_dt_response_metric_sa(dt_metrics, d_name,
                                            normalization_type = "GR",
                                            metric = "x_mean") 
  expect_is(dt_response, "data.table")
  expect_named(dt_response, c(meta_col, "GR_gDR_x_mean"))
  
  sel_met <- c("xc50", "x_mean", "x_max")
  dt_response <-
    prep_dt_response_metric_sa(dt_metrics, d_name,
                               metric = sel_met)
  expect_is(dt_response, "data.table")
  expect_named(dt_response, c(meta_col, sprintf("RV_gDR_%s", sel_met)))
  
  # testing assertions
  expect_error(prep_dt_response_metric_sa(dt_metrics = unlist(dt_metrics),
                                          d_name = d_name),
               "Assertion on 'dt_metrics' failed: Must be a data.table")
  expect_error(prep_dt_response_metric_sa(dt_metrics = dt_metrics,
                                          d_name = 1),
               "Assertion on 'd_name' failed: Must be of type 'string'")
  expect_error(prep_dt_response_metric_sa(dt_metrics = dt_metrics,
                                          d_name = "drug_xx"),
               "Assertion on 'd_name' failed: Must be element of set")
  expect_error(prep_dt_response_metric_sa(dt_metrics = dt_metrics,
                                          d_name = d_name, 
                                          normalization_type = "str"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(prep_dt_response_metric_sa(dt_metrics = dt_metrics,
                                          d_name = d_name, 
                                          metric = 123),
               "Assertion on 'metric' failed: Must be of type 'character'")
  expect_error(prep_dt_response_metric_sa(dt_metrics = dt_metrics,
                                          d_name = d_name, 
                                          metric = "str"),
               "Assertion on 'metric' failed: Must be a subset of")
  expect_error(prep_dt_response_metric_sa(dt_metrics = dt_metrics,
                                          d_name = d_name, 
                                          fit_source = 123),
               "Assertion on 'fit_source' failed: Must be of type 'string'")
})

test_that("prep_dt_response_dose_sa works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix_small")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  dt_average <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Averaged")
  d_name <- "drug_004"
  ls_conc <- sprintf("%s", unique(dt_average[[conc]]))
  res <- data.table::dcast(
    dt_average[get(drug_name) == d_name & normalization_type == "RV", ],
    formula = get(cellline_name) ~ get(conc),
    metric = "x")
  
  dt_response <- prep_dt_response_dose_sa(dt_average, d_name) # default
  expect_is(dt_response, "data.table")
  expect_named(dt_response, c(meta_col, sprintf("RV_gDR_x_%s", ls_conc)))
  expect_equal(NROW(dt_response), NROW(res))
  
  dt_response <- prep_dt_response_dose_sa(dt_average, d_name,
                                          normalization_type = "GR",
                                          metric = "x_std")
  expect_is(dt_response, "data.table")
  expect_named(dt_response, c(meta_col, sprintf("GR_gDR_x_std_%s", ls_conc)))
  
  # testing assertions
  expect_error(prep_dt_response_dose_sa(dt_average = unlist(dt_average),
                                        d_name = d_name),
               "Assertion on 'dt_average' failed: Must be a data.table")
  expect_error(prep_dt_response_dose_sa(dt_average = dt_average,
                                        d_name = 1),
               "Assertion on 'd_name' failed: Must be of type 'string'")
  expect_error(prep_dt_response_dose_sa(dt_average = dt_average,
                                        d_name = "drug_xx"),
               "Assertion on 'd_name' failed: Must be element of set")
  expect_error(prep_dt_response_dose_sa(dt_average = dt_average,
                                        d_name = d_name, 
                                        normalization_type = "str"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(prep_dt_response_dose_sa(dt_average = dt_average,
                                        d_name = d_name, 
                                        metric = "str"),
               "Assertion on 'metric' failed: Must be element of set")
  expect_error(prep_dt_response_dose_sa(dt_average = dt_average,
                                        d_name = d_name, 
                                        fit_source = 123),
               "Assertion on 'fit_source' failed: Must be of type 'string'")
})

test_that("prep_dt_response_sa works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix_small")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  dt_average <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Averaged")
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Metrics")
  d_name <- "drug_004"
  
  sel_met <- c("xc50", "x_mean", "x_max")
  res_met <- dt_metrics[get(drug_name) == d_name & normalization_type == "RV", ]
  ls_conc <- sprintf("%s", unique(dt_average[[conc]]))
  res_dos <- data.table::dcast(
    dt_average[get(drug_name) == d_name & normalization_type == "RV", ],
    formula = get(cellline_name) ~ get(conc),
    metric = "x")
  
  dt_response_sa <- prep_dt_response_sa(dt_average, dt_metrics, d_name) # default
  expect_is(dt_response_sa, "data.table")
  expect_named(dt_response_sa, 
               c(meta_col, sprintf("RV_gDR_%s", sel_met), sprintf("RV_gDR_x_%s", ls_conc)))
  expect_equal(NROW(dt_response_sa), NROW(res_met))
  expect_equal(NROW(dt_response_sa), NROW(res_dos))
  
  # testing assertions
  expect_error(prep_dt_response_sa(dt_average = unlist(dt_average),
                                   dt_metrics = dt_metrics,
                                   d_name = d_name),
               "Assertion on 'dt_average' failed: Must be a data.table")
  expect_error(prep_dt_response_sa(dt_average = dt_average,
                                   dt_metrics = unlist(dt_metrics),
                                   d_name = d_name),
               "Assertion on 'dt_metrics' failed: Must be a data.table")
  expect_error(prep_dt_response_sa(dt_average = dt_average,
                                   dt_metrics = dt_metrics,
                                   d_name = 1),
               "Assertion on 'd_name' failed: Must be of type 'string'")
  expect_error(prep_dt_response_sa(dt_average = dt_average,
                                   dt_metrics = dt_metrics,
                                   d_name = "drug_xx"),
               "Assertion on 'd_name' failed: Must be element of set")
  expect_error(prep_dt_response_sa(dt_average = dt_average,
                                   dt_metrics = dt_metrics,
                                   d_name = d_name, 
                                   normalization_type = "str"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(prep_dt_response_sa(dt_average = dt_average,
                                   dt_metrics = dt_metrics,
                                   d_name = d_name, 
                                   fit_source = 123),
               "Assertion on 'fit_source' failed: Must be of type 'string'")
})
