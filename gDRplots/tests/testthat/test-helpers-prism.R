context("Test helpers-prism")

conc <- gDRutils::get_env_identifiers("concentration")
drug_name <- gDRutils::get_env_identifiers("drug_name")
drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
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
  dt_response <- prep_dt_response_metric_sa(dt_metrics, d_name,
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

test_that("prep_dt_response_scores works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix_small")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_scores <- gDRutils::convert_se_assay_to_dt(se = se,
                                                assay_name = "scores")
  d_name <- "drug_004"
  d_name2 <- "drug_026"
  
  res <- 
    dt_scores[get(drug_name) == d_name & get(drug_name_2) == d_name2 & normalization_type == "RV", ]
  
  dt_response <- prep_dt_response_scores(dt_scores, d_name, d_name2) # default
  expect_is(dt_response, "data.table")
  expect_named(dt_response, c(meta_col, "RV_gDR_hsa_score"))
  expect_equal(NROW(dt_response), NROW(res))
  
  dt_response <- prep_dt_response_scores(dt_scores, d_name, d_name2,
                                         normalization_type = "GR",
                                         metric = "bliss_score") 
  expect_is(dt_response, "data.table")
  expect_named(dt_response, c(meta_col, "GR_gDR_bliss_score"))
  
  sel_met <- c("hsa_score", "bliss_score")
  dt_response <- prep_dt_response_scores(dt_scores, d_name, d_name2,
                                         metric = sel_met)
  expect_is(dt_response, "data.table")
  expect_named(dt_response, c(meta_col, sprintf("RV_gDR_%s", sel_met)))
  
  # testing assertions
  expect_error(prep_dt_response_scores(dt_scores = unlist(dt_scores),
                                       d_name = d_name,
                                       d_name2 = d_name2),
               "Assertion on 'dt_scores' failed: Must be a data.table")
  expect_error(prep_dt_response_scores(dt_scores = dt_scores,
                                       d_name = 1,
                                       d_name2 = d_name2),
               "Assertion on 'd_name' failed: Must be of type 'string'")
  expect_error(prep_dt_response_scores(dt_scores = dt_scores,
                                       d_name = "drug_xx",
                                       d_name2 = d_name2),
               "Assertion on 'd_name' failed: Must be element of set")
  expect_error(prep_dt_response_scores(dt_scores = dt_scores,
                                       d_name = d_name,
                                       d_name2 = 1),
               "Assertion on 'd_name2' failed: Must be of type 'string'")
  expect_error(prep_dt_response_scores(dt_scores = dt_scores,
                                       d_name = d_name,
                                       d_name2 = "drug_yy"),
               "Assertion on 'd_name2' failed: Must be element of set")
  expect_error(prep_dt_response_scores(dt_scores = dt_scores,
                                       d_name = d_name,
                                       d_name2 = d_name2,
                                       normalization_type = "str"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(prep_dt_response_scores(dt_scores = dt_scores,
                                       d_name = d_name,
                                       d_name2 = d_name2,
                                       metric = 123),
               "Assertion on 'metric' failed: Must be of type 'character'")
  expect_error(prep_dt_response_scores(dt_scores = dt_scores,
                                       d_name = d_name,
                                       d_name2 = d_name2,
                                       metric = "str"),
               "Assertion on 'metric' failed: Must be a subset of")
  expect_error(prep_dt_response_scores(dt_scores = dt_scores,
                                       d_name = d_name,
                                       d_name2 = d_name2,
                                       fit_source = 123),
               "Assertion on 'fit_source' failed: Must be of type 'string'")
})

test_that("prep_dt_response_metric_diff works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix_small")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Metrics")
  d_name <- "drug_004"
  d_name2 <- "drug_026"
  
  subset <- 
    dt_metrics[get(drug_name) == d_name & get(drug_name_2) == d_name2 & 
                 normalization_type == "RV", ][!is.na(cotrt_value)]
  # Inf -> 10^dt_metrics[["maxlog10Concentration"]] # nolint
  subset[, xc50 := ifelse(is.infinite(xc50), 10^maxlog10Concentration, xc50)]
  
  meta_col_str <- paste(c("rId", "cId", gDRutils::get_env_identifiers("cellline_name")), collapse = " + ")
  
  dcast_formula <- as.formula(paste(meta_col_str, "~ cotrt_value + source"))
  
  
  res <- data.table::dcast(subset, formula = dcast_formula, 
                           value.var = "xc50")
  ls_col_diff <- setdiff(colnames(res), meta_col)
  ls_col_diff_fin <- ls_col_diff[!grepl("0_", ls_col_diff)]
  sel_met <- c("xc50", "x_mean", "x_max")
  comb <- c(c("cotrt_zero", "cotrt", "cotrt_diff"))
  
  dt_response <- prep_dt_response_metric_diff(dt_metrics, d_name, d_name2) # default
  expect_is(dt_response, "data.table")
  expect_true(all(
    names(dt_response) %in% 
      c(meta_col, drug_name, drug_name_2,
        do.call(paste0, expand.grid(sprintf("RV_gDR_xc50_%s_", comb), ls_col_diff_fin)))
  ))
  expect_equal(NROW(dt_response), NROW(res))
  expect_equal(dt_response[["RV_gDR_xc50_cotrt_zero_0.001_col_fittings"]], 
               res[["0_col_fittings"]])
  expect_equal(dt_response[["RV_gDR_xc50_cotrt_0.001_col_fittings"]], 
               res[["0.001_col_fittings"]])
  expect_equal(dt_response[["RV_gDR_xc50_cotrt_0.00316_row_fittings"]], 
               res[["0.00316_row_fittings"]]) 
  
  dt_response <-
    prep_dt_response_metric_diff(dt_metrics, d_name, d_name2,
                                 normalization_type = "GR",
                                 metric = c("xc50", "x_mean", "x_max"))
  
  expect_is(dt_response, "data.table")
  expect_is(dt_response, "data.table")
  expect_true(all(
    names(dt_response) %in% 
      c(meta_col, drug_name, drug_name_2,
        do.call(paste0, expand.grid(sprintf("GR_gDR_%s", sel_met), sprintf("_%s_", comb), ls_col_diff_fin)))
  ))
  expect_equal(NROW(dt_response), NROW(res))
  
})
#nolint start
# test_that("prep_dt_depmap_feat works as expected", {
#   # TODO in GDR-2710
# })
# 
# test_that("prep_dt_depmap_meta works as expected", {
#   # TODO in GDR-2710
# })
# 
# test_that("prep_dt_assoc works as expected", {
#   # TODO in GDR-2710
# })
#nolint end
