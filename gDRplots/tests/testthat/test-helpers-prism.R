context("Test helpers-prism")

conc <- gDRutils::get_env_identifiers("concentration")
drug_name <- gDRutils::get_env_identifiers("drug_name")
drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
cellline_name <- gDRutils::get_env_identifiers("cellline_name")
drug_moa_2 <- gDRutils::get_env_identifiers("drug_moa2")
meta_col <- c("rId", "cId", cellline_name)

test_that("prep_dt_response_metric_sa works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix_small")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Averaged")
  dt_metrics_capped <-
    gDRutils::cap_assay_infinities(conc_assay_dt = dt_average,
                                   assay_dt = dt_metrics,
                                   experiment_name = gDRutils::get_supported_experiments("sa"),
                                   capping_fold = 5)
  d_name <- "drug_004"
  
  res <- dt_metrics[get(drug_name) == d_name & normalization_type == "RV", ]
  
  dt_response <- prep_dt_response_metric_sa(dt_metrics, d_name) # default
  expect_is(dt_response, "data.table")
  expect_named(dt_response, c(meta_col, "RV_gDR_log10_xc50"))
  expect_equal(NROW(dt_response), NROW(res))
  
  d_name <- "drug_021"
  dt_response <- prep_dt_response_metric_sa(dt_metrics, d_name, metric = "xc50")
  dt_response_capped <- prep_dt_response_metric_sa(dt_metrics_capped, d_name, metric = "xc50")
  expect_identical(dt_response[, c(meta_col)], dt_response_capped[, c(meta_col)])
  expect_true(all(is.infinite(dt_response$RV_gDR_log10_xc50)))
  expect_false(all(is.infinite(dt_response_capped$RV_gDR_log10_xc50)))
  
  dt_response <- prep_dt_response_metric_sa(dt_metrics, d_name,
                                            normalization_type = "GR",
                                            metric = c("xc50", "x_mean"))
  expect_is(dt_response, "data.table")
  expect_named(dt_response, c(meta_col, "GR_gDR_log10_xc50", "GR_gDR_x_mean"))
  
  sel_met <- c("xc50", "x_mean", "x_max")
  dt_response <- prep_dt_response_metric_sa(dt_metrics, d_name,
                                            metric = sel_met)
  expect_is(dt_response, "data.table")
  expect_named(dt_response, c(meta_col, "RV_gDR_log10_xc50", sprintf("RV_gDR_%s", sel_met[2:3])))
  
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
  meta_col_str <- paste(c("rId", "cId", gDRutils::get_env_identifiers("cellline_name")), collapse = " + ")
  dcast_formula <- as.formula(paste(meta_col_str, "~ cotrt_value + dilution_drug"))
  
  res <- data.table::dcast(subset, formula = dcast_formula, 
                           value.var = "xc50")
  ls_col_diff <- setdiff(colnames(res), meta_col)
  ls_col_diff_fin <- ls_col_diff[!grepl("0_", ls_col_diff)]
  sel_met <- c("xc50", "x_mean", "x_max")
  comb <- c(c("cotrt_zero", "cotrt", "cotrt_diff"))
  
  # scenario: default
  dt_response <- 
    prep_dt_response_metric_diff(dt_metrics = dt_metrics,
                                 d_name = d_name,
                                 d_name2 = d_name2) # default
  
  expect_is(dt_response, "data.table")
  expect_true(all(
    names(dt_response) %in% 
      c(meta_col, drug_name, drug_name_2,
        do.call(paste0, expand.grid(sprintf("RV_gDR_log10_xc50_%s_", comb), ls_col_diff_fin)))
  ))
  expect_equal(NROW(dt_response), NROW(res))
  expect_equal(dt_response[["RV_gDR_log10_xc50_cotrt_zero_0.001_drug_1"]], 
               log10(res[["0_drug_1"]]))
  expect_equal(dt_response[["RV_gDR_log10_xc50_cotrt_0.001_drug_1"]], 
               log10(res[["0.001_drug_1"]]))
  expect_equal(dt_response[["RV_gDR_log10_xc50_cotrt_0.00316_drug_2"]], 
               log10(res[["0.00316_drug_2"]]))
  
  # scenario: "GR" and list of metrics
  dt_response_GR <-
    prep_dt_response_metric_diff(dt_metrics = dt_metrics,
                                 d_name = d_name,
                                 d_name2 = d_name2,
                                 normalization_type = "GR",
                                 metric = c("x_mean", "x_max"))
  
  expect_is(dt_response_GR, "data.table")
  expect_true(all(
    names(dt_response_GR) %in% 
      c(meta_col, drug_name, drug_name_2,
        do.call(paste0, expand.grid(sprintf("GR_gDR_%s", sel_met), sprintf("_%s_", comb), ls_col_diff_fin)))
  ))
  expect_equal(NROW(dt_response_GR), NROW(res))
  
  # scenario: capped values
  dt_average <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Averaged")
  dt_metrics_capped <- 
    gDRutils::cap_assay_infinities(conc_assay_dt = dt_average,
                                   assay_dt = dt_metrics,
                                   experiment_name = gDRutils::get_supported_experiments("combo"))
  
  dt_response_cap <- prep_dt_response_metric_diff(dt_metrics = dt_metrics_capped, 
                                                  d_name = d_name, 
                                                  d_name2 = d_name2) # default
  expect_is(dt_response_cap, "data.table")
  expect_true(all(
    names(dt_response_cap) %in% 
      c(meta_col, drug_name, drug_name_2,
        do.call(paste0, expand.grid(sprintf("RV_gDR_log10_xc50_%s_", comb), ls_col_diff_fin)))
  ))
  expect_equal(NROW(dt_response_cap), NROW(res))
  expect_equal(dt_response_cap[["RV_gDR_log10_xc50_cotrt_zero_0.001_drug_1"]], 
               log10(res[["0_drug_1"]]))
  expect_equal(dt_response_cap[["RV_gDR_log10_xc50_cotrt_0.001_drug_1"]], 
               log10(res[["0.001_drug_1"]]))
  
  ls_col_inf <- names(dt_response)[
    vapply(names(dt_response), function(nm) all(is.infinite(dt_response[[nm]])), logical(1))]
  expect_false(all(
    vapply(ls_col_inf, function(nm) all(is.infinite(dt_response_cap[[nm]])), logical(1))))
  ls_col_na <- names(dt_response)[
    vapply(names(dt_response), function(nm) all(is.na(dt_response[[nm]])), logical(1))]
  expect_false(all(
    vapply(ls_col_inf, function(nm) all(is.na(dt_response_cap[[nm]])), logical(1))))
  ls_col_equal <- names(dt_response)[!names(dt_response) %in% c(ls_col_inf, ls_col_na)]
  ls_col_equal <- ls_col_equal[!grepl("_drug_2", ls_col_equal)]
  expect_equal(dt_response_cap[, ls_col_equal, with = FALSE], 
               dt_response[, ls_col_equal, with = FALSE])
  
  # scenario: additional column
  dt_response_addcol <- prep_dt_response_metric_diff(dt_metrics = dt_metrics_capped, 
                                                     d_name = d_name, 
                                                     d_name2 = d_name2,
                                                     additional_cols = drug_moa_2)
  expect_is(dt_response_addcol, "data.table")
  expect_equal(setdiff(names(dt_response_addcol), names(dt_response)), drug_moa_2)
  expect_equal(dt_response_addcol[, -drug_moa_2, with = FALSE], dt_response_cap)
  
  # scenario: cell lines diff (cgs case)
  dt_response_cl_diff <- prep_dt_response_metric_diff(dt_metrics = dt_metrics_capped,
                                                      d_name = NULL,
                                                      d_name2 = NULL,
                                                      cellline1 = "cellline_GB",
                                                      cellline2 = "cellline_HB",
                                                      additional_cols = drug_moa_2)
  expect_is(dt_response_cl_diff, "data.table")
  expect_true("xc50_cellline_diff" %in% names(dt_response_cl_diff))
  expect_equal(sum(endsWith(names(dt_response_cl_diff), "c1")), 1)
  expect_equal(sum(endsWith(names(dt_response_cl_diff), "c2")), 1)
  
  
  # testing assertions
  expect_error(prep_dt_response_metric_diff(dt_metrics = unlist(dt_metrics),
                                            d_name = d_name,
                                            d_name2 = d_name2),
               "Assertion on 'dt_metrics' failed: Must be a data.table")
  expect_error(prep_dt_response_metric_diff(dt_metrics = dt_metrics,
                                            d_name = 1,
                                            d_name2 = d_name2),
               "Assertion on 'd_name' failed: Must be of type 'string'")
  expect_error(prep_dt_response_metric_diff(dt_metrics = dt_metrics,
                                            d_name = "drug_xx",
                                            d_name2 = d_name2),
               "Assertion on 'd_name' failed: Must be element of set")
  expect_error(prep_dt_response_metric_diff(dt_metrics = dt_metrics,
                                            d_name = d_name,
                                            d_name2 = d_name2,
                                            normalization_type = "str"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(prep_dt_response_metric_diff(dt_metrics = dt_metrics,
                                            d_name = d_name, 
                                            metric = 123),
               "Assertion on 'metric' failed: Must be of type 'character'")
  expect_error(prep_dt_response_metric_diff(dt_metrics = dt_metrics,
                                            d_name = d_name,
                                            d_name2 = d_name2,
                                            metric = "str"),
               "Assertion on 'metric' failed: Must be a subset of")
  expect_error(prep_dt_response_metric_diff(dt_metrics = dt_metrics,
                                            d_name = d_name,
                                            d_name2 = d_name2,
                                            fit_source = 123),
               "Assertion on 'fit_source' failed: Must be of type 'string'")
  expect_error(prep_dt_response_metric_diff(dt_metrics = dt_metrics,
                                            d_name = d_name,
                                            d_name2 = d_name2,
                                            additional_cols = 123),
               "Assertion on 'additional_cols' failed: Must be element of set")
})


test_that("prep_dt_depmap_meta works as expected", {
  id_col <- c("ModelID", "CCLEName") 
  test_meta_data_path <- system.file("testdata/Model.csv", package = "gDRplots")
  tab_model <- data.table::fread(test_meta_data_path)
  
  obj_meta_1 <- prep_dt_depmap_meta(meta_data_path = test_meta_data_path) # default
  expect_is(obj_meta_1, "list")
  expect_is(obj_meta_1[[1]], "data.table")
  expect_equal(names(obj_meta_1), c("dt_depmap", "selected_feat_meta_col"))
  expect_equal(obj_meta_1$selected_feat_meta_col, "PatientRace")
  expect_true(all(id_col %in% names(obj_meta_1$dt_depmap)))
  expect_true(all(unique(tab_model[["PatientRace"]]) %in% names(obj_meta_1$dt_depmap)))
  expect_true(all(vapply(obj_meta_1$dt_depmap[, .SD, .SDcols = -id_col], is.numeric, logical(1))))
  
  meta_2 <- "OncotreeLineage"
  obj_meta_2 <- prep_dt_depmap_meta(test_meta_data_path,
                                    metadata_col = meta_2)
  expect_is(obj_meta_2, "list")
  expect_equal(obj_meta_2$selected_feat_meta_col, meta_2)
  expect_true(all(unique(tab_model[[meta_2]]) %in% names(obj_meta_2$dt_depmap)))
  expect_true(all(vapply(obj_meta_2$dt_depmap[, .SD, .SDcols = -id_col], is.numeric, logical(1))))
  
  # scenario: origin column is numeric type
  meta_3 <- "Age"
  obj_meta_3 <- prep_dt_depmap_meta(meta_data_path = test_meta_data_path,
                                    metadata_col = meta_3)
  expect_is(obj_meta_3, "list")
  expect_equal(obj_meta_3$selected_feat_meta_col, meta_3)
  ls_col_3 <- vapply(unique(tab_model[[meta_3]]), change_NA_into_char, FUN.VALUE = character(1))
  expect_true(all(c(id_col, ls_col_3) %in% names(obj_meta_3$dt_depmap)))
  expect_true(all(vapply(obj_meta_3$dt_depmap[, .SD, .SDcols = -id_col], is.numeric, logical(1))))
  
  # scenario: origin column is logial type
  meta_4 <- "SourceDetail"
  obj_meta_4 <- prep_dt_depmap_meta(meta_data_path = test_meta_data_path,
                                    metadata_col = meta_4)
  expect_is(obj_meta_4, "list")
  expect_equal(obj_meta_4$selected_feat_meta_col, meta_4)
  ls_col_4 <- vapply(unique(tab_model[[meta_4]]), change_NA_into_char, FUN.VALUE = character(1))
  expect_true(all(c(id_col, ls_col_4) %in% names(obj_meta_4$dt_depmap)))
  expect_true(all(vapply(obj_meta_4$dt_depmap[, .SD, .SDcols = -id_col], is.numeric, logical(1))))
  
  # scenario: origin column is factor type
  meta_5 <- "Sex"
  obj_meta_5 <- prep_dt_depmap_meta(meta_data_path = test_meta_data_path,
                                    metadata_col = meta_5)
  expect_is(obj_meta_5, "list")
  expect_equal(obj_meta_5$selected_feat_meta_col, meta_5)
  ls_col_5 <- vapply(unique(tab_model[[meta_5]]), change_NA_into_char, FUN.VALUE = character(1))
  expect_true(all(c(id_col, ls_col_5) %in% names(obj_meta_5$dt_depmap)))
  expect_true(all(vapply(obj_meta_5$dt_depmap[, .SD, .SDcols = -id_col], is.numeric, logical(1))))
  
  # scenario: empty string in values
  meta_6 <- "TreatmentStatus"
  obj_meta_6 <- prep_dt_depmap_meta(meta_data_path = test_meta_data_path,
                                    metadata_col = meta_6)
  expect_is(obj_meta_6, "list")
  expect_equal(obj_meta_6$selected_feat_meta_col, meta_6)
  ls_col_6 <- vapply(unique(tab_model[[meta_6]]), change_NA_into_char, FUN.VALUE = character(1))
  expect_true(all(c(id_col, ls_col_6[ls_col_6 != ""]) %in% names(obj_meta_6$dt_depmap)))
  expect_true(all(vapply(obj_meta_6$dt_depmap[, .SD, .SDcols = -id_col], is.numeric, logical(1))))
  
  # scenario: id columns selected
  meta_7 <- id_col[1]
  obj_meta_7 <- prep_dt_depmap_meta(meta_data_path = test_meta_data_path,
                                    metadata_col = meta_7)
  expect_is(obj_meta_7, "list")
  expect_equal(obj_meta_7$selected_feat_meta_col, meta_7)
  expect_true(all(id_col %in% names(obj_meta_7$dt_depmap)))
  expect_equal(NCOL(obj_meta_7$dt_depmap), NROW(id_col))
  
  expect_error(prep_dt_depmap_meta(123), 
               "Assertion on 'meta_data_path' failed: Must be of type 'string'")
  expect_error(prep_dt_depmap_meta("testdata/meta_data.qs"), 
               "Assertion on 'File ext must be csv' failed: Must be TRUE")
  expect_error(prep_dt_depmap_meta("testdata/meta_data.csv"), 
               "Assertion on 'meta_data_path' failed: File does not exist")
  expect_error(prep_dt_depmap_meta(test_meta_data_path, metadata_col = 123),
               "Assertion on 'metadata_col' failed: Must be of type 'string'")
  expect_error(prep_dt_depmap_meta(test_meta_data_path, metadata_col = "some_meta"),
               "failed: Must be a subset of")
})

test_that("prep_dt_depmap_feat works as expected", {
  id_col <- c("ModelID", "CCLEName") 
  test_meta_data_path <- system.file("testdata/Model.csv", package = "gDRplots")
  tab_model <- data.table::fread(test_meta_data_path)
  
  test_feat_data_path <- system.file("testdata", package = "gDRplots")
  
  obj_feat_1 <- prep_dt_depmap_feat(feat_data_path = test_feat_data_path,
                                    meta_data_path = test_meta_data_path) # default
  expect_is(obj_feat_1, "list")
  expect_is(obj_feat_1[[1]], "data.table")
  expect_equal(names(obj_feat_1), c("dt_depmap", "selected_feat_meta_col"))
  expect_equal(obj_feat_1$selected_feat_meta_col, "CRISPRGeneEffect")
  expect_true(all(id_col %in% names(obj_feat_1$dt_depmap)))
  expect_true(all(vapply(obj_feat_1$dt_depmap[, .SD, .SDcols = -id_col], is.numeric, logical(1))))
  
  obj_feat_2 <- prep_dt_depmap_feat(feat_data_path = test_feat_data_path,
                                    meta_data_path = test_meta_data_path,
                                    feature_set = "OmicsSomaticMutationsMatrixHotspot")
  expect_is(obj_feat_2, "list")
  expect_is(obj_feat_2[[1]], "data.table")
  expect_equal(names(obj_feat_2), c("dt_depmap", "selected_feat_meta_col"))
  expect_equal(obj_feat_2$selected_feat_meta_col, "OmicsSomaticMutationsMatrixHotspot")
  expect_true(all(id_col %in% names(obj_feat_2$dt_depmap)))
  expect_true(all(vapply(obj_feat_2$dt_depmap[, .SD, .SDcols = -id_col], is.numeric, logical(1))))
  
  # sceanraio: not supported feature
  expect_message(
    prep_dt_depmap_feat(feat_data_path = test_feat_data_path,
                        meta_data_path = test_meta_data_path,
                        feature_set = "OmicsSignaturesProfile"),
    "The `OmicsSignaturesProfile` feature is not supported.")
  
  expect_error(prep_dt_depmap_feat(feat_data_path = 123, 
                                   meta_data_path = test_meta_data_path), 
               "Assertion on 'feat_data_path' failed: Must be of type 'string'")
  expect_error(prep_dt_depmap_feat(feat_data_path = test_feat_data_path,
                                   meta_data_path = "testdata/meta_data.qs"),
               "Assertion on 'File ext must be csv' failed: Must be TRUE")
  expect_error(prep_dt_depmap_feat(feat_data_path = test_feat_data_path,
                                   meta_data_path = "testdata/meta_data.csv"),
               "Assertion on 'meta_data_path' failed: File does not exist")
  expect_error(prep_dt_depmap_feat(feat_data_path = test_feat_data_path,
                                   meta_data_path = test_meta_data_path, 
                                   feature_set = 123),
               "Assertion on 'feature_set' failed: Must be of type 'string'")
  expect_error(prep_dt_depmap_feat(feat_data_path = test_feat_data_path,
                                   meta_data_path = test_meta_data_path,
                                   feature_set = "some_feat"),
               "Assertion on 'feat_path' failed: File does not exist")
  expect_error(prep_dt_depmap_feat(feat_data_path = test_feat_data_path,
                                   meta_data_path = test_meta_data_path,
                                   with_decoding = "decode"),
               "Assertion on 'with_decoding' failed: Must be of type 'logical flag'")
})


test_that("prep_dt_assoc works as expected", {
  # prep response data
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se_sa <- mae[[gDRutils::get_supported_experiments("sa")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se_sa,
                                                 assay_name = "Metrics")
  d_name <- "drug_002"
  d_name2 <- "drug_026"
  dt_response_met <- 
    prep_dt_response_metric_sa(dt_metrics, d_name,
                               metric = c("xc50", "x_mean", "x_max"))
  
  # prep depmap data
  obj_depmap_feat <- 
    prep_dt_depmap_feat(feat_data_path = system.file("testdata", package = "gDRplots"),
                        meta_data_path = system.file("testdata/Model.csv", package = "gDRplots"),
                        feature_set = "CRISPRGeneEffect")
  obj_depmap_meta <- 
    prep_dt_depmap_meta(meta_data_path = system.file("testdata/Model.csv", package = "gDRplots"))
  
  sel_met_1 <- "RV_gDR_x_mean"
  res_1 <- 
    prep_dt_assoc(dt_response = dt_response_met[, .SD, .SDcols = c("CellLineName", sel_met_1)],
                  dt_depmap = obj_depmap_feat[["dt_depmap"]])
  expect_is(res_1, "list")
  expect_length(res_1, 4)
  expect_equal(res_1[["selected_metric"]], sel_met_1)
  expect_null(res_1[["selected_feat_meta_col"]])
  expect_equal(names(res_1[["dt_assoc"]]), c("feature", "response", "rho", "q_value"))
  expect_equal(NROW(res_1[["dt_assoc"]]),
               NCOL(obj_depmap_feat[["dt_depmap"]][, .SD, .SDcols = -c("CCLEName", "ModelID")]))
  
  sel_met_2 <- "RV_gDR_log10_xc50"
  res_2 <- 
    prep_dt_assoc(dt_response = dt_response_met[, .SD, .SDcols = c("CellLineName", sel_met_2)][1:5, ],
                  dt_depmap = obj_depmap_feat[["dt_depmap"]],
                  selected_feat_meta_col = "XZ_fatures")
  expect_is(res_2, "list")
  expect_length(res_2, 4)
  expect_is(res_2[["dt_assoc"]], "data.table")
  ls_cl_res <- 
    dt_response_met[, .SD, .SDcols = c("CellLineName", sel_met_2)][is.finite(get(sel_met_2))][["CellLineName"]]
  ls_cl_dm <- stats::na.omit(obj_depmap_feat[["dt_depmap"]][CCLEName %in% ls_cl_res])
  expect_equal(NROW(res_2[["dt_assoc"]]), 0) # becaues ls_cl_dm < 6
  expect_null(res_2[["condition_info"]])
  expect_equal(res_2[["selected_feat_meta_col"]], "XZ_fatures")
  expect_equal(res_2[["selected_metric"]], sel_met_2)
  
  sel_met_3 <- "RV_gDR_x_max"
  res_3 <- 
    prep_dt_assoc(dt_response = dt_response_met[, .SD, .SDcols = c("CellLineName", sel_met_3)],
                  dt_depmap = obj_depmap_meta[["dt_depmap"]],
                  selected_feat_meta_col = "meta_xx")
  expect_is(res_3, "list")
  expect_length(res_3, 4)
  expect_equal(res_3[["selected_metric"]], sel_met_3)
  expect_equal(res_3[["selected_feat_meta_col"]], "meta_xx")
  expect_equal(names(res_3[["dt_assoc"]]), c("feature", "response", "rho", "q_value"))
  expect_equal(
    NROW(res_3[["dt_assoc"]]),
    sum(unlist(lapply(
      obj_depmap_meta[["dt_depmap"]][CCLEName %in% dt_response_met[["CellLineName"]], .SD, 
                                     .SDcols = -c("CCLEName", "ModelID")], 
      function(x) sum(x) > 0))))
  
  
  dt_res_not_num <- data.table::copy(dt_response_met[, .SD, .SDcols = c("CellLineName", sel_met_1)])
  dt_res_not_num[[sel_met_1]] <- as.character(dt_res_not_num[[sel_met_1]])
  expect_error(
    prep_dt_assoc(dt_response = dt_res_not_num,
                  dt_depmap = obj_depmap_feat[["dt_depmap"]]),
    "Column in `dt_response` with metric should be numeric.")
  expect_error(
    prep_dt_assoc(dt_response = dt_response_met,
                  dt_depmap = obj_depmap_feat[["dt_depmap"]]),
    "Provide `dt_response` for one metric.")
  dt_depmap_not_num <- data.table::copy(obj_depmap_feat[["dt_depmap"]])
  cols <- names(dt_depmap_not_num)
  dt_depmap_not_num[, (cols) := lapply(.SD, function(i) as.character(i)), .SDcols = cols]
  expect_error(
    prep_dt_assoc(dt_response = dt_response_met[, .SD, .SDcols = c("CellLineName", sel_met_2)],
                  dt_depmap = dt_depmap_not_num),
    "Provide `dt_depmap` with numeric values.")
})


test_that(".prep_dt_OmicsArmLevelCNA works as expected", {
  data_path <- system.file("testdata/OmicsArmLevelCNA.csv", package = "gDRplots")
  tab_raw <- data.table::fread(data_path)
  ls_chro <- grep("^[0-9].*(p$|q$)", names(tab_raw), value = TRUE)
  id_col <- setdiff(names(tab_raw), ls_chro)
  
  res_1 <- .prep_dt_OmicsArmLevelCNA(tab_raw)
  expect_is(res_1, "data.table")
  expect_equal(names(res_1),
               c(id_col, paste0(ls_chro, "_loss"), paste0(ls_chro, "_gain")))
  expect_true(all(
    vapply(ls_chro, function(i) { 
      all(res_1[[id_col]][res_1[[paste0(ls_chro[i], "_loss")]] == 1] == 
            tab_raw[[id_col]][tab_raw[[ls_chro[i]]] == -1])
    }, logical(1))))
  expect_true(all(
    vapply(ls_chro, function(i) { 
      all(res_1[[id_col]][res_1[[paste0(ls_chro[i], "_gain")]] == 1] == 
            tab_raw[[id_col]][tab_raw[[ls_chro[i]]] == 1])
    }, logical(1))))
  
  # scenario: check recording inside prep_dt_depmap_feat
  test_meta_data_path <- system.file("testdata/Model.csv", package = "gDRplots")
  test_feat_data_path <- system.file("testdata", package = "gDRplots")
  
  res_2 <- prep_dt_depmap_feat(feat_data_path = test_feat_data_path,
                               meta_data_path = test_meta_data_path,
                               feature_set = "OmicsArmLevelCNA",
                               with_decoding = TRUE)
  expect_is(res_2, "list")
  expect_is(res_2[[1]], "data.table")
  expect_equal(res_2$selected_feat_meta_col, "OmicsArmLevelCNA")
  expect_equal(res_2[[1]][, .SD, .SDcols = -c("ModelID", "CCLEName")],
               res_1[, .SD, .SDcols = -id_col])
  
  res_3 <- prep_dt_depmap_feat(feat_data_path = test_feat_data_path,
                               meta_data_path = test_meta_data_path,
                               feature_set = "OmicsArmLevelCNA",
                               with_decoding = FALSE)
  expect_is(res_3, "list")
  expect_is(res_3[[1]], "data.table")
  expect_equal(res_3[[1]][, .SD, .SDcols = -c("ModelID", "CCLEName")],
               tab_raw[, .SD, .SDcols = -id_col])
  
  expect_error(.prep_dt_OmicsArmLevelCNA(dt_depmap = as.list(tab_raw)),
               "Assertion on 'dt_depmap' failed: Must be a data.table")
})

