context("Test helpers-prism-plotlist")

test_that("create_PRISM_plot_list_sa works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Averaged")
  d_names <- c("drug_021", "drug_026")
  
  meta_data_path <- system.file("testdata/Model.csv", package = "gDRplots")
  metadata_columns <- c("OncotreeLineage", "Sex", "PatientRace")
  feat_data_path <- system.file("testdata", package = "gDRplots")
  feature_sets <- c("CRISPRGeneEffect", "OmicsSomaticMutationsMatrixHotspot")
  
  res_1 <- create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                     dt_metrics = dt_metrics,
                                     dt_average = dt_average,
                                     meta_data_path = meta_data_path,
                                     feat_data_path = feat_data_path,
                                     feature_sets = feature_sets) # default
  expect_is(res_1, "list")
  expect_length(res_1, NROW(feature_sets))
  expect_length(res_1[[1]], NROW(d_names))
  expect_length(res_1[[1]][[1]], 1) # only RV
  expect_length(res_1[[1]][[1]][[1]],
                NROW(c("xc50", "x_mean", "x_max")) + NROW(unique(dt_average$Concentration))) # metrics
  
  # scenario: only dt_metrics and not available feature
  res_2 <- create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                     dt_metrics = dt_metrics,
                                     dt_average = NULL,
                                     metric = "x_mean",
                                     normalization_type_vec = c("RV", "GR"),
                                     meta_data_path = meta_data_path,
                                     feat_data_path = feat_data_path,
                                     feature_sets = c(feature_sets, "non_available_feature"))
  expect_is(res_2, "list")
  expect_length(res_2, NROW(feature_sets))
  expect_equal(names(res_2), feature_sets)
  expect_length(res_2[[1]], NROW(c("RV", "GR"))) 
  expect_length(res_2[[1]], NROW(d_names))
  expect_length(res_2[[1]][[1]][["RV"]], NROW(c("x_mean")))
  expect_is(res_2[[1]][[1]][["RV"]][[1]], "ggplot")
  
  # scenario: all meta columns and features
  res_3 <- create_PRISM_plot_list_sa(drug_name_vec = c(d_names, "non_available_drug"),
                                     dt_metrics = dt_metrics,
                                     dt_average = dt_average,
                                     normalization_type_vec = "GR",
                                     meta_data_path = meta_data_path,
                                     feat_data_path = feat_data_path,
                                     feature_sets = feature_sets[1],
                                     metadata_columns = metadata_columns)
  expect_is(res_3, "list")
  expect_length(res_3, NROW(feature_sets[1]) + NROW(metadata_columns))
  expect_equal(names(res_3), c(feature_sets[1], metadata_columns))
  expect_length(res_3[[1]], NROW(d_names))
  expect_equal(names(res_3[[1]]), d_names)
  expect_length(res_3[[1]][[1]], 1) # only GR
  expect_length(res_3[[1]][[1]][[1]],
                NROW(c("xc50", "x_mean", "x_max")) + NROW(unique(dt_average$Concentration))) # metrics
  
  # scenario: only dt_average and not available meta
  res_4 <- create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                     dt_metrics = NULL,
                                     dt_average = dt_average,
                                     normalization_type_vec = c("RV", "GR"),
                                     meta_data_path = meta_data_path,
                                     feat_data_path = NULL,
                                     feature_sets = feature_sets,
                                     metadata_columns = c(metadata_columns[1], "non_available_meta"))
  expect_is(res_4, "list")
  expect_length(res_4, NROW(metadata_columns[1]))
  expect_equal(names(res_4), metadata_columns[1])
  expect_length(res_4[[1]], NROW(d_names))
  expect_length(res_4[[1]][[1]], NROW(c("RV", "GR")))
  expect_length(res_4[[1]][[1]][[1]],
                NROW(unique(dt_average$Concentration))) # metrics
  
  # scenario: no data for given drugs
  dt_metrics_GR <- data.table::copy(dt_metrics)[normalization_type == "GR", ]
  expect_message({
    res_5 <- create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                       dt_metrics = dt_metrics_GR,
                                       dt_average = dt_average,
                                       meta_data_path = meta_data_path,
                                       feat_data_path = feat_data_path,
                                       feature_sets = feature_sets)
  }, "There was no data for selected drugs.")
  expect_length(res_5, 0)
  
  # testing assertions
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = 1:3,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = meta_data_path,
                                         feat_data_path = feat_data_path,
                                         feature_sets = feature_sets[1]),
               "Assertion on 'drug_name_vec' failed: Must be of type 'character'")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = unlist(dt_metrics),
                                         dt_average = dt_average,
                                         meta_data_path = meta_data_path,
                                         feat_data_path = feat_data_path,
                                         feature_sets = feature_sets),
               "Assertion on 'dt_metrics' failed: Must be a data.table")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = unlist(dt_average),
                                         meta_data_path = meta_data_path,
                                         feat_data_path = feat_data_path,
                                         feature_sets = feature_sets[1]),
               "Assertion on 'dt_average' failed: Must be a data.table")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = NULL,
                                         dt_average = NULL,
                                         meta_data_path = meta_data_path,
                                         feat_data_path = feat_data_path,
                                         feature_sets = feature_sets[1]),
               "Provide response data - at least one of `dt_metrics` or `dt_average`.")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         metric = NULL,
                                         meta_data_path = meta_data_path,
                                         feat_data_path = feat_data_path,
                                         feature_sets = feature_sets[1]),
               "Assertion on 'metric' failed: Must be a subset of")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         metric = "not_known_metric",
                                         meta_data_path = meta_data_path,
                                         feat_data_path = feat_data_path,
                                         feature_sets = feature_sets[1]),
               "Assertion on 'metric' failed: Must be a subset of")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = 1:3,
                                         feat_data_path = feat_data_path,
                                         feature_sets = feature_sets[1]),
               "Assertion on 'meta_data_path' failed: Must be of type 'string'")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,   
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = "path/meta_data.qs",
                                         feat_data_path = feat_data_path,
                                         feature_sets = feature_sets),
               "Assertion on 'File ext must be csv' failed: Must be TRUE")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,   
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = meta_data_path,
                                         feat_data_path = 1,
                                         feature_sets = feature_sets[1]),
               "Assertion on 'feat_data_path' failed: Must be of type 'string'")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,   
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = meta_data_path,
                                         feat_data_path = "path/that/not_exist",
                                         feature_sets = feature_sets[1]),
               "Assertion on 'feat_data_path' failed: Directory 'path/that/not_exist' does not exist")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = meta_data_path,
                                         feat_data_path = feat_data_path,
                                         feature_sets = 1:2),
               "Assertion on 'feature_sets' failed: Must be of type 'character'")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = meta_data_path,
                                         feat_data_path = feat_data_path,
                                         feature_sets = NULL),
               "Provide consistent values for `feature_sets` and `feat_data_path` for DepMap subset.")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = meta_data_path,
                                         feat_data_path = NULL,
                                         feature_sets = NULL,
                                         metadata_columns = NULL),
               "Provide `feature_sets` or `metadata_columns` for DepMap subset.")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = meta_data_path,
                                         feat_data_path = feat_data_path,
                                         feature_sets = NULL,
                                         metadata_columns = 1:3),
               "Assertion on 'metadata_columns' failed: Must be of type 'character'")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = meta_data_path,
                                         feat_data_path = NULL,
                                         feature_sets = feature_sets,
                                         metadata_columns = NULL),
               "Provide consistent values for `feature_sets` and `feat_data_path` for DepMap subset.")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = meta_data_path,
                                         feat_data_path = feat_data_path,
                                         feature_sets = feature_sets,
                                         clear_taxonomy_info = "str"),
               "Assertion on 'clear_taxonomy_info' failed: Must be of type 'logical flag'")

})

test_that("create_PRISM_plot_list_combo works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix_small")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Metrics")
  dt_scores <- gDRutils::convert_se_assay_to_dt(se = se,
                                                assay_name = "scores")
  d_names <- c("drug_004", "drug_005")
  d_names_2 <- "drug_021"
  
  meta_data_path <- system.file("testdata/Model.csv", package = "gDRplots")
  metadata_columns <- c("OncotreeLineage", "Sex", "PatientRace")
  feat_data_path <- system.file("testdata", package = "gDRplots")
  feature_sets <- c("CRISPRGeneEffect", "OmicsSomaticMutationsMatrixHotspot")
  
  res_1 <- create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                        drug2_name_vec = d_names_2,
                                        dt_metrics = dt_metrics,
                                        dt_scores = dt_scores,
                                        meta_data_path = meta_data_path,
                                        feat_data_path = feat_data_path,
                                        feature_sets = feature_sets) # default
  expect_is(res_1, "list")
  expect_length(res_1, NROW(feature_sets))
  expect_length(res_1[[1]], NROW(expand.grid(d_names, d_names_2, stringsAsFactors = FALSE)))
  expect_length(res_1[[1]][[1]], 1) # only RV
  expect_true(all(vapply(c("xc50", "x_mean", "x_max", "hsa_score", "bliss_score"), function(met) {
    any(grepl(met, names(res_1[[1]][[1]][[1]])))
  }, FUN.VALUE = logical(1)))) # metrics
  
  # scenario: only dt_metrics and not available feature
  res_2 <- create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                        drug2_name_vec = d_names_2,
                                        dt_metrics = dt_metrics,
                                        dt_scores = NULL,
                                        normalization_type_vec = c("RV", "GR"),
                                        meta_data_path = meta_data_path,
                                        feat_data_path = feat_data_path,
                                        feature_sets = c(feature_sets[2], "non_available_feature"))
  expect_is(res_2, "list")
  expect_length(res_2, NROW(feature_sets[2]))
  expect_equal(names(res_2), feature_sets[2])
  expect_length(res_2[[1]], NROW(c("RV", "GR"))) 
  expect_length(res_2[[1]], NROW(expand.grid(d_names, d_names_2, stringsAsFactors = FALSE)))
  expect_false(all(vapply(c("hsa_score", "bliss_score"), function(met) {
    any(grepl(met, names(res_2[[1]][[1]][[1]])))
  }, FUN.VALUE = logical(1)))) # dt_scores = NULL nolint
  expect_is(res_2[[1]][[1]][["RV"]][[1]], "ggplot")
  
  # scenario: all meta columns and features
  res_3 <- create_PRISM_plot_list_combo(drug1_name_vec = c(d_names, "non_available_drug"),
                                        drug2_name_vec = d_names_2,
                                        dt_metrics = dt_metrics,
                                        dt_scores = NULL,
                                        normalization_type_vec = "GR",
                                        metric = "x_max",
                                        meta_data_path = meta_data_path,
                                        feat_data_path = feat_data_path,
                                        feature_sets = feature_sets[1],
                                        metadata_columns = metadata_columns)
  expect_is(res_3, "list")
  expect_length(res_3, NROW(feature_sets[1]) + NROW(metadata_columns))
  expect_equal(names(res_3), c(feature_sets[1], metadata_columns))
  expect_length(res_3[[1]],  NROW(expand.grid(d_names, d_names_2, stringsAsFactors = FALSE)))
  expect_length(res_3[[1]][[1]], 1) # only GR
  expect_true(all(grepl("GR_gDR_x_max", names(res_3[[1]][[1]][[1]])))) # metrics
  
  # scenario: only dt_average and not available meta
  res_4 <- create_PRISM_plot_list_combo(drug1_name_vec = d_names[1],
                                        drug2_name_vec = d_names_2,
                                        dt_metrics = NULL,
                                        dt_scores = dt_scores,
                                        metric_scores = "hsa_score",
                                        normalization_type_vec = "GR",
                                        meta_data_path = meta_data_path,
                                        feat_data_path = NULL,
                                        feature_sets = feature_sets,
                                        metadata_columns = c(metadata_columns[1], "non_available_meta"))
  expect_is(res_4, "list")
  expect_length(res_4, NROW(metadata_columns[1]))
  expect_equal(names(res_4), metadata_columns[1])
  expect_equal(names(res_4[[1]]), sprintf("%s x %s", d_names[1], d_names_2))
  expect_length(res_4[[1]][[1]], NROW(c("GR")))
  expect_true(all(grepl("GR_gDR_hsa_score", names(res_4[[1]][[1]][[1]])))) # metrics
  
  # scenario: no data for given drugs
  dt_metrics_GR <- data.table::copy(dt_metrics)[normalization_type == "GR", ]
  expect_message({
    res_5 <- create_PRISM_plot_list_combo(drug1_name_vec = d_names[1],
                                          drug2_name_vec = d_names_2,
                                          dt_metrics = dt_metrics_GR,
                                          dt_scores = dt_scores,
                                          normalization_type_vec = "RV",
                                          metric = "xc50",
                                          metric_scores = "hsa_score",
                                          meta_data_path = meta_data_path,
                                          feat_data_path = feat_data_path,
                                          feature_sets = feature_sets,
                                          metadata_columns = "non_available_meta")
  }, "There was no data for selected drugs combination.")
  expect_length(res_5, 0)
  
  # testing assertions
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = 1:3,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feat_data_path = feat_data_path,
                                            feature_sets = feature_sets),
               "Assertion on 'drug1_name_vec' failed: Must be of type 'character'")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = c(NA, NA),
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feat_data_path = feat_data_path,
                                            feature_sets = feature_sets),
               "Assertion on 'drug2_name_vec' failed: Contains only missing values")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = data.table::data.table(),
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feat_data_path = feat_data_path,
                                            feature_sets = feature_sets),
               "Assertion on 'dt_metrics' failed: Must have at least 1 rows")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = unlist(dt_scores),
                                            meta_data_path = meta_data_path,
                                            feat_data_path = feat_data_path,
                                            feature_sets = feature_sets),
               "Assertion on 'dt_scores' failed: Must be a data.table")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = NULL,
                                            dt_scores = NULL,
                                            meta_data_path = meta_data_path,
                                            feat_data_path = feat_data_path,
                                            feature_sets = feature_sets),
               "Provide response data - at least one of `dt_metrics` or `dt_scores`.")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = NULL,
                                            metric = NULL,
                                            meta_data_path = meta_data_path,
                                            feat_data_path = feat_data_path,
                                            feature_sets = feature_sets),
               "Assertion on 'metric' failed: Must be a subset of")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = 1:3,
                                            feature_sets = feature_sets),
               "Assertion on 'meta_data_path' failed: Must be of type 'string'")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = "path/meta_data.qs",
                                            feat_data_path = feat_data_path,
                                            feature_sets = feature_sets),
               "Assertion on 'File ext must be csv' failed: Must be TRUE")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feat_data_path = 1,
                                            feature_sets = feature_sets[1]),
               "Assertion on 'feat_data_path' failed: Must be of type 'string'")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feat_data_path = "path/that/not_exist",
                                            feature_sets = feature_sets[1]),
               "Assertion on 'feat_data_path' failed: Directory 'path/that/not_exist' does not exist")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feat_data_path = feat_data_path,
                                            feature_sets = 1:2),
               "Assertion on 'feature_sets' failed: Must be of type 'character'")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feat_data_path = feat_data_path,
                                            feature_sets = NULL),
               "Provide consistent values for `feature_sets` and `feat_data_path` for DepMap subset.")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feat_data_path = NULL,
                                            feature_sets = NULL,
                                            metadata_columns = NULL),
               "Provide `feature_sets` or `metadata_columns` for DepMap subset.")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feat_data_path = feat_data_path,
                                            feature_sets = NULL,
                                            metadata_columns = 1:3),
               "Assertion on 'metadata_columns' failed: Must be of type 'character'")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feat_data_path = NULL,
                                            feature_sets = feature_sets,
                                            metadata_columns = NULL),
               "Provide consistent values for `feature_sets` and `feat_data_path` for DepMap subset.")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feat_data_path = feat_data_path,
                                            feature_sets = feature_sets,
                                            clear_taxonomy_info = NULL),
               "Assertion on 'clear_taxonomy_info' failed: Must be of type 'logical flag'")
})
