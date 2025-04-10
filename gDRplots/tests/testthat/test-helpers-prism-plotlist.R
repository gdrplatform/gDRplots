context("Test helpers-prism-plotlist")

test_that("create_PRISM_plot_list_sa works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Averaged")
  d_names <- c("drug_004", "drug_021")
  
  meta_data_path <- system.file("testdata/Model.csv", package = "gDRplots")
  feature_sets <- c("CRISPRGeneEffect", "OmicsCNGene")
  prefixes <- c("KO_", "CN_")
  metadata_columns <- "OncotreeLineage"
  
  # TODO in GDR-2710
  
  # testing assertions
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = 1:3,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = meta_data_path,
                                         feature_sets = feature_sets,
                                         prefixes = prefixes),
               "Assertion on 'drug_name_vec' failed: Must be of type 'character'")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = unlist(dt_metrics),
                                         dt_average = dt_average,
                                         meta_data_path = meta_data_path,
                                         feature_sets = feature_sets,
                                         prefixes = prefixes),
               "Assertion on 'dt_metrics' failed: Must be a data.table")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = unlist(dt_average),
                                         meta_data_path = meta_data_path,
                                         feature_sets = feature_sets,
                                         prefixes = prefixes),
               "Assertion on 'dt_average' failed: Must be a data.table")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = 1:3,
                                         feature_sets = feature_sets,
                                         prefixes = prefixes),
               "Assertion on 'meta_data_path' failed: Must be of type 'string'")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,   
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = "path/meta_data.qs",
                                         feature_sets = feature_sets,
                                         prefixes = prefixes),
               "Assertion on 'File ext must be csv' failed: Must be TRUE")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = meta_data_path,
                                         feature_sets = 1:2,
                                         prefixes = prefixes),
               "Assertion on 'feature_sets' failed: Must be of type 'character'")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = meta_data_path,
                                         feature_sets = feature_sets[1],
                                         prefixes = prefixes),
               "`prefixes` has to be the same length as `feature_sets`")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = meta_data_path,
                                         feature_sets = feature_sets,
                                         prefixes = 1:2),
               "Assertion on 'prefixes' failed: Must be of type 'character'")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = meta_data_path,
                                         feature_sets = NULL,
                                         prefixes = NULL),
               "Provide `feature_sets` or `metadata_columns` for DepMam subset")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = meta_data_path,
                                         feature_sets = NULL,
                                         prefixes = NULL,
                                         metadata_columns = 1:3),
               "Assertion on 'metadata_columns' failed: Must be of type 'character'")
})

test_that("create_PRISM_plot_list_combo works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Metrics")
  dt_scores <- gDRutils::convert_se_assay_to_dt(se = se,
                                                assay_name = "scores")
  d_names <- c("drug_001", "drug_011")
  d_names_2 <- c("drug_026", "drug_031")
  
  meta_data_path <- system.file("testdata/Model.csv", package = "gDRplots")
  feature_sets <- c("CRISPRGeneEffect", "OmicsCNGene")
  prefixes <- c("KO_", "CN_")
  metadata_columns <- "OncotreeLineage"
  
  # TODO in GDR-2710
  
  # testing assertions
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = 1:3,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feature_sets = feature_sets,
                                            prefixes = prefixes),
               "Assertion on 'drug1_name_vec' failed: Must be of type 'character'")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = c(NA, NA),
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feature_sets = feature_sets,
                                            prefixes = prefixes),
               "Assertion on 'drug2_name_vec' failed: Contains only missing values")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = data.table::data.table(),
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feature_sets = feature_sets,
                                            prefixes = prefixes),
               "Assertion on 'dt_metrics' failed: Must have at least 1 rows")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = unlist(dt_scores),
                                            meta_data_path = meta_data_path,
                                            feature_sets = feature_sets,
                                            prefixes = prefixes),
               "Assertion on 'dt_scores' failed: Must be a data.table")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = 1:3,
                                            feature_sets = feature_sets,
                                            prefixes = prefixes),
               "Assertion on 'meta_data_path' failed: Must be of type 'string'")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = "path/meta_data.qs",
                                            feature_sets = feature_sets,
                                            prefixes = prefixes),
               "Assertion on 'File ext must be csv' failed: Must be TRUE")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feature_sets = 1:2,
                                            prefixes = prefixes),
               "Assertion on 'feature_sets' failed: Must be of type 'character'")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feature_sets = feature_sets[1],
                                            prefixes = prefixes),
               "`prefixes` has to be the same length as `feature_sets`")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feature_sets = feature_sets,
                                            prefixes = 1:2),
               "Assertion on 'prefixes' failed: Must be of type 'character'")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feature_sets = NULL,
                                            prefixes = NULL),
               "Provide `feature_sets` or `metadata_columns` for DepMam subset")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feature_sets = NULL,
                                            prefixes = NULL,
                                            metadata_columns = 1:3),
               "Assertion on 'metadata_columns' failed: Must be of type 'character'")
})