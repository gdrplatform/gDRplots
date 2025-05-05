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
  metadata_columns <- "OncotreeLineage"
  feat_data_path <- file.path(tempdir(), "depmap_data")
  dir.create(feat_data_path, showWarnings = FALSE)
  feat_file <- file.path(feat_data_path, "CRISPRGeneEffect.csv")
  writeLines(text = "CRISPRGeneEffect data", con = feat_file)
  file.exists(feat_file)
  dir.exists(feat_data_path)
  
  on.exit({
    unlink(feat_data_path, recursive = TRUE)
    unlink(feat_file)
  })
  feature_sets <- c("CRISPRGeneEffect", "OmicsCNGene")
  
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
               "Provide consistent values for `feature_sets` and `feat_data_path` for DepMam subset.")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         meta_data_path = meta_data_path,
                                         feat_data_path = NULL,
                                         feature_sets = NULL,
                                         metadata_columns = NULL),
               "Provide `feature_sets` or `metadata_columns` for DepMam subset.")
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
               "Provide consistent values for `feature_sets` and `feat_data_path` for DepMam subset.")
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
  metadata_columns <- "OncotreeLineage"
  
  meta_data_path <- system.file("testdata/Model.csv", package = "gDRplots")
  metadata_columns <- "OncotreeLineage"
  feat_data_path <- file.path(tempdir(), "depmap_data")
  dir.create(feat_data_path, showWarnings = FALSE)
  feat_file <- file.path(feat_data_path, "CRISPRGeneEffect.csv")
  writeLines(text = "CRISPRGeneEffect data", con = feat_file)
  file.exists(feat_file)
  dir.exists(feat_data_path)
  
  on.exit({
    unlink(feat_data_path, recursive = TRUE)
    unlink(feat_file)
  })
  feature_sets <- c("CRISPRGeneEffect", "OmicsCNGene")
  
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
               "Provide consistent values for `feature_sets` and `feat_data_path` for DepMam subset.")
  expect_error(create_PRISM_plot_list_combo(drug1_name_vec = d_names,
                                            drug2_name_vec = d_names_2,
                                            dt_metrics = dt_metrics,
                                            dt_scores = dt_scores,
                                            meta_data_path = meta_data_path,
                                            feat_data_path = NULL,
                                            feature_sets = NULL,
                                            metadata_columns = NULL),
               "Provide `feature_sets` or `metadata_columns` for DepMam subset.")
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
               "Provide consistent values for `feature_sets` and `feat_data_path` for DepMam subset.")
})
