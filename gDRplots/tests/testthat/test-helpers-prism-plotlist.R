context("Test helpers-prism-plotlist")

test_that("create_PRISM_plot_list_sa works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix_small")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Averaged")
  d_names <- c("drug_004", "drug_021")
  feature_sets <- c("CRISPRGeneEffect", "OmicsCNGene")
  prefixes <- c("KO_", "CN_")
  metadata_columns <- "OncotreeLineage"
  
  # testing assertions
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = 1:3,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         feature_sets = feature_sets,
                                         prefixes = prefixes),
               "Assertion on 'drug_name_vec' failed: Must be of type 'character'")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = unlist(dt_average),
                                         feature_sets = feature_sets,
                                         prefixes = prefixes),
               "Assertion on 'dt_average' failed: Must be a data.table")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         feature_sets = 1:2,
                                         prefixes = prefixes),
               "Assertion on 'feature_sets' failed: Must be of type 'character'")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         feature_sets = feature_sets[1],
                                         prefixes = prefixes),
               "`prefixes` has to be the same length as `feature_sets`")
  expect_error(create_PRISM_plot_list_sa(drug_name_vec = d_names,
                                         dt_metrics = dt_metrics,
                                         dt_average = dt_average,
                                         feature_sets = feature_sets,
                                         prefixes = 1:2),
               "Assertion on 'prefixes' failed: Must be of type 'character'")
})