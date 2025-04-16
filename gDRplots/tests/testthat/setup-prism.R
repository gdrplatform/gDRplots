cellline_name <- gDRutils::get_env_identifiers("cellline_name")

# gDR data ----
mae <- gDRutils::get_synthetic_data("combo_matrix")
se_sa <- mae[[gDRutils::get_supported_experiments("sa")]]
dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se_sa,
                                               assay_name = "Metrics")
dt_average <- gDRutils::convert_se_assay_to_dt(se = se_sa,
                                               assay_name = "Averaged")
dt_metrics_capped <-
  gDRutils::cap_assay_infinities(conc_assay_dt = dt_average,
                                 assay_dt = dt_metrics,
                                 experiment_name = gDRutils::get_supported_experiments("sa"),
                                 capping_fold = 5)

se_combo <- mae[[gDRutils::get_supported_experiments("combo")]]
dt_metrics_combo <- gDRutils::convert_se_assay_to_dt(se = se_combo,
                                                     assay_name = "Metrics")
dt_scores <- gDRutils::convert_se_assay_to_dt(se = se_combo,
                                              assay_name = "scores")

d_name <- "drug_002"
d_name2 <- "drug_026"
dt_response_met <- 
  prep_dt_response_metric_sa(dt_metrics, d_name,
                             metric = c("xc50", "x_mean", "x_max"))
dt_response_met_capped <- 
  prep_dt_response_metric_sa(dt_metrics_capped, d_name,
                             metric = c("xc50", "x_mean", "x_max"))
dt_response_dose <- 
  prep_dt_response_dose_sa(dt_average, d_name)

dt_response_score <- 
  prep_dt_response_scores(dt_scores, d_name, d_name2, 
                          metric = c("hsa_score", "bliss_score"))
dt_response_diff <-
  prep_dt_response_metric_diff(dt_metrics_combo, d_name, d_name2,
                               metric = c("xc50", "x_mean", "x_max"))


# fake depmap data ----
cell_lines <- gDRtestData::create_synthetic_cell_lines()[["CellLineName"]]
dt_model <- data.table::data.table(
  ModelID = sprintf("ACH-%06d", seq_along(cell_lines)),
  CCLEName = cell_lines
)
drugs <- gDRtestData::create_synthetic_drugs()[["DrugName"]]

#_feature_sets
dt_depmap_feat <- data.table::data.table(
  CCLEName = cell_lines,
  "XZ_A1QW" = withr::with_seed(42, rnorm(n = NROW(cell_lines), mean = -0.05, sd = 0.11)),
  "XZ_A2GH" = withr::with_seed(42, rnorm(n = NROW(cell_lines), mean = -0.03, sd = 0.13)),
  "XZ_A3OP" = withr::with_seed(42, rnorm(n = NROW(cell_lines), mean = 0.045, sd = 0.10)),
  "XZ_A4RT" = withr::with_seed(42, rnorm(n = NROW(cell_lines), mean = 0.05, sd = 0.10)),
  "XZ_A5BN" = withr::with_seed(42, rnorm(n = NROW(cell_lines), mean = 0.11, sd = 0.13))
)
dt_depmap_feat[CCLEName %in% c("cellline_FD", "cellline_NE"), 
               setdiff(names(dt_depmap_feat), "CCLEName") := NA]
dt_depmap_feat <- merge(dt_model, dt_depmap_feat, by = "CCLEName")
data.table::setkey(dt_depmap_feat, NULL)

obj_depmap_feat <- list(
  dt_depmap = dt_depmap_feat,
  selected_feat_meta_col = "XZ_fatures"
)

dt_depmap_feat_2 <- data.table::data.table(
  CCLEName = cell_lines,
  "NU_X1QW" = withr::with_seed(42, sample(c(0, 1), size = NROW(cell_lines), replace = TRUE)),
  "NU_X2GH" = withr::with_seed(314, sample(c(0, 1), size = NROW(cell_lines), replace = TRUE)),
  "NU_X3OP" = withr::with_seed(271, sample(c(0, 1), size = NROW(cell_lines), replace = TRUE)),
  "NU_X4RT" = withr::with_seed(981, sample(c(0, 1), size = NROW(cell_lines), replace = TRUE)),
  "NU_X5BN" = rep(1, size = NROW(cell_lines))
)
dt_depmap_feat_2[CCLEName %in% c("cellline_FD", "cellline_NE"), 
                 setdiff(names(dt_depmap_feat_2), "CCLEName") := NA]
dt_depmap_feat_2 <- merge(dt_model, dt_depmap_feat_2, by = "CCLEName")
data.table::setkey(dt_depmap_feat_2, NULL)

obj_depmap_feat_2 <- list(
  dt_depmap = dt_depmap_feat_2,
  selected_feat_meta_col = "NU_fatures"
)

#_meta
dt_depmap_meta_lng <- data.table::data.table(
  CCLEName = cell_lines,
  meta_xx = withr::with_seed(42, sample(x = sprintf("meta_%s", c("AA", "BB", "CC")), 
                                        size = NROW(cell_lines), replace = TRUE))
)
dt_depmap_meta_lng[CCLEName %in% c("cellline_OO", "cellline_AA"), ][["meta_xx"]] <- NA
dt_depmap_meta_lng[CCLEName == "cellline_FD", ][["meta_xx"]] <- "meta_DD"
dt_depmap_meta_lng[CCLEName == "cellline_NE", ][["meta_xx"]] <- "longer_than_other_meta_EE"

dt_depmap_meta <- data.table::dcast(data = dt_depmap_meta_lng, 
                                    formula = CCLEName ~ meta_xx, 
                                    fun.aggregate = length)
data.table::setkey(dt_depmap_meta, NULL)
dt_depmap_meta <- merge(dt_model, dt_depmap_meta, by = "CCLEName")

obj_depmap_meta <- list(
  dt_depmap = dt_depmap_meta,
  selected_feat_meta_col = "meta_xx"
)

#_assoc
ls_feat <- setdiff(names(obj_depmap_feat[["dt_depmap"]]), c("ModelID", "CCLEName"))
dist_q <- withr::with_seed(42, rnorm(400, mean = 0.3, sd = 0.05))

dt_assoc_sa <- data.table::data.table(
  feature = ls_feat,
  response = rep("RV_gDR_xc50", NROW(ls_feat)),
  rho = withr::with_seed(42, rnorm(NROW(ls_feat), mean = 0, sd = 0.05)),
  q_value = withr::with_seed(42, sample(dist_q[dist_q < 0.3], NROW(ls_feat), replace = TRUE))
)

obj_assoc_sa <- list(dt_assoc = dt_assoc_sa,
                     condition_info = unique(dt_response_met[["rId"]]),
                     selected_metric = "RV_gDR_xc50",
                     selected_feat_meta_col = "XZ_fatures")

ls_meta <- setdiff(names(obj_depmap_meta[["dt_depmap"]]), c("ModelID", "CCLEName"))
dist_q <- withr::with_seed(42, rnorm(400, mean = 1, sd = 0.1))

dt_assoc_combo <- data.table::data.table(
  feature = ls_meta,
  response = rep("hsa_score", NROW(ls_meta)),
  rho = withr::with_seed(42, rnorm(NROW(ls_meta), mean = 0, sd = 0.05)),
  q_value = withr::with_seed(42, sample(dist_q[dist_q < 0.75], NROW(ls_meta), replace = TRUE))
)

obj_assoc_combo <- list(dt_assoc = dt_assoc_combo,
                        condition_info = unique(dt_response_score[["rId"]]),
                        selected_metric = "hsa_score",
                        selected_feat_meta_col = "meta_xx")
