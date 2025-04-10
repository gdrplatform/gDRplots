context("Test helpers-associations")

cellline_name <- gDRutils::get_env_identifiers("cellline_name")

# data ----
mae <- gDRutils::get_synthetic_data("combo_matrix")
se_sa <- mae[[gDRutils::get_supported_experiments("sa")]]
dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se_sa,
                                               assay_name = "Metrics")
d_name <- "drug_002"
dt_response_met <- 
  prep_dt_response_metric_sa(dt_metrics, d_name,
                             metric = c("xc50", "x_mean", "x_max"))

# fake depmap data
cell_lines <- gDRtestData::create_synthetic_cell_lines()[["CellLineName"]]
dt_model <- data.table::data.table(
  ModelID = sprintf("ACH-%06d", seq_along(cell_lines)),
  CCLEName = cell_lines
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

# inputs
selected_metric <- "RV_gDR_log10_xc50"
Y <- as.matrix(
  dt_response_met[, .SD, .SDcols = c(cellline_name, selected_metric)],
  rownames = "CellLineName"
)
Y_all <- as.matrix(
  dt_response_met[, .SD, .SDcols = -c("rId", "cId")],
  rownames = "CellLineName"
)
Y_vec <- dt_response_met$RV_gDR_x_mean
names(Y_vec) <- dt_response_met[[cellline_name]]
X <- as.matrix(
  dt_depmap_meta[CCLEName %in% dt_response_met[[cellline_name]], .SD, .SDcols = -c("ModelID")], 
  rownames = "CCLEName"
)
res_col_names <- c("feature", "est_beta", "est_beta_se", 
                   "posterior_mean", "posterior_sd", "prob_negative", "prob_positive", 
                   "rho", "p_value", "q_value", "s_value", "lfsr", "lfdr")

# tests ----
test_that("calc_assoc works as expected", {
  
  res_1 <- calc_assoc(X, Y) # default matrix
  expect_is(res_1, "data.table")
  expect_true(all(c(res_col_names, "response") %in% names(res_1)))
  
  res_2 <- calc_assoc(X, Y_vec) # default vector
  expect_is(res_2, "data.table")
  expect_true(all(res_2$feature %in% colnames(X)[colSums(X[names(Y_vec), ]) > 0]))
  
  # scenario: matrix Y has more than one column
  res_3 <- calc_assoc(X, Y_all)
  expect_is(res_3, "data.table")
  

  expect_error(calc_assoc(X = matrix(LETTERS[1:9], nrow = 3), Y),
               "Assertion on 'X' failed: Must store numerics")
  expect_error(calc_assoc(X = matrix(1:9, nrow = 3), Y),
               "Must have names")
  expect_error(calc_assoc(X = X, Y = LETTERS[1:9]),
               "Assertion on 'Y' failed: Must inherit from class")
  expect_error(calc_assoc(X = X, Y = dt_response_met$RV_gDR_x_mean),
               "Must have names")
})

test_that(".calc_assoc_vector works as expected", {
  
  # scenario: vector Y has no variance 
  Y_vec_no_var <- rep(8, NROW(X))
  names(Y_vec_no_var) <- names(Y_vec)
  expect_warning({
    res_4 <- calc_assoc(X, Y = Y_vec_no_var)
  }, "Y has no variance")
  expect_is(res_4, "data.table")
  expect_true(all(res_col_names %in% names(res_4)))
  expect_true(all(res_4[, lapply(.SD, is.na), .SDcols = -c("feature")]))
  expect_true(all(res_4[, lapply(.SD, is.numeric), .SDcols = -c("feature")]))
})

test_that(".calc_assoc_matrix works as expected", {
  
  # scenario: matrix Y has no variance 
  Y_mat_no_var <- matrix(rep(0.88, NROW(X)), ncol = 1,
                         dimnames = list(rownames(Y), selected_metric))
  expect_warning({
    res_5 <- calc_assoc(X, Y = Y_mat_no_var) 
  }, "The following columns in Y have no variance")
  expect_is(res_5, "data.table")
  expect_true(all(c(res_col_names, "response") %in% names(res_5)))
  expect_true(all(res_5[, lapply(.SD, is.na), .SDcols = -c("feature", "response")]))
  expect_true(all(res_5[, lapply(.SD, is.numeric), .SDcols = -c("feature", "response")]))
})
