context("Test helpers-associations")

# data ----
X <- matrix(c(rep(2, 5), rep(1, 6), rep(0, 2), rep(1, 2), rep(0, 3), rep(2, 8), rep(3, 6)), nrow = 8,
            dimnames = list(sprintf("row_%s", 1:8), sprintf("feat_%s", 1:4)))

Y <- matrix(c(10:17), ncol = 1,
            dimnames = list(sprintf("row_%s", 1:8), "met"))
Y_all <- matrix(c(10:17, 110:117, 210:217), ncol = 3,
                dimnames = list(sprintf("row_%s", 1:8), sprintf("met_%s", 1:3)))
Y_vec <- 1:8
names(Y_vec) <- sprintf("row_%s", 1:8)

# small matrix
X_feat <- matrix(c(rep(NA, 2), rep(1, 10), rep(NA, 2), rep(1, 2), rep(NA, 3), rep(2, 5), rep(1, 2), rep(3, 2)), 
                 nrow = 7,
                 dimnames = list(sprintf("row_%s", 1:7), sprintf("feat_%s", 1:4)))
Y_feat <- matrix(c(10:16), ncol = 1,
                 dimnames = list(sprintf("row_%s", 1:7), "met"))

# big matrix
n_big <- 50
X_big <- matrix(
  withr::with_seed(42, sample(c(0, 1, NA), size =  20 * n_big, replace = TRUE, prob = c(0.9, 0.09, 0.01))), 
  nrow = n_big, dimnames = list(sprintf("row_%s", 1:n_big), sprintf("feat_%s", 1:20)))

Y_big <- as.matrix(data.table::data.table(
  met_11 = withr::with_seed(42, rnorm(n = n_big, mean = -0.05, sd = 0.11)),
  met_12 = withr::with_seed(42, rnorm(n = n_big, mean = -0.03, sd = 0.13)),
  met_13 = withr::with_seed(42, sample(c(1:10, NA), n_big, replace = TRUE)),
  met_14 = rep(c(0.15, NA, NA), length.out = n_big)
),)
rownames(Y_big) <- sprintf("row_%02d", 1:n_big)

# association columns ----
res_col_names <- c("feature", "est_beta", "est_beta_se", 
                   "posterior_mean", "posterior_sd", "prob_negative", "prob_positive", 
                   "rho", "p_value", "q_value", "s_value", "lfsr", "lfdr")

# tests ----
test_that("calc_assoc works as expected", {
  res_1 <- calc_assoc(X, Y) # default matrix
  expect_is(res_1, "data.table")
  expect_true(all(c(res_col_names, "response") %in% names(res_1)))
  expect_equal(NROW(res_1), NROW(colnames(X)[colSums(X[rownames(Y), ]) > 0]))
  
  res_2 <- calc_assoc(X, Y_vec) # default vector
  expect_is(res_2, "data.table")
  expect_true(all(res_2$feature %in% colnames(X)[colSums(X[names(Y_vec), ]) > 0]))
  expect_equal(NROW(res_2), NROW(colnames(X)[colSums(X[names(Y_vec), ]) > 0]))
  
  # scenario: matrix Y has more than one column
  res_3 <- calc_assoc(X, Y_all)
  expect_is(res_3, "data.table")
  expect_true(all(c(res_col_names, "response") %in% names(res_3)))
  expect_equal(NROW(res_3), NROW(colnames(X)[colSums(X[rownames(Y_all), ]) > 0]) * NCOL(Y_all))
  expect_true(all(colnames(X) %in% res_3$feature))
  expect_true(all(colnames(Y_all) %in% res_3$response))
  
  expect_error(calc_assoc(X = matrix(LETTERS[1:9], nrow = 3), Y),
               "Assertion on 'X' failed: Must store numerics")
  expect_error(calc_assoc(X = matrix(1:9, nrow = 3), Y),
               "Must have names")
  expect_error(calc_assoc(X = X, Y = LETTERS[1:9]),
               "Assertion on 'Y' failed: Must inherit from class")
  expect_error(calc_assoc(X = X, Y = seq(0.25, 0.95, length.out = 8)),
               "Must have names")
  expect_error(calc_assoc(X = X, Y = seq(0.25, 0.95, length.out = 8)),
               "Must have names")
  expect_error(calc_assoc(X = X, Y = matrix(seq(0.25, 0.95, length.out = 16), ncol = 2)),
               "Must have names")
  
  # scenario: vector Y has no variance
  Y_vec_no_var <- rep(8, NROW(X))
  names(Y_vec_no_var) <- names(Y_vec)
  expect_warning({
    res_4 <- calc_assoc(X = X, Y = Y_vec_no_var)
  }, "Y has no variance")
  expect_is(res_4, "data.table")
  expect_true(all(res_col_names %in% names(res_4)))
  expect_true(all(res_4$feature %in% colnames(X)[colSums(X[names(Y_vec), ]) > 0]))
  expect_true(all(res_4[, lapply(.SD, is.na), .SDcols = -c("feature")]))
  expect_true(all(res_4[, lapply(.SD, is.numeric), .SDcols = -c("feature")]))
  
  # scenario: vector Y has too many NA
  Y_vec_NA <- Y_vec
  Y_vec_NA[1:6] <- NA
  expect_warning({
    expect_error({
      calc_assoc(X = X, Y = Y_vec_NA) 
    }, "Error: all input values are missing")
  }, "Missing values generated in calculation of cor. Likely cause: too many missing entries or zero variance.")
  
  # scenario: vector Y has not matching size
  expect_error({
    calc_assoc(X = X, Y = Y_vec[-1])
  }, "The X and Y dimensions must match.")
  
  # scenario: vector Y has not matching name (the same order)
  Y_vec_wrongname <- Y_vec 
  names(Y_vec_wrongname) <- sprintf("%s_fix", toupper(names(Y_vec_wrongname)))
  res_5 <- calc_assoc(X = X, Y = Y_vec_wrongname)
  expect_equal(res_5, calc_assoc(X, Y = Y_vec))
  
  # scenario: vector Y has wrong order
  Y_vec_wrongorder <- rev(Y_vec)
  res_6 <- calc_assoc(X = X, Y = Y_vec_wrongorder)
  expect_error(expect_equal(res_6, calc_assoc(X, Y = Y_vec)))
  
  # scenario: matrix X has no variance
  X_no_var <- X
  X_no_var[] <- 0.888
  expect_error({
    calc_assoc(X = X_no_var, Y = Y_vec)
  }, "Error: all input values are missing")
  
  # scenario: matrix X has NA columns
  X_NA_col <- X
  X_NA_col[, 2:3] <- NA
  expect_warning({
    res_7 <- calc_assoc(X = X_NA_col, Y = Y_vec)
  }, "NaNs produced")
  expect_true(all(res_7[, lapply(.SD, is.numeric), .SDcols = -c("feature")]))
  expect_false(all(colnames(X_NA_col)[2:3] %in% res_7$feature))
  
  # scenario: matrix X has to many NA
  X_NA <- X[, 1:4]
  X_NA[1:2, 3:4] <- NA
  X_NA[5:6, 1:2] <- NA
  X_NA[8, ] <- NA
  expect_error({
    calc_assoc(X = X_NA, Y = Y_vec)
  }, "all input values are missing")
  
  # scenario: matrix X has not matching size
  expect_error({
    calc_assoc(X = X[1:6, ], Y = Y_vec)
  }, "The X and Y dimensions must match.")
  
  # scenario: matrix Y has no variance
  Y_mat_no_var <- matrix(rep(0.88, NROW(X)), ncol = 1,
                         dimnames = list(rownames(Y), colnames(Y)))
  expect_warning({
    res_8 <- calc_assoc(X, Y = Y_mat_no_var)
  }, "The following columns in Y have no variance")
  expect_is(res_8, "data.table")
  expect_true(all(c(res_col_names, "response") %in% names(res_8)))
  expect_true(all(colnames(X) %in% res_8$feature))
  expect_true(all(res_8[, lapply(.SD, is.na), .SDcols = -c("feature", "response")]))
  expect_true(all(res_8[, lapply(.SD, is.numeric), .SDcols = -c("feature", "response")]))
  
  # scenario: one column in matrix Y has no variance
  Y_all_no_var <- Y_all
  Y_all_no_var[, 3] <- 0.888
  expect_warning({
    res_9 <- calc_assoc(X, Y = Y_all_no_var)
  }, "The following columns in Y have no variance")
  expect_is(res_9, "data.table")
  expect_true(all(c(res_col_names, "response") %in% names(res_9)))
  expect_true(all(res_9[, lapply(.SD, is.numeric), .SDcols = -c("feature", "response")]))
  expect_true(all(res_9[response == colnames(Y_all_no_var)[3],
                        lapply(.SD, is.na), .SDcols = -c("feature", "response")]))
  expect_false(all(res_9[response != colnames(Y_all_no_var)[3],
                         lapply(.SD, is.na), .SDcols = -c("feature", "response")]))
  
  # scenario: only one column in matrix Y has variance
  Y_all_no_var <- Y_all
  Y_all_no_var[, 3] <- 0.888
  Y_all_no_var[, 2] <- 0.8
  expect_warning({
    res_10 <- calc_assoc(X, Y = Y_all_no_var)
  }, "The following columns in Y have no variance")
  expect_is(res_10, "data.table")
  expect_true(all(c(res_col_names, "response") %in% names(res_10)))
  expect_true(all(res_10[, lapply(.SD, is.numeric), .SDcols = -c("feature", "response")]))
  expect_true(all(res_10[response == colnames(Y_all_no_var)[3],
                         lapply(.SD, is.na), .SDcols = -c("feature", "response")]))
  expect_false(all(res_10[response != colnames(Y_all_no_var)[3],
                          lapply(.SD, is.na), .SDcols = -c("feature", "response")]))
  
  # scenario: matrix Y has no matching size
  expect_error({
    calc_assoc(X, Y = Y[1:6, , drop = FALSE])
  }, "The X and Y dimensions must match.")
  expect_error({
    calc_assoc(X, Y = Y_all[1:6, , drop = FALSE])
  }, "The X and Y dimensions must match.")
  
  # scenario: matrix Y has not matching name (the same order)
  Y_wrongname <- Y
  rownames(Y_wrongname) <- sprintf("%s_fix", toupper(rownames(Y_wrongname)))
  res_11 <- calc_assoc(X = X, Y = Y_wrongname)
  expect_is(res_11, "data.table")
  expect_equal(res_11, calc_assoc(X, Y = Y))
  
  # scenario: matrix Y has wrong order
  Y_wrongorder <- Y[rev(seq_along(Y)), , drop = FALSE]
  res_12 <- calc_assoc(X = X, Y = Y_wrongorder)
  expect_is(res_12, "data.table")
  expect_error(expect_equal(res_12, calc_assoc(X, Y = Y)))
  
  # scenario: matrix Y has no colname
  Y_no_colname <- Y
  colnames(Y_no_colname) <- NULL
  res_13 <- calc_assoc(X = X, Y = Y_no_colname)
  expect_is(res_13, "data.table")
  expect_equal(res_13[, -"response", with = FALSE], 
               calc_assoc(X, Y = Y)[, -"response", with = FALSE])
  expect_true(all(sprintf("var_%s", seq_len(NCOL(Y))) %in% res_13$response))
  
  # scenario: matrix Y has too many NA
  Y_NA_row <- Y
  Y_NA_row[1:6] <- NA
  expect_warning({
    expect_error({
      calc_assoc(X, Y = Y_NA_row)
    }, "all input values are missing")
  }, "Missing values generated in calculation of cor. Likely cause: too many missing entries or zero variance.")
  
  # scenario: one column in matrix Y has one column with NA
  Y_all_NA <- Y_all
  Y_all_NA[, 3] <- NA
  expect_error({
    expect_warning({
      calc_assoc(X, Y = Y_all_NA)
    }, "NaNs produced")
  }, "all input values are missing")
  
  # scenario: matrix Y has some NA
  Y_all_NA_col <- Y_all
  Y_all_NA_col[1:2, 3] <- NA
  Y_all_NA_col[5:7, 1:2] <- NA
  expect_error({
    expect_warning({
      calc_assoc(X, Y = Y_all_NA_col)
    }, "Missing values generated in calculation of cor. Likely cause: too many missing entries or zero variance.")
  }, "all input values are missing")
  
  # scenario: matrix X and Y are too small (min row number has to be 6)
  expect_error({
    calc_assoc(X[1:5, ], Y[1:5, , drop = FALSE])
  }, "all input values are missing")
  
  # scenario: X and Y with only one not NA p.val and lack of feature name
  res_14 <- calc_assoc(X = X_feat, Y = Y_feat)
  expect_is(res_14, "data.table")
  expect_true(all(rownames(stats::na.omit(t(X_feat))) %in% res_14$feature))
  expect_true(all(colnames(Y_feat) %in% res_14$response))

  # scenario: X and Y are big matrix
  expect_warning({
    res_15 <- calc_assoc(X = X_big, Y = Y_big)
  }, "The following columns in Y have no variance")
  expect_is(res_15, "data.table")
  expect_true(all(res_15[, lapply(.SD, is.numeric), .SDcols = -c("feature", "response")]))
  expect_true(all(colnames(X_big) %in% res_15$feature))
  expect_true(all(colnames(Y_big) %in% res_15$response))
})

test_that(".order_assoc_result works as expected", {
  tab_assoc <- data.table::data.table(
    betahat = withr::with_seed(42, sample(seq(-4.2, 2.1, 0.35), 4)),
    sebetahat = withr::with_seed(42, sample(seq(0, 1.42, 0.13), 4)),
    NegativeProb = withr::with_seed(42, sample(seq(0.1, 0.9, 0.1), 4)),
    PositiveProb = withr::with_seed(314, sample(seq(0.05, .95, 0.08), 4)),
    lfsr = withr::with_seed(271, sample(seq(0.1, 0.95, 0.03), 4)),
    svalue = withr::with_seed(42, sample(seq(0.01, 0.33, 0.08), 4)),
    lfdr = withr::with_seed(271, sample(seq(0, 1, 0.035), 4)),
    qvalue = withr::with_seed(981, sample(seq(0.01, 0.98, 0.07), 4)),
    PosteriorMean = withr::with_seed(42, sample(seq(-4.2, 2.1, 0.35), 4)) / 2,
    PosteriorSD = withr::with_seed(271, sample(seq(0.95, 1.35, 0.05), 4)),
    dep.var = rep("met_1", 1),
    ind.var = sprintf("feat_%s", 1:4),   
    p.val = withr::with_seed(314, sample(seq(0.01, .55, 0.13), 4)),
    rho = withr::with_seed(42, sample(seq(-1, 1, 0.35), 4))
  )
  
  res <- .order_assoc_result(tab_assoc)
  expect_true(all(names(res) %in% c(res_col_names, "response")))
  
  expect_error(.order_assoc_result(unlist(tab_assoc)),
               "Assertion on 'res_dt' failed: Must be a data.table")
})
