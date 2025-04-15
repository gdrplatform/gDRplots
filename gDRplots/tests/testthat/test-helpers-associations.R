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

# association columns ----
res_col_names <- c("feature", "est_beta", "est_beta_se", 
                   "posterior_mean", "posterior_sd", "prob_negative", "prob_positive", 
                   "rho", "p_value", "q_value", "s_value", "lfsr", "lfdr")

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
  
  expect_error(calc_assoc(X = matrix(LETTERS[1:9], nrow = 3), Y),
               "Assertion on 'X' failed: Must store numerics")
  expect_error(calc_assoc(X = matrix(1:9, nrow = 3), Y),
               "Must have names")
  expect_error(calc_assoc(X = X, Y = LETTERS[1:9]),
               "Assertion on 'Y' failed: Must inherit from class")
  expect_error(calc_assoc(X = X, Y = seq(0.25, 0.95, length.out = 8)),
               "Must have names")
})

test_that(".calc_assoc_vector works as expected", {
  # scenario: vector Y has no variance 
  Y_vec_no_var <- rep(8, NROW(X))
  names(Y_vec_no_var) <- names(Y_vec)
  expect_warning({
    res_1 <- .calc_assoc_vector(X = X, Y = Y_vec_no_var)
  }, "Y has no variance")
  expect_is(res_1, "data.table")
  expect_true(all(res_col_names %in% names(res_1)))
  expect_true(all(res_1[, lapply(.SD, is.na), .SDcols = -c("feature")]))
  expect_true(all(res_1[, lapply(.SD, is.numeric), .SDcols = -c("feature")]))
  
  # scenario: vector Y has too many NA
  Y_vec_NA <- Y_vec
  Y_vec_NA[1:6] <- NA
  expect_warning({
    expect_error({
      .calc_assoc_vector(X = X, Y = Y_vec_NA) 
    }, "Error: all input values are missing")
  }, "Missing values generated in calculation of cor. Likely cause: too many missing entries or zero variance.")
  
  # scenario: vector Y has not matching size
  expect_error({
    .calc_assoc_vector(X = X, Y = Y_vec[-1])
  }, "non-conformable arguments")
  
  # scenario: vector Y has not matching name (the same order)
  Y_vec_wrongname <- Y_vec 
  names(Y_vec_wrongname) <- sprintf("%s_fix", toupper(names(Y_vec_wrongname)))
  res_2 <- .calc_assoc_vector(X = X, Y = Y_vec_wrongname)
  expect_equal(res_2, .calc_assoc_vector(X, Y = Y_vec))
  
  # scenario: vector Y has wrong order
  Y_vec_wrongorder <- rev(Y_vec)
  res_3 <- .calc_assoc_vector(X = X, Y = Y_vec_wrongorder)
  expect_error(expect_equal(res_3, .calc_assoc_vector(X, Y = Y_vec)))
  
  # scenario: vector Y has no name
  Y_vec_noname <- Y_vec 
  names(Y_vec_noname) <- NULL
  expect_error(.calc_assoc_vector(X = X, Y = Y_vec_noname), 
               "Must have names")
  
  # scenario: vector Y in not numeric
  expect_error(.calc_assoc_vector(X = X, Y = LETTERS[1:9]), 
               "Assertion on 'Y' failed: Must be of type 'numeric'")
  
  # scenario: matrix X has no variance
  X_no_var <- X
  X_no_var[] <- 0.888
  expect_error({
    .calc_assoc_vector(X = X_no_var, Y = Y_vec)
  }, "Error: all input values are missing")
  
  # scenario: matrix X has NA
  X_NA_col <- X
  X_NA_col[, 2:3] <- NA
  expect_warning({
    res_4 <- .calc_assoc_vector(X = X_NA_col, Y = Y_vec)
  }, "NaNs produced")
  expect_true(all(res_4[, lapply(.SD, is.numeric), .SDcols = -c("feature")]))
  expect_false(all(colnames(X_NA_col)[2:3] %in% res_4$feature))
  
  # scenario: matrix X has to many NA
  X_NA <- X[, 1:4]
  X_NA[1:2, 3:4] <- NA
  X_NA[5:6, 1:2] <- NA
  X_NA[8, ] <- NA
  expect_error({
    .calc_assoc_vector(X = X_NA, Y = Y_vec)
  }, "all input values are missing")
  
  # scenario: matrix X has not matching size
  expect_error({
    .calc_assoc_vector(X = X[1:6, ], Y = Y_vec)
  }, "non-conformable arguments")
  
})

test_that(".calc_assoc_matrix works as expected", {
  # scenario: matrix Y has no variance 
  Y_mat_no_var <- matrix(rep(0.88, NROW(X)), ncol = 1,
                         dimnames = list(rownames(Y), colnames(Y)))
  expect_warning({
    res_1 <- .calc_assoc_matrix(X, Y = Y_mat_no_var) 
  }, "The following columns in Y have no variance")
  expect_is(res_1, "data.table")
  expect_true(all(c(res_col_names, "response") %in% names(res_1)))
  expect_true(all(res_1[, lapply(.SD, is.na), .SDcols = -c("feature", "response")]))
  expect_true(all(res_1[, lapply(.SD, is.numeric), .SDcols = -c("feature", "response")]))
  
  # scenario: one column in matrix Y has no variance 
  Y_all_no_var <- Y_all
  Y_all_no_var[, 3] <- 0.888
  expect_warning({
    res_2 <- .calc_assoc_matrix(X, Y = Y_all_no_var)
  }, "The following columns in Y have no variance")
  expect_is(res_2, "data.table")
  expect_true(all(c(res_col_names, "response") %in% names(res_2)))
  expect_true(all(res_2[, lapply(.SD, is.numeric), .SDcols = -c("feature", "response")]))
  expect_true(all(res_2[response == colnames(Y_all_no_var)[3], 
                        lapply(.SD, is.na), .SDcols = -c("feature", "response")]))
  expect_false(all(res_2[response != colnames(Y_all_no_var)[3], 
                         lapply(.SD, is.na), .SDcols = -c("feature", "response")]))
  
  # scenario: matrix Y has no matching size
  expect_error({
    .calc_assoc_matrix(X, Y = Y[1:6, , drop = FALSE])
  }, "non-conformable arguments")
  expect_error({
    .calc_assoc_matrix(X, Y = Y_all[1:6, , drop = FALSE])
  }, "non-conformable arguments")
  
  # scenario: matrix Y has not matching name (the same order)
  Y_wrongname <- Y
  rownames(Y_wrongname) <- sprintf("%s_fix", toupper(rownames(Y_wrongname)))
  res_3 <- .calc_assoc_matrix(X = X, Y = Y_wrongname)
  expect_is(res_3, "data.table")
  expect_equal(res_3, .calc_assoc_matrix(X, Y = Y))
  
  # scenario: matrix Y has wrong order
  Y_wrongorder <- Y[rev(seq_along(Y)), , drop = FALSE]
  res_4 <- .calc_assoc_matrix(X = X, Y = Y_wrongorder)
  expect_is(res_4, "data.table")
  expect_error(expect_equal(res_4, .calc_assoc_matrix(X, Y = Y)))
  
  # scenario: matrix Y has no name
  Y_noname <- Y
  rownames(Y_noname) <- NULL
  expect_error(.calc_assoc_matrix(X = X, Y = Y_noname), 
               "Must have names")
  
  # scenario: matrix Y has too many NA
  Y_NA_row <- Y
  Y_NA_row[1:6] <- NA
  expect_warning({
    expect_error({
      .calc_assoc_matrix(X, Y = Y_NA_row)
    }, "all input values are missing")
  }, "Missing values generated in calculation of cor. Likely cause: too many missing entries or zero variance.")
  
  # scenario: one column in matrix Y has one column with NA
  Y_all_NA <- Y_all
  Y_all_NA[, 3] <- NA
  expect_error({
    expect_warning({
      .calc_assoc_matrix(X, Y = Y_all_NA)
    }, "NaNs produced")
  }, "all input values are missing")
  
  # scenario: matrix Y has some NA
  Y_all_NA_col <- Y_all
  Y_all_NA_col[1:2, 3] <- NA
  Y_all_NA_col[5:7, 1:2] <- NA
  expect_error({
    expect_warning({
      .calc_assoc_matrix(X, Y = Y_all_NA_col)
    }, "Missing values generated in calculation of cor. Likely cause: too many missing entries or zero variance.")
  }, "all input values are missing")
  
  # scenario: matrix X and Y are too small
  expect_error({
    .calc_assoc_matrix(X[1:5, ], Y[1:5, , drop = FALSE])
  }, "all input values are missing")
  
  # scenario: X and Y with only one not NA p.val
  X_1 <- matrix(c(rep(NA, 2), rep(1, 10), rep(NA, 2), rep(1, 2), rep(NA, 3), rep(2, 5), rep(1, 2), rep(3, 2)), nrow = 7,
                dimnames = list(sprintf("row_%s", 1:7), sprintf("feat_%s", 1:4)))
  Y_1 <- matrix(c(10:16), ncol = 1,
                dimnames = list(sprintf("row_%s", 1:7), "met"))
  res_5 <- .calc_assoc_matrix(X = X_1, Y = Y_1)
  expect_is(res_5, "data.table")
  
  # scenario: lack of feature name
  X_f <- matrix(c(1:6, rep(1, 6)), nrow = 6,
                dimnames = list(sprintf("row_%s", 1:6), sprintf("feat_%s", 1:2)))
  Y_f <- matrix(c(10:15), ncol = 1,
                dimnames = list(sprintf("row_%s", 1:6), "met"))
  res_6 <- .calc_assoc_matrix(X = X_f, Y = Y_f)
  expect_is(res_6, "data.table")
})

test_that(".order_assoc_result works as expected", {
  tab_assoc <- data.table::data.table(
    betahat = withr::with_seed(42, sample(seq(-4.2, 2.1, 0.35), 4)),
    sebetahat = withr::with_seed(42, sample(seq(0, 1.42, 0.13), 4)),
    NegativeProb = withr::with_seed(42, sample(seq(0.1, 0.9, 0.1), 4)),
    PositiveProb = withr::with_seed(314, sample(seq(0.05, .95, 0.08), 4)),
    lfsr = withr::with_seed(271, sample(seq(0.1, 0.95, 0.03), 4)),
    svalue = withr::with_seed(42, sample(seq(0.01, 0.33, 0.08), 4)),
    lfdr =withr::with_seed(271, sample(seq(0, 1, 0.035), 4)),
    qvalue = withr::with_seed(981, sample(seq(0.01, 0.98, 0.07), 4)),
    PosteriorMean = withr::with_seed(42, sample(seq(-4.2, 2.1, 0.35), 4)) / 2,
    PosteriorSD = withr::with_seed(271, sample(seq(0.95, 1.35, 0.05), 4)),
    dep.var = rep("met_1", 1),
    ind.var = sprintf("feat_%s", 1:4),   
    p.val= withr::with_seed(314, sample(seq(0.01, .55, 0.13), 4)),
    rho = withr::with_seed(42, sample(seq(-1, 1, 0.35), 4))
  )
  
  res <- .order_assoc_result(tab_assoc)
  expect_true(all(names(res) %in% c(res_col_names, "response")))
  
  expect_error(.order_assoc_result(unlist(tab_assoc)),
               "Assertion on 'res_dt' failed: Must be a data.table")
})
