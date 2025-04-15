context("Test helpers-associations")

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
  expect_error(calc_assoc(X = X, Y = dt_response_met$RV_gDR_x_mean),
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
  },"Missing values generated in calculation of cor. Likely cause: too many missing entries or zero variance.")
  
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
               "Assertion on 'Y' failed: Must inherit from class 'numeric'")
  
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
  X_NA[1:2, 4] <- NA
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
                         dimnames = list(rownames(Y), selected_metric))
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
})
