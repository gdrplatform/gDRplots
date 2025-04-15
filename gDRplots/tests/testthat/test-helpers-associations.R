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
  
  # scenario: one column in matrix Y has no variance 
  Y_all_no_var <- Y_all
  Y_all_no_var[, 3] <- 0.888
  expect_warning({
    res_6 <- calc_assoc(X, Y = Y_all_no_var)
  }, "The following columns in Y have no variance")
  expect_is(res_6, "data.table")
  expect_true(all(c(res_col_names, "response") %in% names(res_6)))
  expect_true(all(res_6[, lapply(.SD, is.numeric), .SDcols = -c("feature", "response")]))
})
