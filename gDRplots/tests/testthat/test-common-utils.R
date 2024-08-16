test_that("create_log_seq works as expected", {
  sequence <- c(1, 1.18920711500272, 1.4142135623731, 1.68179283050743, 2)
  
  # sequence is generated
  expect_equal(create_log_seq(1, 2, 5), sequence)
  # sequence has proper length
  expect_equal(length(create_log_seq(1, 2, 5)), 5)
  # sequence has proper limits
  expect_equal(create_log_seq(1, 2, 5)[1], 1)
  expect_equal(create_log_seq(1, 2, 5)[length(create_log_seq(1, 2, 5))], 2)
  # sequence grows linearly in log domain
  s <- create_log_seq(1, 2, 5)
  sLog <- log10(s)
  differences <- diff(sLog)
  differences <- signif(differences, 10)
  expect_true(length(unique(differences)) == 1L)
  
  expect_error(create_log_seq(NULL, 2, 5))
  expect_error(create_log_seq(-2, 2, 5))
  expect_error(create_log_seq(1, -5, 5))
  expect_error(create_log_seq(1, "str", 5))
  expect_error(create_log_seq(1, 2, TRUE))
  expect_error(create_log_seq(1, 2, 0))
})