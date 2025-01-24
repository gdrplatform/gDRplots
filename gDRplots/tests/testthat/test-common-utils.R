context("Test common utils")

test_that("compute_distances works as expected", {
  set.seed(1234)
  
  # matrix without row names
  y <- matrix(LETTERS[1:9], nrow = 3, ncol = 3, dimnames = list(letters[1:3]))
  expect_error(compute_distances(y),
               "Assertion on 'x' failed: Must be of type 'numeric'")
  y <- matrix(seq_len(9), nrow = 3, ncol = 3)
  rownames(y) <- letters[seq_len(nrow(y))]
  expect_error(compute_distances(y, method = "str"),
               "Assertion on 'method' failed: Must be element of set")
  expect_error(compute_distances(y, use = "str"),
               "Assertion on 'use' failed: Must be element of set")
  expect_error(compute_distances(y, stand = "str"),
               "Assertion on 'stand' failed: Must be of type 'logical flag', not 'character'.")
  expect_error(compute_distances(y, dummy = "str"),
               "Assertion on 'dummy' failed: Must be of type 'number', not 'character'.")
  
  # all rows with variance
  x <- matrix(seq_len(9), nrow = 3, ncol = 3)
  rownames(x) <- letters[seq_len(nrow(x))]
  expect_equal(as.numeric(compute_distances(x)), rep(0, 3))
  expect_equal(as.numeric(compute_distances(x, stand = TRUE)), rep(1, 3))
  diag(x) <- NA
  expect_equal(as.numeric(compute_distances(x)), rep(0, 3))
  expect_equal(as.numeric(compute_distances(x, dummy = NA)), as.numeric(rep(NA, 3)))
  expect_equal(as.numeric(compute_distances(x, dummy = -Inf)), rep(-Inf, 3))
  
  xx <- matrix(-4:4, nrow = 3, ncol = 3)
  rownames(xx) <- letters[seq_len(nrow(xx))]
  expect_equal(as.numeric(compute_distances(xx)), rep(0, 3))
  expect_equal(as.numeric(compute_distances(xx, stand = TRUE)), rep(1, 3))
  expect_equal(as.numeric(compute_distances(xx, method = "manhattan")), c(3, 6, 3)) 
  
  # single row with no variance
  a <- t(matrix(c(rep(1, 3), sample(seq_len(100), 15)), nrow = 3, ncol = 6))
  rownames(a) <- letters[seq_len(nrow(a))]
  expect_equal(sort(unique(as.numeric(compute_distances(a)))), seq(0, 2, 0.5))
  
  # single row with variance
  b <- t(matrix(c(rep(1, 15), sample(seq_len(100), 3)), nrow = 3, ncol = 6))
  rownames(b) <- letters[seq_len(nrow(b))]
  expect_equal(as.numeric(compute_distances(b)), rep(1, 15))
  
  # single row with variance - small
  y <- matrix(c(1, 1, 2, 1), nrow = 2, ncol = 2)
  rownames(y) <- letters[seq_len(nrow(y))]
  expect_equal(as.numeric(compute_distances(y)), 1)
  
  # all rows with variance
  z <- matrix(sample(seq_len(100), 9), nrow = 3, ncol = 3)
  rownames(z) <- letters[seq_len(nrow(z))]
  expect_equal(as.numeric(compute_distances(z)), c(1.5, 1.5, 0))
  
  # all rows with no variance
  y <- matrix(1, nrow = 3, ncol = 3)
  rownames(y) <- letters[seq_len(nrow(y))]
  expect_equal(as.numeric(compute_distances(y)), rep(1, 3))
  
})

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