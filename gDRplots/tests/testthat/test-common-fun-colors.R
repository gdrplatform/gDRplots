context("Test common fun colors")


test_that("is_valid_color works as expected", {
  color_names <- c("#33cc33", "#d6f5d6", "#0a290a", "#F9B42DFF", "#714D6932", "#C2F970DC")
  expect_true(all(vapply(color_names, is_valid_color, logical(1))))
  
  color_names <- c("darkblue", "yellow", "tomato")
  expect_true(all(vapply(color_names, is_valid_color, logical(1))))
  
  color_names <- c("nice pink", "RED", "#0a290", "#C2F970D")
  expect_false(all(vapply(color_names, is_valid_color, logical(1))))
  
  expect_error(is_valid_color(1), "Must be of type 'string'")
  expect_error(is_valid_color(NULL), "Must be of type 'string'")
  expect_error(is_valid_color(NA), "Assertion on 'col_name' failed: May not be NA.")
})

test_that("get_iso_colors", {
  ### expected values
  gic <- get_iso_colors()
  expect_true(length(gic) > 2)
  expect_identical("character", class(gic))
  gic2 <- get_iso_colors(formals(get_iso_colors)[[1]][[3]])
  expect_true(any(gic != gic2))
  expect_identical(length(gic), length(gic2))
  
  ### errors
  expect_error(get_iso_colors("inv_param"), "'arg' should be one of ")
})

test_that("assert_RGB_format", {
  color_vector <- c(25, 56, 189)
  expect_equal(assert_RGB_format(color_vector), NULL)
  color_vector <- c(201, 128, 352)
  expect_error(assert_RGB_format(color_vector), 
               "Some value is greater than 255. Not valid RGB format.")
})