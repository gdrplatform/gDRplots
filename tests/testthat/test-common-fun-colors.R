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

test_that("get_iso_colors works as expected", {
  ### expected values
  gic <- get_iso_colors()
  expect_true(length(gic) > 2)
  expect_identical(class(gic), "character")
  gic2 <- get_iso_colors(formals(get_iso_colors)[[1]][[3]])
  expect_true(any(gic != gic2))
  expect_identical(length(gic), length(gic2))

  ### errors
  expect_error(get_iso_colors("inv_param"), "'arg' should be one of ")
})

test_that("assert_RGB_format as expected", {
  color_vector <- c(25, 56, 189)
  expect_equal(assert_RGB_format(color_vector), NULL)
  color_vector <- c(201, 128, 352)
  expect_error(assert_RGB_format(color_vector),
               "Some value is greater than 255. Not valid RGB format.")
})

test_that("get_qual_colors works as expected", {
  max_len <- sum(RColorBrewer::brewer.pal.info[
    RColorBrewer::brewer.pal.info$category == "qual" &
      RColorBrewer::brewer.pal.info$colorblind == TRUE, ]$maxcolors)

  expect_equal(NROW(get_qual_colors()), max_len)
  expect_equal(NROW(unique(get_qual_colors())), max_len)
  expect_equal(get_qual_colors(0), "#000000")
  expect_equal(get_qual_colors(5), c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"))
  expect_equal(NROW(get_qual_colors(42)), 42)
  expect_equal(NROW(unique(get_qual_colors(42))), 42)

  expect_error(get_qual_colors("one"),
               "Assertion on 'n' failed: Must be of type 'single integerish value'")
  expect_error(get_qual_colors(c(2, 3)),
               "Assertion on 'n' failed: Must have length 1.")
  expect_error(get_qual_colors(2.3),
               "Assertion on 'n' failed: Must be of type 'single integerish value'")
  expect_error(get_qual_colors(-1),
               "Assertion on 'n' failed: Element 1 is not >= 0.")
})


test_that("is_color_dark works as expected", {
  color_names <- c("#33cc33", "#d6f5d6", "#0a290a")
  expect_false(is_color_dark(color_names[1]))
  expect_false(is_color_dark(color_names[2]))
  expect_true(is_color_dark(color_names[3]))

  color_names <- c("darkblue", "yellow", "tomato")
  expect_true(is_color_dark(color_names[1]))
  expect_false(is_color_dark(color_names[2]))
  expect_false(is_color_dark(color_names[3]))

  expect_error(is_color_dark(1), "Must be of type 'string'")
  expect_error(is_color_dark("nice pink"), "Must be a valid color name")
})


test_that("get_col_luminance works as expected", {
  color_names <- c("#33cc33", "#d6f5d6", "#0a290a")
  expect_equal(get_col_luminance(color_names[1]), 0.429871, tolerance = 1e-5)
  expect_equal(get_col_luminance(color_names[2]), 0.839746, tolerance = 1e-5)
  expect_equal(get_col_luminance(color_names[3]), 0.016340, tolerance = 1e-5)

  color_names <- c("darkblue", "yellow", "tomato")
  expect_equal(get_col_luminance(color_names[1]), 0.0186408, tolerance = 1e-5)
  expect_equal(get_col_luminance(color_names[2]), 0.9278000, tolerance = 1e-5)
  expect_equal(get_col_luminance(color_names[3]), 0.3238907, tolerance = 1e-5)

  expect_error(get_col_luminance(1), "Must be of type 'string'")
  expect_error(get_col_luminance("nice pink"), "Must be a valid color name")
})
