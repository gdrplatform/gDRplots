context("Test fun utils")

test_that("isColDark works as expected", {
  color_names <- c("#33cc33", "#d6f5d6", "#0a290a")
  expect_false(isColDark(color_names[1]))
  expect_false(isColDark(color_names[2]))
  expect_true(isColDark(color_names[3]))
  
  color_names <- c("darkblue", "yellow", "tomato")
  expect_true(isColDark(color_names[1]))
  expect_false(isColDark(color_names[2]))
  expect_false(isColDark(color_names[3]))
  
  expect_error(isColDark(1), "Must be of type 'string'")
  expect_error(isColDark("nice pink"), "Must be valid color name")
})


test_that("getColLuminance works as expected", {
  color_names <- c("#33cc33", "#d6f5d6", "#0a290a")
  expect_equal(getColLuminance(color_names[1]), 0.429871, tolerance = 1e-5)
  expect_equal(getColLuminance(color_names[2]), 0.839746, tolerance = 1e-5)
  expect_equal(getColLuminance(color_names[3]), 0.016340, tolerance = 1e-5)
  
  color_names <- c("darkblue", "yellow", "tomato")
  expect_equal(getColLuminance(color_names[1]), 0.0186408, tolerance = 1e-5)
  expect_equal(getColLuminance(color_names[2]), 0.9278000, tolerance = 1e-5)
  expect_equal(getColLuminance(color_names[3]), 0.3238907, tolerance = 1e-5)
  
  expect_error(getColLuminance(1), "Must be of type 'string'")
  expect_error(getColLuminance("nice pink"), "Must be valid color name")
})


test_that("isValidColor works as expected", {
  color_names <- c("#33cc33", "#d6f5d6", "#0a290a", "#F9B42DFF", "#714D6932", "#C2F970DC")
  expect_true(all(vapply(color_names, isValidColor, logical(1))))
  
  color_names <- c("darkblue", "yellow", "tomato")
  expect_true(all(vapply(color_names, isValidColor, logical(1))))
  
  color_names <- c("nice pink", "RED", "#0a290", "#C2F970D")
  expect_false(all(vapply(color_names, isValidColor, logical(1))))
  
  expect_error(isValidColor(1), "Must be of type 'string'")
  expect_error(isValidColor(NULL), "Must be of type 'string'")
  expect_error(isValidColor(NA), "Assertion on 'col_name' failed: May not be NA.")
})

test_that("colorToHex works as expected", {
  color_names <- c("orange", "darkblue", "green", "lavenderblush", "gray66", "slategray2", "tomato")
  expect_identical(
    vapply(color_names, colorToHex, character(1), USE.NAMES = FALSE),
    c("#FFA500", "#00008B", "#00FF00", "#FFF0F5", "#A8A8A8", "#B9D3EE", "#FF6347"))
  
  expect_error(colorToHex(1), "Must be of type 'string'")
  expect_error(colorToHex(NULL), "Must be of type 'string'")
  expect_error(colorToHex("pinki"), "Must be valid color name")
})
