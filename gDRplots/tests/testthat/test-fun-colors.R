context("Test fun colors")

test_that("brew_palette works as expected", {
  pal_name <- sample(
    x = c("Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3"), size = 1)
  pal_col <- RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[pal_name, ]$maxcolors, pal_name)
  
  n_s <- 2
  small_pal <- brew_palette(n_s, pal_name)
  expect_equal(small_pal, pal_col[1:n_s])
  
  n_n <- NROW(pal_col) - 1
  normal_pal <- brew_palette(n_n, pal_name)
  expect_equal(normal_pal, pal_col[1:n_n])
  
  n_l <- ceiling(NROW(pal_col) * 1.5)
  long_pal <- brew_palette(n_l, pal_name)
  expect_equal(long_pal, rep(pal_col, length.out = n_l))
  
  shuffled_normal_pal <- brew_palette(n_n, pal_name, shuffle = TRUE)
  expect_false(identical(normal_pal, shuffled_normal_pal))
  expect_identical(sort(normal_pal), sort(shuffled_normal_pal))
  
  expect_error(brew_palette(n = "str", name = "Accent"), 
               "Assertion on 'n' failed: Must be of type 'number', not 'character'.")
  expect_error(brew_palette(n = 0, name = "Accent"),
               "Assertion on 'n' failed: Element 1 is not >= 1.")
  expect_error(brew_palette(n = 3, name = 1),
               "Assertion on 'name' failed: Must be of type 'string', not 'double'")
  expect_error(brew_palette(n = 3, name = "str"),
               "Assertion on 'name' failed: Must be element of set")
  expect_error(brew_palette(n = 3, name = "Accent", shuffle = 1),
               "Assertion on 'shuffle' failed: Must be of type 'logical', not 'double'.")
})


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
  expect_error(get_col_luminance("nice pink"), "Must be valid color name")
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