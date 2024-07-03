context("Test utils")

# Example ggplot object

p1 <- ggplot2::ggplot(ggplot2::mpg, ggplot2::aes(displ, hwy, color = class)) + ggplot2::geom_point()

# Example pheatmap object
mat <- matrix(rnorm(100), 10, 10)
p2 <- pheatmap::pheatmap(mat)

test_that("estimate_plot_size works for ggplot object", {
  size <- estimate_plot_size(p1)
  expect_named(size, c("width", "height"))
  expect_equal(size[["width"]], 13.5)
  expect_equal(size[["height"]], 9.5)
})

test_that("estimate_plot_size works for pheatmap object", {
  size <- estimate_plot_size(p2)
  expect_named(size, c("width", "height"))
  expect_equal(size[["width"]], 15)
  expect_equal(size[["height"]], 11)
})

test_that("estimate_plot_size handles invalid inputs", {
  invalid_plot <- list()
  expect_error(estimate_plot_size(invalid_plot),
               "Assertion on 'plt' failed: Must inherit from class 'ggplot'/'pheatmap', but has class 'list'.")
  
  invalid_base_width <- -5
  expect_error(estimate_plot_size(p1, base_width = invalid_base_width),
               "Assertion on 'base_width' failed: Element 1 is not >= 0.")
  
  invalid_base_height <- -5
  expect_error(estimate_plot_size(p1, base_height = invalid_base_height),
               "Assertion on 'base_height' failed: Element 1 is not >= 0.")
  
  invalid_scale_factor <- -0.5
  expect_error(estimate_plot_size(p1, scale_factor = invalid_scale_factor),
               "Assertion on 'scale_factor' failed: Element 1 is not >= 0.")
})
