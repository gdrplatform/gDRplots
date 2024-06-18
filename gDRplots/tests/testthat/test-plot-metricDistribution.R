context("Test metricDistribution plot")

m <- 9
n <- 5

synthetic_data <- gDRutils::gen_synthetic_data(m, n)
drug_names <- synthetic_data$drug_names
cell_names <- synthetic_data$cell_names
dt <- synthetic_data$dt

axis_limits <- "fixed"
points <- "off"
show_tick_labels <- TRUE
type <- "box"
var_grp <- "none"
var_col <- "none"
var_x <- "Cell Line Name"
var_y <- "GR_AOC"

test_that("'get_axis_min'", {
  my_data <- c(0.11, 0.01, -0.08)
  expect_error(get_axis_min(list()),
               regexp = "Must be of type 'numeric', not 'list'",
               fixed = TRUE)
  expect_error(get_axis_min(my_data, nearest_shift = c(1, 2)),
               regexp = "Assertion on 'nearest_shift' failed: Must have length 1")
  expect_identical(get_axis_min(my_data), -0.1)
  expect_identical(get_axis_min(my_data, nearest_shift = 0.5), -0.5)
  
  expect_error(get_axis_max(list()),
               regexp = "Must be of type 'numeric', not 'list'",
               fixed = TRUE)
  expect_error(get_axis_max(my_data, nearest_shift = c(1, 2)),
               regexp = "Assertion on 'nearest_shift' failed: Must have length 1")
  expect_identical(get_axis_max(my_data), 0.2)
  expect_identical(get_axis_max(my_data, nearest_shift = 0.5), 0.5)
})

# plotly_metric_distribution tests
test_that("check output type and data",  {
  p <- plotly_metric_distribution(dt, var_x, var_y, var_grp, var_col, type, points, axis_limits, show_tick_labels)
  expect_type(p, "list")
  expect_identical(p$x$attrs[[1]]$type,
                   "box")
  expect_identical(p$x$attrs[[1]]$x,
                   rep(cell_names, each = m))
  expect_identical(p$x$attrs[[1]]$y,
                   dt$GR_AOC)
})

test_that("check input arguments", {
  expect_error(
    plotly_metric_distribution(n, var_x, var_y, var_grp, var_col, type, points, axis_limits, show_tick_labels),
    "Assertion on 'data' failed: Must be a data.table, not double."
  )
  expect_error(
    plotly_metric_distribution(dt, n, var_y, var_grp, var_col, type, points, axis_limits, show_tick_labels),
    "Assertion on 'var_x' failed: Must be of type 'string', not 'double'."
  )
  expect_error(
    plotly_metric_distribution(dt, var_x, n, var_grp, var_col, type, points, axis_limits, show_tick_labels),
    "Assertion on 'var_y' failed: Must be of type 'string', not 'double'."
  )
  expect_error(
    plotly_metric_distribution(dt, var_x, var_y, n, var_col, type, points, axis_limits, show_tick_labels),
    "Assertion on 'var_grp' failed: Must be of type 'string', not 'double'."
  )
  expect_error(
    plotly_metric_distribution(dt, var_x, var_y, var_grp, n, type, points, axis_limits, show_tick_labels),
    "Assertion on 'var_col' failed: Must be of type 'string', not 'double'."
  )
  expect_error(
    plotly_metric_distribution(dt, var_x, var_y, var_grp, var_col, n, points, axis_limits, show_tick_labels),
    "Assertion on 'type' failed: Must be of type 'string', not 'double'."
  )
  expect_error(
    plotly_metric_distribution(dt, var_x, var_y, var_grp, var_col, type, n, axis_limits, show_tick_labels),
    "Assertion on 'points' failed: Must be of type 'string', not 'double'."
  )
  expect_error(
    plotly_metric_distribution(dt, var_x, var_y, var_grp, var_col, type, points, n, show_tick_labels),
    "'arg' must be NULL or a character vector"
  )
})
