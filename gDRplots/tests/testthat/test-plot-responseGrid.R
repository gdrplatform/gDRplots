context("Test responseGrid plot")

m <- 9
n <- 5

synthetic_data <- gDRutils::gen_synthetic_data(m, n)
drug_names <- synthetic_data$drug_names
cell_names <- synthetic_data$cell_names
dt <- synthetic_data$dt

var_y <- "GR value"

prepared_curves <- prepareCurves(dt)
logarithmic_sequence <- rep(exp(seq(log(1e-3), log(50e+0), length.out = 100)), 19)

# plotly_response_grid tests
test_that("check output type, data, values", {
  p <- plotly_response_grid(prepared_curves, var_y)
  expect_type(p, "list")
  expect_equal(p$x$attrs[[1]]$x, logarithmic_sequence)
  expect_identical(p$x$attrs[[1]]$type, "scatter")
})

test_that("expect error with wrong input", {
  expect_error(plotly_response_grid(),
               "argument \"data\" is missing, with no default")
  expect_error(
    plotly_response_grid(n, var_y),
    "Assertion on 'data' failed: Must be a data.table, not double")
  expect_error(
    plotly_response_grid(prepared_curves, n),
    "Assertion on 'var_y' failed: Must be of type 'string', not 'double'")
  expect_error(
    plotly_response_grid(prepared_curves, var_y, range_x = n),
    "Assertion on 'range_x' failed: Must have length 2, but has length 1.")
  expect_error(
    plotly_response_grid(prepared_curves, var_y, title = n),
    "Assertion on 'title' failed: Must be of type 'string', not 'double'")
})
