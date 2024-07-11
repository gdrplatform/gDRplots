context("Test metricRanking plot")

test_that("check output type and data",  {
  testdata <- gDRutils::get_testdata()
  drug_names <- testdata$drug_names
  cell_line_names <- testdata$cell_line_names
  dt <- testdata$dt
  
  axis_limits <- "fixed"
  show_tick_labels <- TRUE
  title <- "Cell Line Name"
  var_col <- "Drug MOA"
  var_grp <- "Tissue"
  var_x <- "Drug Name"
  var_y <- "x_AOC"
  
  levels <- c(
    "drug_010",
    "drug_002",
    "drug_006",
    "drug_009",
    "drug_003",
    "drug_008",
    "drug_005",
    "drug_004",
    "drug_007",
    "drug_011"
  )
  
  plt <- plotly_metric_ranking(dt, var_x, var_y, var_col, var_grp, title, show_tick_labels)
  expect_type(plt, "list")
  expect_identical(plt$x$attrs[[1]]$type, "bar")
  expect_identical(sort(plt$x$attrs[[1]]$y), sort(dt$x_AOC))
  # skipping test for now because can't be reproduced locally
  # expect_identical(sort(p$x$attrs[[1]]$x),
  #                  sort(factor(dt$`Drug Name`, levels = levels)))
  
})

test_that("check output type and data",  {
  m <- 9
  n <- 5
  
  synthetic_data <- gDRutils::gen_synthetic_data(m, n)
  drug_names <- synthetic_data$drug_names
  cell_names <- synthetic_data$cell_names
  dt <- synthetic_data$dt
  
  axis_limits <- "fixed"
  show_tick_labels <- TRUE
  title <- "Cell Line Name"
  var_col <- "none"
  var_grp <- "none"
  var_x <- "Drug Name"
  var_y <- "GR_AOC"
  
  plt <- plotly_metric_ranking(dt, var_x, var_y, var_col, var_grp, title, show_tick_labels)
  expect_type(plt, "list")
  expect_identical(plt$x$attrs[[1]]$type, "bar")
  expect_identical(plt$x$attrs[[1]]$y, dt$GR_AOC)
  expect_identical(plt$x$attrs[[1]]$x,
                   factor(rep(drug_names, n), levels = drug_names))

  
  expect_error(
    plotly_metric_ranking(n, var_x, var_y, var_col, var_grp, title, show_tick_labels),
    "Assertion on 'data' failed: Must be a data.table, not double."
  )
  expect_error(
    plotly_metric_ranking(dt, n, var_y, var_col, var_grp, title, show_tick_labels),
    "Assertion on 'var_x' failed: Must be of type 'string', not 'double'."
  )
  expect_error(
    plotly_metric_ranking(dt, var_x, n, var_col, var_grp, title, show_tick_labels),
    "Assertion on 'var_y' failed: Must be of type 'string', not 'double'."
  )
  expect_error(
    plotly_metric_ranking(dt, var_x, var_y, n, var_grp, title, show_tick_labels),
    "Assertion on 'var_col' failed: Must be of type 'string', not 'double'."
  )
  expect_error(
    plotly_metric_ranking(dt, var_x, var_y, var_col, n, title, show_tick_labels),
    "Assertion on 'var_grp' failed: Must be of type 'string', not 'double'."
  )
  expect_error(
    plotly_metric_ranking(dt, var_x, var_y, var_col, var_grp, n, show_tick_labels),
    "Assertion on 'title' failed: Must be of type 'string', not 'double'."
  )
})
