context("Test fun plotly")

test_that("gDR_plotly_config works as expected", {
  
  plt <- plotly::plot_ly(x = rnorm(50), y = rnorm(50), type = "scatter", mode = "markers")
  
  plt_1 <- gDR_plotly_config(plt)
  expect_false(plt_1$x$config$displaylogo)
  expect_false(plt_1$x$config$showSendToCloud)
  expect_false(plt_1$x$config$showEditInChartStudio)
  expect_false(plt_1$x$config$showLink)
  expect_false(plt_1$x$config$sendData)
  expect_true(plt_1$x$config$responsive)
  expect_equal(unlist(plt_1$x$config$modeBarButtons),
               c("zoom2d", "resetViews", "autoScale2d",
                 "toggleSpikelines", "toggleHover", "hoverCompareCartesian"))
  expect_false(plt_1$x$config$editable)
  
  
  plt_2 <- gDR_plotly_config(plt, editable = TRUE)
  expect_true(plt_2$x$config$editable)
  
  plt_3 <- gDR_plotly_config(plt, edits = get_plotly_edits())
  expect_false(plt_3$x$config$editable)
  expect_equal(unlist(plt_3$x$config$edits),
               c(axisTitleText = TRUE, titleText = TRUE))
  
  expect_error(gDR_plotly_config(plt = "str"))
  expect_error(gDR_plotly_config(plt = plt, editable = "str"))
})