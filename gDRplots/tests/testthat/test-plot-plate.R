context("Test plot-plate")


# Create test data
test_data <- withr::with_seed(
  123,
  data.table::data.table(
    WellColumn = rep(1:12, each = 8),
    WellRow = rep(LETTERS[1:8], times = 12),
    Gnumber = c(rep("untreated", 48), sample(1:5, size = 48, replace = TRUE)),
    Gnumber_2 = c(rep("untreated", 48), sample(1:5, size = 48, replace = TRUE)),
    Concentration = runif(96, min = 0, max = 100),
    ReadoutValue = runif(96, min = 0, max = 100),
    clid = "CellLineA",
    Barcode = rep(c("A", "B"), 48)
  )
)

test_that("plot_plate_stack_info works correctly", {
  # Test if the function returns a list of ggplot objects
  plots <- plot_plate_stack_info(test_data)
  expect_length(plots, 2)
  expect_true(all(vapply(plots, inherits, what = "ggplot", FUN.VALUE = logical(1))))
  
  # Test if it handles data with a single barcode correctly
  single_barcode_data <- test_data[Barcode == "A"]
  single_plots <- plot_plate_stack_info(single_barcode_data)
  expect_length(single_plots, 1)
  expect_true(all(vapply(single_plots, inherits, what = "ggplot", FUN.VALUE = logical(1))))
  
  # Test lack of barcode
  nobarcode_data <- data.table::copy(test_data)[, Barcode := NULL]
  empty_plot <- plot_plate_stack_info(nobarcode_data)
  expect_length(ggplot2::ggplot_build(empty_plot)$data[[1]], 0)
  
  # Test for combo data
  combo_data <- data.table::copy(test_data)
  combo_data[, Concentration_2 := rep(withr::with_seed(123, runif(48, min = 0, max = 10)), 2)]
  combo_plot <- plot_plate_stack_info(combo_data)
  expect_length(combo_plot, NROW(unique(combo_data$Barcode)))
  expect_equal(NROW(unique(ggplot2::ggplot_build(combo_plot[[1]])$data[[1]]$fill)),
               NROW(unique(combo_data[Barcode == "A", ][["Concentration"]])))
  expect_equal(NROW(unique(ggplot2::ggplot_build(combo_plot[[1]])$data[[3]]$fill)),
               NROW(unique(combo_data[Barcode == "A", ][["Concentration_2"]])))
  
  # Test if it handles empty data
  empty_data <- test_data[0]
  empty_plots <- plot_plate_stack_info(empty_data)
  expect_length(empty_plots, 0)
  
  # Test if it handles missing required columns
  incomplete_data <- test_data[, .(WellColumn, WellRow, Concentration, Gnumber, Barcode)]
  expect_error(plot_plate_stack_info(incomplete_data), "but is missing elements")
})

test_that("plot_plate works correctly", {
  # Test if the function returns a list of ggplot objects
  plots <- plot_plate(test_data, "Gnumber")
  expect_length(plots, 2)
  expect_true(all(vapply(plots, inherits, what = "ggplot", FUN.VALUE = logical(1))))
  
  # Test if the function handles data with a single barcode correctly
  single_barcode_data <- test_data[Barcode == "A"]
  single_plots <- plot_plate(single_barcode_data, "Gnumber")
  expect_length(single_plots, 1)
  expect_true(all(vapply(single_plots, inherits, what = "ggplot", FUN.VALUE = logical(1))))
  
  # Test if the function handle numeric columns (rounding for Concentration)
  numeric_plot_1 <- plot_plate(single_barcode_data, "Concentration")
  heatmap_lbls <- ggplot2::ggplot_build(numeric_plot_1[[1]])$data[[2]]$label
  expect_equal(max(nchar(gsub("(\\d)+\\.", "", heatmap_lbls))), 5)
  expect_true(is.factor(heatmap_lbls))
  
  # Test if the function handle numeric columns (rounding for other than Concentration)
  numeric_plot_2 <- plot_plate(single_barcode_data, "ReadoutValue")
  heatmap_lbls <- ggplot2::ggplot_build(numeric_plot_2[[1]])$data[[2]]$label
  expect_equal(max(nchar(gsub("(\\d)+\\.", "", heatmap_lbls))), 5)
  expect_false(is.factor(heatmap_lbls))
  
  # Test lack of barcode
  nobarcode_data <- data.table::copy(test_data)[, Barcode := NULL]
  empty_plot <- plot_plate(nobarcode_data, column_name = "Gnumber")
  expect_length(ggplot2::ggplot_build(empty_plot)$data[[1]], 0)
  
  # Test if it handles empty data
  empty_data <- test_data[0]
  empty_plots <- plot_plate(empty_data, "Gnumber")
  expect_length(empty_plots, 0)
  
  # Test if it handles missing required columns
  incomplete_data <- test_data[, .(WellColumn, WellRow, Concentration, Barcode)]
  expect_error(plot_plate(incomplete_data, "Gnumber"), "Assertion")
})

test_that("generate_color_mappings works correctly", {
  # Test if the function generates color mappings correctly
  untrt_tag <- "untreated"
  color_mapping <- generate_color_mappings(test_data, untrt_tag)
  expect_true("untreated" %in% names(color_mapping))
  expect_true(length(color_mapping) > 1)
  
  # Test if the function handles no untreated tags correctly
  test_data_no_untrt <- test_data[Gnumber != "untreated"]
  color_mapping_no_untrt <- generate_color_mappings(test_data_no_untrt, untrt_tag)
  expect_true("untreated" %in% names(color_mapping_no_untrt))
})
