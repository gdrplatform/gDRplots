# Create test data
set.seed(123)
test_data <- data.table::data.table(
  WellColumn = rep(1:12, each = 8),
  WellRow = rep(LETTERS[1:8], times = 12),
  Gnumber = c(rep("untreated", 48), sample(1:5, size = 48, replace = TRUE)),
  Gnumber_2 = c(rep("untreated", 48), sample(1:5, size = 48, replace = TRUE)),
  Concentration = runif(96, min = 0, max = 100),
  ReadoutValue = runif(96, min = 0, max = 100),
  Barcode = rep(c("A", "B"), 48)
)

test_that("plot_plate_data works correctly", {
  # Test if the function returns a list of ggplot objects
  plots <- plot_plate_data(test_data)
  expect_length(plots, 2)
  expect_true(all(sapply(plots, inherits, what = "ggplot")))
  
  # Test if it handles data with a single barcode correctly
  single_barcode_data <- test_data[Barcode == "A"]
  single_plots <- plot_plate_data(single_barcode_data)
  expect_length(single_plots, 1)
  expect_true(all(sapply(single_plots, inherits, what = "ggplot")))
  
  # Test if it handles empty data
  empty_data <- test_data[0]
  empty_plots <- plot_plate_data(empty_data)
  expect_length(empty_plots, 0)
  
  # Test if it handles missing required columns
  incomplete_data <- test_data[, .(WellColumn, WellRow, Concentration, Gnumber, Barcode)]
  expect_error(plot_plate_data(incomplete_data), "but is missing elements")
})

test_that("plot_plate_single_data works correctly", {
  # Test if the function returns a list of ggplot objects
  plots <- plot_plate_single_data(test_data, "Gnumber")
  expect_length(plots, 2)
  expect_true(all(sapply(plots, inherits, what = "ggplot")))
  
  # Test if it handles data with a single barcode correctly
  single_barcode_data <- test_data[Barcode == "A"]
  single_plots <- plot_plate_single_data(single_barcode_data, "Gnumber")
  expect_length(single_plots, 1)
  expect_true(all(sapply(single_plots, inherits, what = "ggplot")))
  
  # Test if it handles empty data
  empty_data <- test_data[0]
  empty_plots <- plot_plate_single_data(empty_data, "Gnumber")
  expect_length(empty_plots, 0)
  
  # Test if it handles missing required columns
  incomplete_data <- test_data[, .(WellColumn, WellRow, Concentration, Barcode)]
  expect_error(plot_plate_single_data(incomplete_data, "Gnumber"), "Assertion")
})

test_that("filter_data_by_barcode works correctly", {
  # Test if the function filters data correctly
  filtered_data <- filter_data_by_barcode(test_data, "A", "Barcode")
  expect_equal(nrow(filtered_data), 48)
  expect_equal(unique(filtered_data$Barcode), "A")
  
  # Test if the function handles no matching barcode
  no_match_data <- filter_data_by_barcode(test_data, "C", "Barcode")
  expect_equal(nrow(no_match_data), 0)
  
  # Test if the function handles incorrect barcode identifier
  expect_error(filter_data_by_barcode(test_data, "A", "non_existent_column"))
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
