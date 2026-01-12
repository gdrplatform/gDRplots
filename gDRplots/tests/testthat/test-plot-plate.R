context("Test plot-plate")


# --- 1. Create Robust Test Data ---
test_data <- withr::with_seed(
  123,
  {
    n <- 96
    dt <- data.table::data.table(
      WellColumn = sprintf("%02d", rep(1:12, each = 8)),
      WellRow = rep(LETTERS[1:8], times = 12),
      Barcode = rep(c("Plate_1", "Plate_2"), each = 48),
      clid = "CellLineA"
    )
    dt[, is_ctrl := rep(c(rep(TRUE, 6), rep(FALSE, 42)), 2)]
    
    dt[, Gnumber := ifelse(is_ctrl, "vehicle", "Drug_A")]
    dt[, Concentration := ifelse(is_ctrl, 0, runif(.N, 0.001, 10))]
    
    dt[, ReadoutValue := ifelse(is_ctrl, 
                                rnorm(.N, mean = 1000, sd = 50), 
                                rnorm(.N, mean = 500, sd = 100))]
    
    dt[, is_ctrl := NULL]
    dt
  }
)

test_that("plot_plate_stack_info works correctly", {

  plots <- plot_plate_stack_info(test_data)
  expect_length(plots, 2)
  expect_true(all(vapply(plots, inherits, what = "ggplot", FUN.VALUE = logical(1))))
  
  single_barcode_data <- test_data[Barcode == "Plate_1"]
  single_plots <- plot_plate_stack_info(single_barcode_data)
  expect_length(single_plots, 1)
  expect_true(inherits(single_plots[[1]], "ggplot"))
  
  nobarcode_data <- data.table::copy(test_data)[, Barcode := NULL]
  expect_error(plot_plate_stack_info(nobarcode_data), " Names must include the elements")
  
  combo_data <- data.table::copy(test_data)
  combo_data[, Gnumber_2 := "Drug_B"]
  combo_data[, Concentration_2 := ifelse(Gnumber == "vehicle", 0, runif(.N, 0.001, 5))]
  
  combo_plot <- plot_plate_stack_info(combo_data)
  
  expect_length(combo_plot, 2)
  p_built <- ggplot2::ggplot_build(combo_plot[[1]])

  expect_true("GeomPoint" %in% class(p_built$plot$layers[[2]]$geom))

  layers_classes <- vapply(p_built$plot$layers, function(x) class(x$geom)[1], character(1))
  point_layers <- which(layers_classes == "GeomPoint")
  expect_gte(length(point_layers), 2)
  
  empty_data <- test_data[0]
  empty_plots <- plot_plate_stack_info(empty_data)
  expect_length(empty_plots, 0)
  
  incomplete_data <- test_data[, .(WellColumn, WellRow, Concentration, Gnumber, Barcode)] 
  expect_error(plot_plate_stack_info(incomplete_data), "Names must include the elements")
  
  expect_error(plot_plate_stack_info(test_data, ctrl_fail_threshold = 1.5), "Element 1 is not <= 1")
  expect_error(plot_plate_stack_info(test_data, items_per_line = 0), "Element 1 is not >= 1")
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
