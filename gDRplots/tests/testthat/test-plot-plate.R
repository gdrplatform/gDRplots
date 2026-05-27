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

test_that("plot_plate_stack_info works as expected", {
  # Test basic list return
  plots <- plot_plate_stack_info(test_data)
  expect_length(plots, 2)
  expect_true(all(vapply(plots, inherits, what = "ggplot", FUN.VALUE = logical(1))))
  expect_named(plots, c("Plate_1", "Plate_2"))

  # Test single barcode input
  single_barcode_data <- test_data[Barcode == "Plate_1"]
  single_plots <- plot_plate_stack_info(single_barcode_data)
  expect_length(single_plots, 1)
  expect_true(inherits(single_plots[[1]], "ggplot"))

  # Test Combo Data Integration
  combo_data <- data.table::copy(test_data)
  combo_data[, Gnumber_2 := "Drug_B"]
  combo_data[, Concentration_2 := ifelse(Gnumber == "vehicle", 0, runif(.N, 0.001, 5))]

  combo_plot <- plot_plate_stack_info(combo_data)
  expect_length(combo_plot, 2)

  p_built <- ggplot2::ggplot_build(combo_plot[[1]])

  # Check for presence of point layers
  layers_classes <- vapply(p_built$plot$layers, function(x) class(x$geom)[1], character(1))
  point_layers <- which(layers_classes == "GeomPoint")
  # Expect at least 2 point layers (Drug 1 + Drug 2)
  expect_gte(length(point_layers), 2)

  # Test Parameter Validation
  expect_error(plot_plate_stack_info(test_data, ctrl_fail_threshold = 1.5), "Element 1 is not <= 1")
  expect_error(plot_plate_stack_info(test_data, items_per_line = 0), "Element 1 is not >= 1")

  # Test missing columns
  incomplete_data <- test_data[, .(WellColumn, WellRow, Concentration, Gnumber, Barcode)]
  expect_error(plot_plate_stack_info(incomplete_data), "Names must include the elements")
})

test_that("plot_single_plate_stack_info works as expected", {

  # Extract data for a single plate
  single_plate_data <- test_data[Barcode == "Plate_1"]

  # 1. Basic Functionality
  p <- plot_single_plate_stack_info(single_plate_data)
  expect_true(inherits(p, "ggplot"))

  p_built <- ggplot2::ggplot_build(p)
  layers_classes <- vapply(p_built$plot$layers, function(x) class(x$geom)[1], character(1))
  expect_true("GeomPoint" %in% layers_classes)
  expect_true("GeomTile" %in% layers_classes)

  # 2. Combo Data Support
  combo_data <- data.table::copy(single_plate_data)
  combo_data[, Gnumber_2 := "Drug_B"]
  combo_data[, Concentration_2 := ifelse(Gnumber == "vehicle", 0, runif(.N, 0.001, 5))]

  p_combo <- plot_single_plate_stack_info(combo_data)
  p_combo_built <- ggplot2::ggplot_build(p_combo)

  # Check for multiple point layers (Drug 1 and Drug 2)
  layers_classes_combo <- vapply(p_combo_built$plot$layers, function(x) class(x$geom)[1], character(1))
  expect_gte(sum(layers_classes_combo == "GeomPoint"), 2)

  # 3. Custom Parameters
  # Test with valid parameters (should not error)
  expect_silent(plot_single_plate_stack_info(single_plate_data, ctrl_fail_threshold = 0.8, n_sd = 2))

  # 4. Implicit Data Calculations (Rank & Color)
  # Remove rank and ensure it plots (function should calculate rank internally)
  dt_norank <- data.table::copy(single_plate_data)
  # rank_1 is calculated inside the function if missing, so we ensure input doesn't have it
  if ("rank_1" %in% names(dt_norank)) {
    dt_norank[, rank_1 := NULL]
  }

  p_calc <- plot_single_plate_stack_info(dt_norank)
  expect_true(inherits(p_calc, "ggplot"))

  # 5. Multiple Shapes (Cell Lines)
  multi_shape_data <- data.table::copy(single_plate_data)
  mid_point <- floor(NROW(multi_shape_data) / 2)
  multi_shape_data[seq_len(mid_point), clid := "CellLine_A"]
  multi_shape_data[(mid_point + 1):.N, clid := "CellLine_B"]

  p_shape <- plot_single_plate_stack_info(multi_shape_data)
  p_shape_built <- ggplot2::ggplot_build(p_shape)

  # Find the main point layer (usually the first GeomPoint)
  point_layer_idx <- which(vapply(p_shape_built$plot$layers, function(x) class(x$geom)[1],
                                  character(1)) == "GeomPoint")[1]
  point_data <- p_shape_built$data[[point_layer_idx]]

  # Check that multiple shapes are rendered
  expect_gt(length(unique(point_data$shape)), 1)
})

test_that("plot_plate works as expected", {
  # Test if the function returns a list of ggplot objects
  plots <- plot_plate(test_data, "Gnumber")
  expect_length(plots, 2)
  expect_true(all(vapply(plots, inherits, what = "ggplot", FUN.VALUE = logical(1))))

  # Test if the function handles data with a single barcode correctly
  single_barcode_data <- test_data[Barcode == "Plate_1"]
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
  # Depending on implementation, ggplot build data might be empty list or 0 rows
  p_built <- ggplot2::ggplot_build(empty_plot)
  expect_length(p_built$data[[1]], 0)

  # Test if it handles empty data
  empty_data <- test_data[0]
  empty_plots <- plot_plate(empty_data, "Gnumber")
  expect_length(empty_plots, 0)

  # Test if it handles missing required columns
  incomplete_data <- test_data[, .(WellColumn, WellRow, Concentration, Barcode)]
  expect_error(plot_plate(incomplete_data, "Gnumber"), "Assertion")
})

test_that("generate_color_mappings works as expected", {
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

test_that(".calc_dose_rank works as expected", {
  # Test basic ranking
  vals <- c(0.1, 0.5, 0.1, 1.0, 0.5)
  expect_equal(.calc_dose_rank(vals), c("1", "2", "1", "3", "2"))

  # Test with 0 and NA
  vals_ctrl <- c(0, 0.5, NA, 0.1)
  expect_equal(.calc_dose_rank(vals_ctrl), c("-", "2", "-", "1"))

  # Test rounding behavior
  vals_round <- c(0.1234567, 0.1234568)
  expect_equal(.calc_dose_rank(vals_round), c("1", "1"))

  # Test empty/all zero
  expect_equal(.calc_dose_rank(c(0, 0)), c("-", "-"))
})

test_that(".format_dose_list works as expected", {
  doses <- c(0.001, 0.01, 0.1, 0.5, 1.0, 10.0)

  # Test normal splitting (4 items per line)
  res_4 <- .format_dose_list(doses, n_per_line = 4)
  expect_match(res_4, "1: 0.001 | 2: 0.01", fixed = TRUE)
  expect_match(res_4, "<br>", fixed = TRUE)

  # Test single line (all fit)
  res_all <- .format_dose_list(doses, n_per_line = 10)
  expect_false(grepl("<br>", res_all))

  # Test strict splitting (1 item per line)
  res_1 <- .format_dose_list(c(0.1, 0.2), n_per_line = 1)
  expect_equal(res_1, "1: 0.1<br>2: 0.2")

  # Test empty input
  expect_equal(.format_dose_list(numeric(0), 4), "")
})
