# Plot a single plate's stack info

Generates the QC and dose-rank heatmap for a single plate subset.
Handles all necessary data transformations (factor conversion, ranking,
color mapping) internally.

## Usage

``` r
plot_single_plate_stack_info(
  dt_subset,
  plate_id = NULL,
  drug_color_mapping = NULL,
  ctrl_fail_threshold = 0.6,
  n_sd = 1,
  use_sd_threshold = TRUE,
  items_per_line = 6
)
```

## Arguments

- dt_subset:

  `data.table`. Subset of the plate data corresponding to a single
  barcode/plate.

- plate_id:

  `character` (Optional). The identifier (barcode) of the plate. If
  `NULL`, attempts to extract it from the data.

- drug_color_mapping:

  `character` (Optional). Named vector of color codes. If `NULL`,
  generated automatically.

- ctrl_fail_threshold:

  `numeric` (Default: 0.6). Flags control wells with a readout below
  `Mean(Controls) * ctrl_fail_threshold`.

- n_sd:

  `numeric` (Default: 1). Sets the upper limit for flagging high signals
  using the formula `Mean(Controls) + (n_sd * SD(Controls))`.

- use_sd_threshold:

  `logical` (Default: TRUE). If `TRUE`, uses the dynamic SD-based limit
  defined by `n_sd`. If `FALSE`, uses a fixed limit of
  `Mean(Controls) * 1.1`.

- items_per_line:

  `integer`. Number of doses to show per line in the legend. Default 6.

## Value

A `ggplot` object.

## Author

Bartosz Czech <czech.bartosz@external.gene.com>

## Examples

``` r
conc_series <- c(0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10)

test_data <- data.table::data.table(
   WellColumn = rep(1:12, each = 8),
   WellRow = rep(LETTERS[1:8], times = 12),
   clid = "CellLineA",
   Barcode = "Plate_1"
)

invisible(test_data[, Concentration := rep(conc_series, length.out = .N)])
invisible(test_data[, Gnumber := ifelse(Concentration == 0, "vehicle", "Drug_A")])
invisible(test_data[Gnumber == "vehicle", Concentration := 0])

invisible(test_data[, ReadoutValue := ifelse(Gnumber == "vehicle",
                                    rnorm(.N, 1000, 50),
                                    rnorm(.N, 1000 * (1 / (1 + Concentration)), 50))])

library(ggtext)
plot_single_plate_stack_info(test_data)


invisible(test_data[, Concentration_2 := rep(rev(conc_series), length.out = .N)])
invisible(test_data[, Gnumber_2 := ifelse(Concentration_2 == 0, "vehicle", "Drug_B")])

invisible(test_data[Gnumber_2 == "vehicle", Concentration_2 := 0])

plot_single_plate_stack_info(test_data)

```
