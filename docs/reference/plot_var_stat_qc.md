# Lollipop plot for metric single-agent data to control quality of the data

Lollipop plot for metric single-agent data to control quality of the
data

## Usage

``` r
plot_var_stat_qc(
  dt_assay,
  cl_name,
  metric = "r2",
  normalization_type = "GR",
  with_table = FALSE
)
```

## Arguments

- dt_assay:

  data.table representing data from the `Metrics` assay, outputted by
  `gDRutils::convert_se_assay_to_dt(se, "Metrics")` and single-agent
  `SummarizedExperiment`

- cl_name:

  string cell line name to be plotted (Cell Line Name)

- metric:

  string with variable name to be plotted; it has to be in `dt_assay`

- normalization_type:

  string with normalization_types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- with_table:

  logical whether table with metric values should be shown next to the
  plot

## Value

`ggplot` object containing lollipop plot with stat value for each drug

## Author

Bartosz Czech <czech.bartosz@external.gene.com>

## Examples

``` r
mae <- gDRutils::get_synthetic_data("small")
se <- mae[[gDRutils::get_supported_experiments("sa")]]

dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
cl_name <- dt_metrics[["CellLineName"]][1]

plot_var_stat_qc(dt_assay = dt_metrics,
                 cl_name = cl_name)


plot_var_stat_qc(dt_assay = dt_metrics,
                 cl_name = cl_name,
                 metric = "r2",
                 normalization_type = "RV")


plot_var_stat_qc(dt_assay = dt_metrics,
                 cl_name = cl_name,
                 metric = "x_AOC",
                 normalization_type = "RV",
                 with_table = TRUE)

```
