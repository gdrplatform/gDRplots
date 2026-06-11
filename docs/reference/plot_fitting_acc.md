# Visualization for the quality control of the fitting for single-agent data

Visualization for the quality control of the fitting for single-agent
data

## Usage

``` r
plot_fitting_acc(dt_assay, cl_name, normalization_type = "GR")
```

## Arguments

- dt_assay:

  data.table representing data from the `Metrics` assay, outputted by
  `gDRutils::convert_se_assay_to_dt(se, "Metrics")` and single-agent
  `SummarizedExperiment`

- cl_name:

  string cell line name to be plotted (Cell Line Name)

- normalization_type:

  string with normalization_types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

## Value

`ggplot` object containing panel with lollipop plots with r2 and rss
values for each drug

## Author

Bartosz Czech <czech.bartosz@external.gene.com>

## Examples

``` r
mae <- gDRutils::get_synthetic_data("small")
se <- mae[[gDRutils::get_supported_experiments("sa")]]

dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
cl_name <- dt_metrics[["CellLineName"]][1]

plot_fitting_acc(dt_assay = dt_metrics,
                 cl_name = cl_name,
                 normalization_type = "RV")
```
