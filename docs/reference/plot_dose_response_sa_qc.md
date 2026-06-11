# Plot drug response curves for single-agent data to control quality of the data

Plot drug response curves for single-agent data to control quality of
the data

## Usage

``` r
plot_dose_response_sa_qc(
  dt_metrics,
  dt_average,
  cl_name,
  d_name,
  normalization_type = "GR",
  fit_source = "gDR"
)
```

## Arguments

- dt_metrics:

  data.table representing data from the `Metrics` assay, outputted by
  `gDRutils::convert_se_assay_to_dt(se, "Metrics")` and single-agent
  `SummarizedExperiment`

- dt_average:

  data.table representing data from the `Averaged` assay, outputted by
  `gDRutils::convert_se_assay_to_dt(se, "Averaged")` and single-agent
  `SummarizedExperiment`

- cl_name:

  string cell line name to be plotted (Cell Line Name)

- d_name:

  string vector with drug name to be plotted (Drug Name)

- normalization_type:

  string with normalization_types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- fit_source:

  string source name for metrics

## Value

`ggplot` object containing plot of dose-response curves (observed and
fitted values)

## Examples

``` r
mae <- gDRutils::get_synthetic_data("small")
se <- mae[[gDRutils::get_supported_experiments("sa")]]

dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
cl_name <- dt_metrics[["CellLineName"]][1]
d_name <- dt_metrics[["DrugName"]][1]

plot_dose_response_sa_qc(dt_metrics = dt_metrics,
                         dt_average = dt_average,
                         cl_name = cl_name,
                         d_name = d_name)

```
