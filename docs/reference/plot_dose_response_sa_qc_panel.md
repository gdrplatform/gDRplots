# Plot panel with drug response curves for single-agent data to control quality of the data

Plot panel with drug response curves for single-agent data to control
quality of the data

## Usage

``` r
plot_dose_response_sa_qc_panel(
  dt_metrics,
  dt_average,
  cl_name,
  d_names = NULL,
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

- d_names:

  character vector with drug names to be plotted (Drug Name); if NULL -
  all available drugs will be plotted

- normalization_type:

  string with normalization_types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- fit_source:

  string source name for metrics

## Value

`ggplot` object with panel with plots of dose-response curves for
selected cell line by drugs (observed and fitted values)

## Examples

``` r
mae <- gDRutils::get_synthetic_data("small")
se <- mae[[1]]

dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
cl_name <- dt_metrics[["CellLineName"]][1]
d_names <- unique(dt_metrics[["DrugName"]])[1:5]

plot_dose_response_sa_qc_panel(dt_metrics = dt_metrics,
                               dt_average = dt_average,
                               cl_name = cl_name)


plot_dose_response_sa_qc_panel(dt_metrics = dt_metrics,
                               dt_average = dt_average,
                               cl_name = cl_name,
                               d_names = d_names)

```
