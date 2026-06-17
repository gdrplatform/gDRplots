# Prep table with metric values for single-agent experiment

Prep table with metric values for single-agent experiment

## Usage

``` r
prep_dt_response_metric_sa(
  dt_metrics,
  d_name,
  normalization_type = "RV",
  metric = "xc50",
  fit_source = "gDR"
)
```

## Arguments

- dt_metrics:

  `data.table` representing data from the `Metrics` assay, outputted by
  [`gDRutils::convert_se_assay_to_dt`](https://gdrplatform.github.io/gDRstyle/reference/convert_se_assay_to_dt.html)
  and single-agent `SummarizedExperiment`

- d_name:

  string with drug name to be plotted (identifiers `DrugName`)

- normalization_type:

  string with normalization types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- metric:

  string name of the metric; one of: "xc50" ("GR50" or "IC50" -
  respectively depending on `normalization_type`), "x_max" ("GR Max" or
  "E Max") or "x_mean" ("GR Mean" or "RV Mean")

- fit_source:

  string source name for metrics

## Value

`data.table` with selected metric, input to
[`prep_dt_assoc`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_assoc.md)

## Examples

``` r
mae <- gDRutils::get_synthetic_data("combo_matrix_small")
se <- mae[[gDRutils::get_supported_experiments("sa")]]
dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
                                               assay_name = "Metrics")
d_name <- "drug_004"
dt_response <- prep_dt_response_metric_sa(dt_metrics, d_name)
dt_response <-
  prep_dt_response_metric_sa(dt_metrics, d_name,
                             metric = c("xc50", "x_mean", "x_max"))

dt_averaged <- gDRutils::convert_se_assay_to_dt(se = se,
                                                assay_name = "Averaged")
dt_metrics_capped <-
  gDRutils::cap_assay_infinities(conc_assay_dt = dt_averaged,
                                 assay_dt = dt_metrics,
                                 experiment_name = gDRutils::get_supported_experiments("sa"),
                                 capping_fold = 5)
d_name <- "drug_026"
dt_response <- prep_dt_response_metric_sa(dt_metrics, d_name)
dt_response <-
  prep_dt_response_metric_sa(dt_metrics, d_name,
                             metric = c("xc50", "x_mean", "x_max"))
```
