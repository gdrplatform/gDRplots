# Prep matrix with metric value based on the Metrics assay

Prep matrix with metric value based on the Metrics assay

## Usage

``` r
prep_pheatmap_matrix(
  dt_response,
  normalization_type = "GR",
  metric = "xc50",
  fit_source = "gDR",
  experiment_type = gDRutils::get_supported_experiments("sa")
)
```

## Arguments

- dt_response:

  `data.table` representing data from the `Metrics` assay, outputted by
  [`gDRutils::convert_se_assay_to_dt`](https://gdrplatform.github.io/gDRstyle/reference/convert_se_assay_to_dt.html)
  and single-agent `SummarizedExperiment` or `data.table` representing
  data from the `scores` assay, outputted by
  [`gDRutils::convert_se_assay_to_dt`](https://gdrplatform.github.io/gDRstyle/reference/convert_se_assay_to_dt.html)
  and combo `SummarizedExperiment`

- normalization_type:

  string with normalization types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- metric:

  string name of the metric; one of: "xc50" ("GR50" or "IC50" -
  respectively depending on `normalization_type`), "x_max" ("GR Max" or
  "E Max") or "x_mean" ("GR Mean" or "RV Mean"); but the values from any
  numeric column can be displayed.

- fit_source:

  string source name for metrics

- experiment_type:

  string with experiment name; one of: "single-agent", "combination" or
  "co-dilution"

## Value

matrix with values for selected metric with `CellLinName` in the rows
and `DrugName` (or combination of `DrugName` and `DrugName_2`) in the
columns

## Author

Janina Smoła <janina.smola@contractors.roche.com>

## Examples

``` r
mae <- gDRutils::get_synthetic_data("combo_matrix")
se <- mae[[gDRutils::get_supported_experiments("sa")]]
dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
                                               assay_name = "Metrics")

mat <- prep_pheatmap_matrix(dt_response = dt_metrics,
                            experiment_type = gDRutils::get_supported_experiments("sa"))

mae <- gDRutils::get_synthetic_data("combo_matrix")
se <- mae[[gDRutils::get_supported_experiments("combo")]]
dt_scores <- gDRutils::convert_se_assay_to_dt(se = se,
                                              assay_name = "scores")

mat <- prep_pheatmap_matrix(dt_response = dt_scores,
                            metric = "hsa_score",
                            normalization_type = "RV",
                            experiment_type = gDRutils::get_supported_experiments("combo"))
```
