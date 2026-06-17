# Prep table with metric values for combination experiment

Prep table with metric values for combination experiment

## Usage

``` r
prep_dt_response_metric_diff(
  dt_metrics,
  d_name,
  d_name2,
  cellline1 = NULL,
  cellline2 = NULL,
  normalization_type = "RV",
  metric = "xc50",
  fit_source = "gDR",
  additional_cols = NULL
)
```

## Arguments

- dt_metrics:

  `data.table` representing data from the `Metrics` assay, outputted by
  [`gDRutils::convert_se_assay_to_dt`](https://gdrplatform.github.io/gDRstyle/reference/convert_se_assay_to_dt.html)
  and combo `SummarizedExperiment`

- d_name:

  string representing the drug name to be plotted (identifier
  `DrugName`). If set to NULL, the function will return a table for all
  available DrugName

- d_name2:

  string representing the drug name to be plotted (identifier
  `DrugName_2`). If set to NULL, the function will return a table for
  all available DrugName_2

- cellline1:

  string representing the first cell line name. If set to NULL, the
  function will return a table for all available cell lines.

- cellline2:

  string representing the second cell line name. If set to NULL, the
  function will return a table for all available cell lines.

- normalization_type:

  string with normalization types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- metric:

  string name of the metric; one of: "xc50" ("GR50" or "IC50" -
  respectively depending on `normalization_type`), "x_max" ("GR Max" or
  "E Max"), "x_mean" ("GR Mean" or "RV Mean") "ec50" ("GEC50" or "EC50")
  or "x_AOC_range" ("GR AOC within set range" or "RV AOC within set
  range")

- fit_source:

  string source name for metrics

- additional_cols:

  character vector with additional cols that should be included in the
  output

## Value

`data.table` with selected metric, input to
[`prep_dt_assoc`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_assoc.md)

## Examples

``` r
mae <- gDRutils::get_synthetic_data("combo_matrix_small")
se <- mae[[gDRutils::get_supported_experiments("combo")]]
dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
                                               assay_name = "Metrics")
d_name <- "drug_004"
d_name2 <- "drug_026"
dt_response <-
  prep_dt_response_metric_diff(dt_metrics, d_name, d_name2,
  metric = c("xc50", "x_mean", "x_max"))

cellline1 <- "cellline_GB"
cellline2 <- "cellline_HB"

dt_response <-
  prep_dt_response_metric_diff(dt_metrics, d_name = NULL, d_name2 = NULL,
  cellline1, cellline2,
  metric = c("xc50", "x_mean", "x_max"))
```
