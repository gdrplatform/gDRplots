# Prep table with metric values by doses for single-agent experiment

Prep table with metric values by doses for single-agent experiment

## Usage

``` r
prep_dt_response_dose_sa(
  dt_average,
  d_name,
  normalization_type = "RV",
  metric = "x",
  fit_source = "gDR"
)
```

## Arguments

- dt_average:

  `data.table` representing data from the `Averaged` assay, outputted by
  [`gDRutils::convert_se_assay_to_dt`](https://gdrplatform.github.io/gDRstyle/reference/convert_se_assay_to_dt.html)
  and `SummarizedExperiment` with chosen data type: single-agent or
  combo

- d_name:

  string with drug name to be plotted (identifiers `DrugName`)

- normalization_type:

  string with normalization types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- metric:

  string name of the metric; one of: "x" (value of "GR" or "RV" itself -
  respectively depending on `normalization_type`), or "x_std" (standard
  deviation)

- fit_source:

  string source name for metrics

## Value

`data.table` with selected metric, input to
[`prep_dt_assoc`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_assoc.md)

## Examples

``` r
mae <- gDRutils::get_synthetic_data("combo_matrix_small")
se <- mae[[gDRutils::get_supported_experiments("sa")]]
dt_average <- gDRutils::convert_se_assay_to_dt(se = se,
                                               assay_name = "Averaged")
d_name <- "drug_004"
dt_response <- prep_dt_response_dose_sa(dt_average, d_name)
```
