# Prep table with metric values for combination experiment

Prep table with metric values for combination experiment

## Usage

``` r
prep_dt_response_scores(
  dt_scores,
  d_name,
  d_name2,
  normalization_type = "RV",
  metric = "hsa_score",
  fit_source = "gDR"
)
```

## Arguments

- dt_scores:

  `data.table` representing data from the `scores` assay, outputted by
  [`gDRutils::convert_se_assay_to_dt`](https://gdrplatform.github.io/gDRstyle/reference/convert_se_assay_to_dt.html)
  and combo `SummarizedExperiment`

- d_name:

  string with drug name to be plotted (identifiers `DrugName`)

- d_name2:

  string with drug name to be plotted (identifiers `DrugName_2`)

- normalization_type:

  string with normalization types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- metric:

  string name of the combo metric; one of: "hsa_score"("Bliss Excess GR"
  or "Bliss Excess RV" - respectively depending on
  `normalization_type`), "bliss_score" ("Bliss Score GR" or "Bliss Score
  RV")

- fit_source:

  string source name for metrics

## Value

`data.table` with selected metric, input to
[`prep_dt_assoc`](prep_dt_assoc.md)

## Examples

``` r
mae <- gDRutils::get_synthetic_data("combo_matrix_small")
se <- mae[[gDRutils::get_supported_experiments("combo")]]
dt_scores <- gDRutils::convert_se_assay_to_dt(se = se,
                                              assay_name = "scores")
d_name <- "drug_004"
d_name2 <- "drug_026"
dt_response <- prep_dt_response_scores(dt_scores, d_name, d_name2)
dt_response <-
  prep_dt_response_scores(dt_scores, d_name, d_name2,
                          metric = c("hsa_score", "bliss_score"))
```
