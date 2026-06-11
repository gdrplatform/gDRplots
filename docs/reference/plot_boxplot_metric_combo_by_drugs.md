# Plot box plots for metric for combo data grouped by drug names

Plot box plots for metric for combo data grouped by drug names

## Usage

``` r
plot_boxplot_metric_combo_by_drugs(
  dt_scores,
  normalization_type = "GR",
  metric = "hsa_score",
  fit_source = "gDR",
  colored_pts_flag = FALSE,
  colors_vec = NULL
)
```

## Arguments

- dt_scores:

  `data.table` representing data from the `scores` assay, outputted by
  `gDRutils::convert_se_assay_to_dt(se, "scores")` and combo
  `SummarizedExperiment`

- normalization_type:

  string with normalization types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- metric:

  string name of the combo metric; one of: "hsa_score"("Bliss Excess GR"
  or "Bliss Excess RV" - respectively depending on `normalization_type`)
  or "bliss_score" ("Bliss Score GR" or "Bliss Score RV")

- fit_source:

  string source name for metrics

- colored_pts_flag:

  a logical flag whether the points should be colored by grouped
  variable - for `group_var` equal `"CellLineName"` points will be
  colored by `"DrugName"` and similarly vice versa

- colors_vec:

  character vector with colors (name or hex value) to color boxplots;
  for `grouped_flag` set as `FALSE` only first from vector will be used

## Value

`ggplot` object containing boxplots for selected combo metric grouped by
drug names

## Author

Janina Smoła <janina.smola@contractors.roche.com>

## Examples

``` r
mae <- gDRutils::get_synthetic_data("combo_matrix")
se <- mae[[gDRutils::get_supported_experiments("combo")]]

dt_scores <- gDRutils::convert_se_assay_to_dt(se = se,
                                              assay_name = "scores")

plot_boxplot_metric_combo_by_drugs(dt_scores)


plot_boxplot_metric_combo_by_drugs(dt_scores,
                                   normalization_type = "RV",
                                   colors_vec = "gold")


plot_boxplot_metric_combo_by_drugs(dt_scores,
                                   metric = "bliss_score",
                                   normalization_type = "RV",
                                   colored_pts_flag = TRUE)

```
