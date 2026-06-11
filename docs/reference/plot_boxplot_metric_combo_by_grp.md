# Plot box plots for metric for combo data grouped by selected variable

Plot box plots for metric for combo data grouped by selected variable

## Usage

``` r
plot_boxplot_metric_combo_by_grp(
  dt_scores,
  selection_var,
  selection_name,
  group_var,
  group_names = NULL,
  normalization_type = "GR",
  metric = "hsa_score",
  fit_source = "gDR",
  named_n = 5,
  named_n_mode = "top",
  grouped_flag = FALSE,
  colors_vec = NULL
)
```

## Arguments

- dt_scores:

  `data.table` representing data from the `scores` assay, outputted by
  `gDRutils::convert_se_assay_to_dt(se, "scores")` and combo
  `SummarizedExperiment`

- selection_var:

  string name of selected main variable - one value from column
  `"CellLineName"` or `"DrugName"`

- selection_name:

  string name of selected variable value from column `selection_var` to
  filter data for plotting

- group_var:

  string name of group variable; should to be not numeric variable from
  `dt_metrics` different than `selection_var` and not containing unique
  values for each row

- group_names:

  character vector with names to subset from column `group_var`; if
  `NULL` then all values will be plotted

- normalization_type:

  string with normalization types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- metric:

  string name of the combo metric; one of: "hsa_score"("Bliss Excess GR"
  or "Bliss Excess RV" - respectively depending on `normalization_type`)
  or "bliss_score" ("Bliss Score GR" or "Bliss Score RV")

- fit_source:

  string source name for metrics

- named_n:

  number of points to label based on the highest or lowest `metric`
  values; if `group_var` is `"DrugName"`, points are labeled by
  `"CellLineName"` and similarly vice versa

- named_n_mode:

  string determines whether the labels are applied to the highest or
  lowest values of `metric`; one of: `"top"` or `"bottom"`

- grouped_flag:

  logical flag whether the boxplots should be colored by `group_var`

- colors_vec:

  character vector with colors (name or hex value) to color boxplots;
  for `grouped_flag` set as `FALSE` only first from vector will be used

## Value

`ggplot` object containing boxplots for selected combo metric grouped by
selected variable

## Author

Janina Smoła <janina.smola@contractors.roche.com>

## Examples

``` r
mae <- gDRutils::get_synthetic_data("combo_matrix")
se <- mae[[gDRutils::get_supported_experiments("combo")]]

dt_scores <- gDRutils::convert_se_assay_to_dt(se = se,
                                              assay_name = "scores")
invisible(dt_scores[, Tissue_grp := data.table::fifelse(Tissue == "tissue_w",
                                                        "tissue_w",
                                                        "tissue_other")])

plot_boxplot_metric_combo_by_grp(dt_scores,
                                 selection_var = "DrugName",
                                 selection_name = c("drug_001", "drug_021"),
                                 group_var = "Tissue")


plot_boxplot_metric_combo_by_grp(dt_scores,
                                 selection_var = "DrugName",
                                 selection_name = c("drug_001", "drug_021"),
                                 group_var = "Tissue_grp",
                                 grouped_flag = TRUE,
                                 colors_vec = c("darkblue", "deeppink"))


plot_boxplot_metric_combo_by_grp(dt_scores,
                                 selection_var = "DrugName",
                                 selection_name = c("drug_001", "drug_021"),
                                 group_var = "Tissue",
                                 grouped_flag = TRUE)


plot_boxplot_metric_combo_by_grp(dt_scores,
                                 selection_var = "DrugName",
                                 selection_name = c("drug_001", "drug_021"),
                                 group_var = "Tissue",
                                 colors_vec = c("darkblue", "orange"))

```
