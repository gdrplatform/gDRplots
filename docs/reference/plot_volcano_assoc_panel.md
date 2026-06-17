# Plot panel with volcano plot and according to the data type - scatter plots or box plots

Plot panel with volcano plot and according to the data type - scatter
plots or box plots

## Usage

``` r
plot_volcano_assoc_panel(
  dt_response,
  dt_depmap,
  selected_metric,
  selected_feat_meta_col
)
```

## Arguments

- dt_response:

  `data.table` with the experimental response data (rows are samples)
  for one metric outputted by one of functions:
  [`prep_dt_response_metric_sa`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_response_metric_sa.md),
  [`prep_dt_response_dose_sa`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_response_dose_sa.md),
  [`prep_dt_response_scores`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_response_scores.md)
  or
  [`prep_dt_response_metric_diff`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_response_metric_diff.md),
  must have at least a column with `CellLineName` and a numeric column
  with metric values.

- dt_depmap:

  `data.table` with dependent variables data loaded from DepMap where
  rows are samples, columns are features/metadata levels; one of: data
  for one feature outputted by
  [`prep_dt_depmap_feat`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_depmap_feat.md)
  or data or data for one metadata outputted by
  [`prep_dt_depmap_meta`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_depmap_meta.md)

- selected_metric:

  string name of the metric in `dt_response`

- selected_feat_meta_col:

  string with name of selected feature from `dt_depmap` or the name of
  the selected metadata from `dt_depmap` - respectively

## Value

A named list with elements:

- `assoc_data` table with association data

- `panel` `ggplot` object containing a panel with volcano plot and
  depending on data type: a scatter plots with correlation for top 4
  variables or boxplots for variable levels

## Author

Janina Smoła <janina.smola@contractors.roche.com>
