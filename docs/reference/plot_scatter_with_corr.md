# Plot scatter with correlation

Plot scatter with correlation

## Usage

``` r
plot_scatter_with_corr(
  dt_response,
  dt_depmap,
  selected_feat,
  selected_feat_meta_col = NULL
)
```

## Arguments

- dt_response:

  `data.table` with experimental response data (rows are samples) for
  one metric outputted by one of functions:
  [`prep_dt_response_metric_sa`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_response_metric_sa.md),
  [`prep_dt_response_dose_sa`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_response_dose_sa.md),
  [`prep_dt_response_scores`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_response_scores.md)
  or
  [`prep_dt_response_metric_diff`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_response_metric_diff.md),

- dt_depmap:

  `data.table` with dependent variables data loaded from DepMap - for
  one feature; (rows are samples, columns are features). outputted by
  [`prep_dt_depmap_feat`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_depmap_feat.md)

- selected_feat:

  string with name of selected feature from `dt_depmap`

- selected_feat_meta_col:

  string with the name of a feature column in DepMap (will be used as a
  plot title)

## Value

`ggplot` object containing a scatter plot with correlation

## Author

Janina Smoła <janina.smola@contractors.roche.com>
