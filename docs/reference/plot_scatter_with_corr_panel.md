# Plot panel with scatter with correlation

Plot panel with scatter with correlation

## Usage

``` r
plot_scatter_with_corr_panel(
  dt_response,
  dt_depmap,
  selected_feats,
  selected_feat_meta_col = NULL,
  ncol = NULL
)
```

## Arguments

- dt_response:

  `data.table` with experimental response data (rows are samples) for
  one metric outputted by one of functions:
  [`prep_dt_response_metric_sa`](prep_dt_response_metric_sa.md),
  [`prep_dt_response_dose_sa`](prep_dt_response_dose_sa.md),
  [`prep_dt_response_scores`](prep_dt_response_scores.md) or
  [`prep_dt_response_metric_diff`](prep_dt_response_metric_diff.md),

- dt_depmap:

  `data.table` with dependent variables data loaded from DepMap - for
  one feature; (rows are samples, columns are features). outputted by
  [`prep_dt_depmap_feat`](prep_dt_depmap_feat.md)

- selected_feats:

  character vector with names of selected features from `dt_depmap`

- selected_feat_meta_col:

  string with the name of a feature column in DepMap (will be used as a
  plot title)

- ncol:

  number of plot column in panel

## Value

`ggplot` object containing panel of scatter plot with correlation for
selected features

## Author

Janina Smoła <janina.smola@contractors.roche.com>
