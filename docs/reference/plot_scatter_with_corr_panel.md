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
  [`prep_dt_response_metric_sa`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_response_metric_sa.md),
  [`prep_dt_response_dose_sa`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_response_dose_sa.md),
  [`prep_dt_response_scores`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_response_scores.md)
  or
  [`prep_dt_response_metric_diff`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_response_metric_diff.md),

- dt_depmap:

  `data.table` with dependent variables data loaded from DepMap - for
  one feature; (rows are samples, columns are features). outputted by
  [`prep_dt_depmap_feat`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_depmap_feat.md)

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

## Examples

``` r
mae_prism <- gDRutils::get_synthetic_data("prism")
se_prism <- mae_prism[[gDRutils::get_supported_experiments("sa")]]
dt_metrics <- gDRutils::convert_se_assay_to_dt(se_prism, "Metrics")
dt_response <- prep_dt_response_metric_sa(
  dt_metrics = dt_metrics,
  d_name = unique(dt_metrics[["DrugName"]])[1],
  normalization_type = "RV", metric = "x_mean")
feat_data_path <- system.file("depmap_data", package = "gDRtestData")
meta_data_path <- system.file("depmap_data/Model.csv.gz", package = "gDRtestData")
obj_depmap_feat <- prep_dt_depmap_feat(feat_data_path, meta_data_path,
                                       feature_set = "OmicsCNGene")
plot_scatter_with_corr_panel(
  dt_response = dt_response,
  dt_depmap = obj_depmap_feat[["dt_depmap"]],
  selected_feats = names(obj_depmap_feat[["dt_depmap"]])[3:5])

```
