# Prep table with calculated linear associations

Prep table with calculated linear associations

## Usage

``` r
prep_dt_assoc(dt_response, dt_depmap, selected_feat_meta_col = NULL)
```

## Arguments

- dt_response:

  `data.table` with experimental response data (rows are samples) for
  one metric

- dt_depmap:

  `data.table` with dependent variables data load from DepMap. (rows are
  samples, columns are features or meta); outputted by one of
  [`prep_dt_depmap_feat`](prep_dt_depmap_feat.md) or
  [`prep_dt_depmap_meta`](prep_dt_depmap_meta.md)

- selected_feat_meta_col:

  string name of feature/meta column in DepMap

## Value

A named list with elements, that may be input to
[`plot_volcano_assoc`](plot_volcano_assoc.md)

- `dt_assoc` `data.table` with calculated association values between
  feature/meta of DepMap and selected metric,

- `condition_info` string describing experiment condition (drugs),

- `selected_feat_meta_col` string name of feature/meta.

## Author

Janina Smoła <janina.smola@contractors.roche.com>
