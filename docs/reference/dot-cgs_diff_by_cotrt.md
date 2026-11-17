# Reshape co-treatment differences to one row per drug and co-treatment concentration

`prep_dt_response_metric_diff` returns one difference column per metric
and co-treatment concentration (e.g.
`RV_gDR_x_mean_cotrt_diff_0.03_drug_1`). A cell line screened at one
concentration has `NA` in the columns of the others, so the columns are
stacked into one column per metric, with the concentration in column
`cotrt`. Rows without any value at a concentration are dropped.

## Usage

``` r
.cgs_diff_by_cotrt(metrics_diff)
```

## Arguments

- metrics_diff:

  data.table returned by `prep_dt_response_metric_diff`.

## Value

data.table with the identifier columns, `cotrt` and one column per
metric.
