# Run GSEA on one metric of one cell line (or cell line comparison)

Drugs without a finite value of the metric are left out, since `fgsea`
accepts only finite statistics. When no drug is left, an empty result is
returned.

## Usage

``` r
.cgs_fgsea(dt, metric, moa_list, name_col, moa_col, unit_name)
```

## Arguments

- dt:

  data.table with the metric column, drug names and MOAs.

- metric:

  string, the metric column.

- moa_list:

  named list of drug names per MOA.

- name_col, moa_col:

  strings, the drug name and MOA columns.

- unit_name:

  string, the cell line (or comparison) used in messages.

## Value

data.table with the `fgsea` result and a `median` column; `leadingEdge`
holds the leading-edge drugs as a comma-separated string.
