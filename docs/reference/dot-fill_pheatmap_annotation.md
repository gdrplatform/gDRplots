# Prep annotation data.table acc to metric matrix for pheatmap::pheatmat

Prep annotation data.table acc to metric matrix for pheatmap::pheatmat

## Usage

``` r
.fill_pheatmap_annotation(
  dt_anno,
  mat_with_metric,
  anno_var = gDRutils::get_env_identifiers("cellline_name")
)
```

## Arguments

- dt_anno:

  `data.table` that specifies the annotations shown on left side of the
  heatmap or shown above the heatmap - depending on the `anno_var`. Each
  row defines the features for a specific row. The rows in the data and
  in the annotation are matched using corresponding names from the
  required `anno_var` column.

- mat_with_metric:

  numeric matrix with metric values; must have named rows and columns

- anno_var:

  string with variable describing annotation dimension: one of:
  `CellLineName` for rows or `DrugName` for column.

## Value

`data.table` with annotation updated to `mat_with_metric`

## Author

Janina Smoła <janina.smola@contractors.roche.com>
