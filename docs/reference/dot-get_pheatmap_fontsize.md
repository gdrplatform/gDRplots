# Get fontsize for rownames or colnames in pheatmap::pheatmap

Get fontsize for rownames or colnames in pheatmap::pheatmap

## Usage

``` r
.get_pheatmap_fontsize(
  matrix,
  dimension = c("row", "col"),
  threshold_count = 40L
)
```

## Arguments

- matrix:

  numeric matrix with metric values.

- dimension:

  character value, either "row" or "col", indicating whether to compute
  fontsize for rows or columns.

- threshold_count:

  integer value of the number of rows/columns for which the font size
  remains standard.

## Value

numeric value of font size.

## Author

Janina Smoła <janina.smola@contractors.roche.com>
