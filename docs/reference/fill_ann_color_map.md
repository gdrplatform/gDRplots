# Fill missing values in the color map for annotation

Fill missing values in the color map for annotation

## Usage

``` r
fill_ann_color_map(dt_ann, map_ann)
```

## Arguments

- dt_ann:

  `data.table` with the annotations

- map_ann:

  `list` with the annotations

## Value

list with color mapping for the annotations with missing items filled in

## See also

[`pheatmap_with_anno_sa`](https://gdrplatform.github.io/gDRplots/reference/pheatmap_with_anno_sa.md)
[`pheatmap_with_anno_combo`](https://gdrplatform.github.io/gDRplots/reference/pheatmap_with_anno_combo.md)

## Author

Janina Smoła <janina.smola@contractors.roche.com>

## Examples

``` r
annotation_manual <- data.table::data.table(
  CellLineName = c("cellline_AA", "cellline_EA", "cellline_IB", "cellline_MC", "cellline_BC"),
  mut_A = c(0, 0, 1, 2, 3),
  mut_B = c("yes", "yes", "no", NA, NA),
  mut_C = c("AA", "AA", "AB", NA, "B")
)

annotation_map <- list(
  mut_A = c("1" = "coral", "0" = "cadetblue"),
  mut_B = c("yes" = "black", "no" = "grey90", "not_checked" = "lightblue"),
  mut_C = c("AA" = "red", "AB" = "orange", "B" = "yellow")
)

fill_ann_color_map(dt_ann = annotation_manual, map_ann = annotation_map)
#> $mut_B
#>       yes        no        NA 
#>   "black"  "grey90" "darkred" 
#> 
#> $mut_C
#>       AA       AB        B       NA 
#>    "red" "orange" "yellow"  "black" 
#> 
```
