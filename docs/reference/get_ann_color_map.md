# Create color map for annotation

Create color map for annotation

## Usage

``` r
get_ann_color_map(dt_ann)
```

## Arguments

- dt_ann:

  `data.table` with the annotations

## Value

list with color mapping for the annotations

## See also

[`pheatmap_qc`](https://gdrplatform.github.io/gDRplots/reference/pheatmap_qc.md)

## Author

Janina Smoła <janina.smola@contractors.roche.com>

## Examples

``` r
mae <- gDRutils::get_synthetic_data("small")
#> Loading required namespace: MultiAssayExperiment
se <- mae[[gDRutils::get_supported_experiments("sa")]][2:5, ]
#> Loading required namespace: BumpyMatrix
dt_average <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "Averaged")
dt_ann <- dt_average[,.SD, .SDcols = c("Tissue", "ReferenceDivisionTime")]

get_ann_color_map(dt_ann)
#> $Tissue
#>  tissue_x  tissue_y  tissue_z 
#> "#1B9E77" "#D95F02" "#7570B3" 
#> 
#> $ReferenceDivisionTime
#>        26        30        34        38        42        46        50        54 
#> "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666" "#A6CEE3" "#1F78B4" "#B2DF8A" 
#>        58        62 
#> "#33A02C" "#FB9A99" 
#> 
```
