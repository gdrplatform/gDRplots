# Extract drug grid and feature metadata from file names

Extract drug grid and feature metadata from file names

## Usage

``` r
.get_info_from_names(file_names, normalization_type = "RV")
```

## Arguments

- file_names:

  character vector of file names following the pattern
  `prefix__feat_drug_NORM_...`

- normalization_type:

  string normalization type ("GR" or "RV")

## Value

named list with `drug_grid` and `feat_meta`

## Examples

``` r
# \donttest{
gDRplots:::.get_info_from_names("plt__metricA_DrugX_RV_extra", "RV")
#> $drug_grid
#> [1] "DrugX"
#> 
#> $feat_meta
#> [1] "metricA"
#> 
# }
```
