# Get color palette for the isobologram levels

Get color palette for the isobologram levels

## Usage

``` r
.get_iso_colors(iso_levels)
```

## Arguments

- iso_levels:

  character vector with isobologram levels

## Value

gDR palette for isoline given in `iso_levels`

## Examples

``` r
# \donttest{
ls_iso_lvl <- c("0.25", "0.5", "0.75")
gDRplots:::.get_iso_colors(ls_iso_lvl)
#>      0.25       0.5      0.75 
#> "#F7AA18" "#ED6412" "#C52B06" 
# }
```
