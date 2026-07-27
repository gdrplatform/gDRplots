# Get color palette for the dose response curves for combination data

Get color palette for the dose response curves for combination data

## Usage

``` r
.get_combo_curves_colors(ls_conc_2)
```

## Arguments

- ls_conc_2:

  factor vector with values for `Concentration_2`

## Value

gDR palette for Concentration_2 given in `ls_conc_2`

## Examples

``` r
# \donttest{
ls_conc <- factor(c("0.001", "0.01", "1"))
gDRplots:::.get_combo_curves_colors(ls_conc)
#>     0.001      0.01         1 
#> "#57042C" "#FB9C05" "#D7F434" 
# }
```
