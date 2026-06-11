# get_iso_colors

get_iso_colors

## Usage

``` r
get_iso_colors(normalization_type = c("RV", "GR"))
```

## Arguments

- normalization_type:

  charvec normalization_types expected in the data

## Value

named charvec with iso colors

## Examples

``` r
get_iso_colors()
#>         0      0.05       0.1      0.15       0.2      0.25       0.3      0.35 
#> "#46aaff" "#46a3f5" "#469ceb" "#4694e2" "#468dd8" "#4686ce" "#467fc4" "#4677ba" 
#>       0.4      0.45       0.5      0.55       0.6      0.65       0.7      0.75 
#> "#4670b0" "#4669a6" "#46629d" "#465b93" "#465389" "#464c7f" "#464575" "#463e6b" 
#>       0.8      0.85       0.9      0.95         1 
#> "#463661" "#462f58" "#46284e" "#462144" "#461a3a" 
```
