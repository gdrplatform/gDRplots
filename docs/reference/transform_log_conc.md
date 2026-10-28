# Transform concentrations values with log10

Transform concentrations values with log10

## Usage

``` r
transform_log_conc(conc_vec)
```

## Arguments

- conc_vec:

  numeric vector with concentration values

## Value

numeric vector with log10 concentration values; log10 for concentration
equal 0 (-Inf) is replaced with one step less in the dose dilution. When
all concentrations are 0, all positions are 0.

## Examples

``` r
# \donttest{
vec <- c(0, 0.003, 0.01, 0.03)
gDRplots:::transform_log_conc(vec)
#> [1] -3.045757 -2.522879 -2.000000 -1.522879
vec <- c(0, 0, 0)
gDRplots:::transform_log_conc(vec)
#> [1] 0 0 0
# }
```
