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
equal 0 (-Inf) is replaced with one step less in the dose dilution

## Examples

``` r
if (FALSE) { # \dontrun{
vec <- c(0, 0.003, 0.01, 0.03)
transform_log_conc(vec)
} # }
```
