# Calculate the luminance of a color

Calculate the luminance of a color

## Usage

``` r
get_col_luminance(col_name)
```

## Arguments

- col_name:

  string representing a valid color

## Value

single element numeric vector

## Examples

``` r
get_col_luminance("blue")
#> [1] 0.0722
get_col_luminance("red")
#> [1] 0.2326
get_col_luminance("#000000")
#> [1] 0
get_col_luminance("#906090")
#> [1] 0.166325
get_col_luminance("#906090F2")
#> [1] 0.166325
```
