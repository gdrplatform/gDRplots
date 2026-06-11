# Determine whether or not a color name is valid

A name of color is valid when either is a color name listed by
[`grDevices::colors`](https://rdrr.io/r/grDevices/colors.html) or a
hexadecimal string of the form `#rrggbb`

## Usage

``` r
is_valid_color(col_name)
```

## Arguments

- col_name:

  string representing a valid color

## Value

logical flag

## Examples

``` r
is_valid_color("darkblue")
#> [1] TRUE
is_valid_color("#FF8C00")
#> [1] TRUE
is_valid_color("#FF8C00DC")
#> [1] TRUE
is_valid_color("RED")
#> [1] FALSE
```
