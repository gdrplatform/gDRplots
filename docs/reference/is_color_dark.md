# Determine whether or not a color is dark

Determine whether or not a color is dark

## Usage

``` r
is_color_dark(col_name)
```

## Arguments

- col_name:

  string representing a valid color

## Value

logical flag

## Examples

``` r
is_color_dark("blue")
#> [1] TRUE
is_color_dark("red")
#> [1] FALSE
is_color_dark("#000000")
#> [1] TRUE
```
