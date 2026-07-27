# Safely convert a character vector to numeric

Converts only values that match a numeric pattern to numeric; others
become NA.

## Usage

``` r
.safe_as_numeric(x)
```

## Arguments

- x:

  character vector to convert

## Value

numeric vector

## Examples

``` r
gDRplots:::.safe_as_numeric(c("0.01", "0.1", "untreated", NA))
#> [1] 0.01 0.10   NA   NA
```
