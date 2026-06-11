# create a log-sequence

Create a sequence of numbers growing in log-domain.

## Usage

``` r
create_log_seq(start, end, length)
```

## Arguments

- start, end:

  numeric, lower and upper margins of the sequence

- length:

  integer, resulting sequence length

## Value

A numeric vector, see `Details`.

## Details

The result is a numeric vector of length `length`. Differences between
items are constant in logarithmic domain and therefore geometrically
increase in linear domain.

## Examples

``` r
create_log_seq(1, 2, 5)
#> [1] 1.000000 1.189207 1.414214 1.681793 2.000000
```
