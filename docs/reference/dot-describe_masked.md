# Describe the masked wells removed when the SE was created

Describe the masked wells removed when the SE was created

## Usage

``` r
.describe_masked(masked)
```

## Arguments

- masked:

  list with `n_wells` and, when known, `n_pairs` (drug x cell line pairs
  affected) and `n_flag_missing` (masked wells without a masked flag),
  or `NULL`.

## Value

A string, or `NULL` when `masked` is `NULL`.
