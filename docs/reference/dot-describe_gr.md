# Describe how GR values were computed in normalization

Describe how GR values were computed in normalization

## Usage

``` r
.describe_gr(gr)
```

## Arguments

- gr:

  list with `n_pairs`, a named vector counting the drug x cell line
  pairs per GR method (`day0`, `ref_div_time`, `none`), and optionally
  `no_gr_cell_lines`; or `NULL`.

## Value

A string, or `NULL` when `gr` is `NULL`.
