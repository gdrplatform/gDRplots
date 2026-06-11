# Get n-top linear associations

Currently, n-top linear associations are selected based on values of
`q_value` (column sorted increasing) and values of `rho` (column sorted
decreasing)

## Usage

``` r
.get_n_top_asssoc(dt_assoc, n_top = 4)
```

## Arguments

- dt_assoc:

  `data.table` with the calculated linear association between DepMap and
  metrics outputted by `calc_assoc`

- n_top:

  number of requested top linear associations

## Value

a vector with name of the n-top linear associations; when `n_top` will
be higher than number of available features - only available will be
returned.

## Examples

``` r
if (FALSE) { # \dontrun{
Y <- matrix(seq(0.5, 2, length.out = 50), nrow = 50,
            dimnames = list(sprintf("row_%s", 1:50), "met_1"))
X <- matrix(
  withr::with_seed(42, sample(c(NA, seq(0.35, 23.5, 1.25)),
                              size = 50*20, replace = TRUE)),
  nrow = 50,
  dimnames = list(sprintf("row_%s", 1:50), sprintf("feat_%s", 1:20)))
tab_assoc <- calc_assoc(X, Y)
.get_n_top_asssoc(tab_assoc)
} # }
```
