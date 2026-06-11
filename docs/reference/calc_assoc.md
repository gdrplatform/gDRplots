# Calculate linear associations

Calculate the linear model associations between dependent variables and
response variable(s) of interest.

## Usage

``` r
calc_assoc(X, Y)
```

## Arguments

- X:

  `matrix` dependent variables data matrix (rows are samples, columns
  are features). Must have the same number of rows as matrix `Y` or
  equal to length of vector `Y`

- Y:

  `vector` or `matrix` experimental response data (rows are samples).
  When `Y` is a matrix must have the same number of rows as matrix `X`;
  when `y` is a vector - its length has to be equal to number of rows in
  matrix `X`.

## Value

`data.table` with calculated linear associations

## Note

inspired by the `calc_assoc` function written by James Hawley

## See also

[`cdsrmodels::lin_associations`](https://rdrr.io/pkg/cdsrmodels/man/lin_associations.html)

## Examples

``` r
X <- matrix(rep(1:13, length.out = 42), nrow = 6,
            dimnames = list(sprintf("row_%s", 1:6), sprintf("feat_%s", 1:7)))
Y <- matrix(c(10:15, 110:115, 210:215), ncol = 3,
            dimnames = list(sprintf("row_%s", 1:6), sprintf("met_%s", 1:3)))
tab_assoc <- calc_assoc(X, Y)
```
