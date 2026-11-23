# Compute pairwise linear associations between columns of X and Y

Adapted from `lin_associations` in the
[cdsrmodels](https://github.com/broadinstitute/cdsr_models) package (MIT
License; authors William Colgan and Andrew Boghossian) to remove a
GitHub-only package dependency. The license text and copyright notice
are in `inst/licenses/cdsrmodels-MIT.txt`
(`system.file("licenses", "cdsrmodels-MIT.txt", package = "gDRplots")`).
Uses [`WGCNA::cor`](https://rdrr.io/pkg/WGCNA/man/cor.html) for fast
pairwise correlation on large matrices and
[`ashr::ash`](https://rdrr.io/pkg/ashr/man/ash.html) for empirical-Bayes
shrinkage of effect sizes.

## Usage

``` r
.lin_associations(
  X,
  Y,
  n.min = 4L,
  shrinkage = TRUE,
  alpha = 0,
  MHC_direction = NULL
)
```

## Arguments

- X:

  `matrix` of independent variables (rows = samples, cols = features).

- Y:

  `matrix` or `vector` of response variables (rows = samples).

- n.min:

  integer; minimum number of degrees of freedom (finite paired
  observations minus 2) required to compute a p-value (default 4, i.e.
  at least 6 paired observations).

- shrinkage:

  logical; apply `ashr` shrinkage (default `TRUE`).

- alpha:

  numeric; `ashr` alpha parameter (default 0).

- MHC_direction:

  character; `"x"` or `"y"` — direction of multiple-hypothesis
  correction (Benjamini-Hochberg). With `"y"`, p-values are adjusted
  across the features of `X` separately for each column of `Y`; with
  `"x"`, across the columns of `Y` separately for each feature of `X`.
  Defaults to `"x"` when `length(Y) >= length(X)`, otherwise `"y"` (same
  default as
  [`cdsrmodels::lin_associations`](https://rdrr.io/pkg/cdsrmodels/man/lin_associations.html)).

## Value

Named list with elements `N`, `rho`, `beta`, `beta.se`, `p.val`,
`q.val`, and `res.table` (a `data.frame` from `ashr`).
