# Compute pairwise linear associations between columns of X and Y

Inline reimplementation of `lin_associations` from the
[cdsrmodels](https://github.com/cancerdatasci/cdsrmodels) package (MIT
License, Broad Institute / DepMap Portal) to remove the GitHub-only
package dependency is replaced. Uses
[`WGCNA::cor`](https://rdrr.io/pkg/WGCNA/man/cor.html) for fast pairwise
correlation on large matrices and
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

  integer; minimum number of finite paired observations required to
  compute a p-value (default 4).

- shrinkage:

  logical; apply `ashr` shrinkage (default `TRUE`).

- alpha:

  numeric; `ashr` alpha parameter (default 0).

- MHC_direction:

  character; `"x"` or `"y"` — direction of multiple-hypothesis
  correction. Defaults to `"y"` when `ncol(Y) >= ncol(X)`, otherwise
  `"x"`.

## Value

Named list with elements `N`, `rho`, `beta`, `beta.se`, `p.val`,
`q.val`, and `res.table` (a `data.frame` from `ashr`).
