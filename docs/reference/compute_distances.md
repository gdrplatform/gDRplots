# Compute distance between rows of a matrix.

This function was lifted from the package `factoextra` and slightly
adjusted. The default method was changed to "spearman" and an option was
added to replace missing values in the resulting distance matrix with an
arbitrary value. This option defaults to 0 so that the function can be
called by `iheatmapr::add_row_clustering` and
`iheatmapr::add_col_clustering`.

## Usage

``` r
compute_distances(
  x,
  method = "spearman",
  use = "pairwise.complete.obs",
  stand = FALSE,
  dummy = 0,
  ...
)
```

## Arguments

- x:

  numeric matrix

- method:

  character string specifying the distance measure to be used; must be
  on of: "euclidean", "maximum", "manhattan", "canberra", "binary",
  "minkowski", "pearson", "spearman" or "kendall"

- use:

  character string specifying the method for computing covariances in
  the presence of missing values. Only used when "pearson", "spearman"
  or "kendall" chosen as a distance measure

- stand:

  logical flag specifying whether the data should be standardized; if
  TRUE, columns are converted to z scores using `scale`

- dummy:

  value to substitute for missing values; defaults to 0

- ...:

  arguments passed to internally called functions

## Value

Object of class `dist`. NA and NaN are substituted by the value of
`dummy`.

## Control

The function offers different distance measures but since it is called
internally, it will always use default arguments. As a result, should
you want to do anything other than compute Spearman correlation distance
and replace NAs with 0, you must edit this function definition
accordingly. At the moment (December 2020) this is sufficient. If the
control were to be given to the user, remove defaults for this
definition and create a wrapper within `plotly_metric_clustering`.

The function was originally named `computeDistances`

## See also

\[factoextra::get_dist()\]

## Author

Alboukadel Kassambara <alboukadel.kassambara@gmail.com>

## Examples

``` r
x <- matrix(1:9, nrow = 3, ncol = 3)
rownames(x) <- letters[seq(NROW(x))]
compute_distances(x)
#>   a b
#> b 0  
#> c 0 0
```
