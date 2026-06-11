# Compute value of param cluster_rows or cluster_cols in pheatmap::pheatmap

The `cluster_rows` and `cluster_cols` parameters pheatmap::pheatmap may
take values: - boolean determining if rows/columns should be clustered -
`hclust` object, this function allows to calculate proper value taking
into account matrix, additional condition and selected function used to
compute the distance in `hclust` object

## Usage

``` r
.get_pheatmap_cluster_param(
  mat_to_cluster,
  distfun = stats::dist,
  additional_condition = TRUE
)
```

## Arguments

- mat_to_cluster:

  numeric matrix to be clustered; cluster dimension must be named

- distfun:

  function used to compute the distance (dissimilarity) between rows;
  defaults to [`stats::dist`](https://rdrr.io/r/stats/dist.html) using
  euclidean euclidean.

- additional_condition:

  additional logical flag whether rows/columns should be clustered

## Value

logical flag determining if rows should be clustered or `hclust` object.

## Details

To compute the correct value when clustering columns - use the
transposed source matrix as `mat_to_cluster`

## See also

[`compute_distances`](compute_distances.md)

## Author

Janina Smoła <janina.smola@contractors.roche.com>

## Examples

``` r
if (FALSE) { # \dontrun{
mat <- matrix(1:24, nrow = 4)
rownames(mat) <- sprintf("row_%s", 1:4)
colnames(mat) <- sprintf("col_%s", 1:6)
.get_pheatmap_cluster_param(mat)
.get_pheatmap_cluster_param(t(mat))
.get_pheatmap_cluster_param(t(mat), distfun = compute_distances)

mat[2,2] <- NA
mat[2,1] <- Inf
.get_pheatmap_cluster_param(mat)
.get_pheatmap_cluster_param(mat, distfun = compute_distances)
.get_pheatmap_cluster_param(t(mat), distfun = compute_distances)
add_cond <- NCOL(mat) > 10
.get_pheatmap_cluster_param(mat, distfun = compute_distances, additional_condition = add_cond)
} # }
```
