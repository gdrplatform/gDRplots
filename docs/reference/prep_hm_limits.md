# Calculate limit for combo heatmap with gDR assumptions

Calculate limit for combo heatmap with gDR assumptions

## Usage

``` r
prep_hm_limits(
  num_vec,
  scale_type = "smooth",
  normalization_type = "GR",
  symmetric = FALSE
)
```

## Arguments

- num_vec:

  numeric vector

- scale_type:

  string, `"smooth"` for response values or `"excess"` for differences

- normalization_type:

  string with normalization_types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- symmetric:

  logical indicating if limits should be symmetric around 0

## Value

capped limits (min and max) for given numeric vector

## Examples

``` r
# \donttest{
vec <- c(-0.1, -0.3, 0, 0.5, Inf, NA)
gDRplots:::prep_hm_limits(vec)
#> [1] -0.3  1.0
gDRplots:::prep_hm_limits(vec, scale_type = "excess", symmetric = TRUE)
#> [1] -0.5  0.5
# }
```
