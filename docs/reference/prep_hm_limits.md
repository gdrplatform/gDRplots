# Calculate limit for combo heatmap with gDR assumptions

Calculate limit for combo heatmap with gDR assumptions

## Usage

``` r
prep_hm_limits(
  num_vec,
  metric = "smooth",
  normalization_type = "GR",
  symmetric = FALSE
)
```

## Arguments

- num_vec:

  numeric vector

- metric:

  string name of the combo metric; one of: "smooth" ("Smooth GR" or
  "Smooth RV" - respectively depending on `normalization_type`)
  "hsa_excess" ("Bliss Excess GR" or "Bliss Excess RV") or
  "bliss_excess" ("Bliss Excess GR" or "Bliss Excess RV")

- normalization_type:

  string with normalization_types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- symmetric:

  logical indicating if limits should be symmetric around 0

## Value

capped limits (min and max) for given numeric vector

## Examples

``` r
if (FALSE) { # \dontrun{
vec <- c(-0.1, -0.3, 0, 0.5, Inf, NA)
prep_hm_limits(vec)
prep_hm_limits(vec, metric = "hsa_excess", symmetric = TRUE)
} # }
```
