# Get Legend Title

Get Legend Title

## Usage

``` r
get_hm_title(metric = "xc50", normalization_type = "GR", dataset_name = NULL)
```

## Arguments

- metric:

  string name of the metric one of: "xc50"("GR50" or "IC50" -
  respectively depending on `normalization_type`), "x_max" ("GR Max" or
  "E Max") or x_mean" ("GR Mean" or "RV Mean")

- normalization_type:

  string with normalization types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- dataset_name:

  string name of dataset

## Value

character title for heatmap

## Examples

``` r
get_hm_title(dataset_name = "Dateset DX123",
             metric = "x_mean",
             normalization_type = "GR")
#> [1] "Dateset DX123 (GR Mean)"

get_hm_title(metric = "xc50",
             normalization_type = "GR")
#> [1] "log10(GR50)"
```
