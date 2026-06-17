# Load DepMap merged data for one selected metadata

Load DepMap merged data for one selected metadata

## Usage

``` r
prep_dt_depmap_meta(meta_data_path, metadata_col = "PatientRace")
```

## Arguments

- meta_data_path:

  string with path to metadata file describing all cancer models/cell
  lines

## Value

A named list with elements, that may be input to
[`prep_dt_assoc`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_assoc.md)

- `dt_depmap` `data.table` with feature data from DepMap (wide format),

- `selected_feat_meta_col` string name of metadata column.

## Examples

``` r
if (FALSE) { # \dontrun{
meta_data_path <- file.path(".", "Model.csv")
dt_depmap_meta <- prep_dt_depmap_meta(meta_data_path)
} # }
```
