# Help function

This function retrieves information about the feature name and drug name
from a file name, which has the schema
\<chunk_name\>\_\_\<feature_name\>\_\<drug_name\>\_\<gDR_metric_name\>
Assumption: - \<feature_name\> is assumed not to contain any underscores
(usually written in camel case) - \<gDR_metric_name\> starts with
normalization type that is `RV` or `GR`

## Usage

``` r
.get_info_from_name(file_name, normalization_type = "RV")
```

## Arguments

- file_name:

  A string with file name

- normalization_type:

  string with normalization types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

## Value

A named list with elements:

- `drug_grid` a string with drug name

- `feat_meta` a string with omic name (feature or meta)

## Author

Janina Smoła <janina.smola@contractors.roche.com>

## Examples

``` r
if (FALSE) { # \dontrun{
f_n <- "name_chunk__FEAT_DRUG_ABC_RV_gDR_log10_xc50.xlsx"
.get_info_from_name(f_n)
} # }
```
