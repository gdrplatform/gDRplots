# Load DepMap merged data for one selected feature

Load DepMap merged data for one selected feature

## Usage

``` r
prep_dt_depmap_feat(
  feat_data_path,
  meta_data_path,
  feature_set = "CRISPRGeneEffect",
  with_decoding = FALSE
)
```

## Arguments

- feat_data_path:

  string with path to the directory containing the molecular feature set
  file to load from DepMap.

- meta_data_path:

  string with path to metadata file describing all cancer models/cell
  lines which are referenced by a dataset contained within the DepMap
  portal. It is usually a file named `Model.csv` or `Model.csv.gz`.

- feature_set:

  string containing the name of the molecular feature set to load from
  DepMap. This name should also correspond to the file containing the
  feature data (without the extension, which is assumed to be `csv` or
  `csv.gz`)

- with_decoding:

  logical whether the feature OmicsArmLevelCNA,
  OmicsSomaticMutationsMatrixHotspot and
  OmicsSomaticMutationsMatrixDamaging should be encoded into a 0-1
  scheme

## Value

A named list with elements, that may be input to
[`prep_dt_assoc`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_assoc.md)

- `dt_depmap` `data.table` with feature data from DepMap (wide format),.

- `selected_feat_meta_col` string name of feature.

## Examples

``` r
if (FALSE) { # \dontrun{
feat_data_path <- file.path(".", "depmapdata")
meta_data_path <- file.path(".", "Model.csv")
dt_depmap_feat <- prep_dt_depmap_feat(feat_data_path = feat_data_path,
                                      meta_data_path = meta_data_path)
} # }
```
