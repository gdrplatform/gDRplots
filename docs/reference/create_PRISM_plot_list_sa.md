# Create a nested list of plots for PRISM data with single-agent metrics

Create a nested list of plots for PRISM data with single-agent metrics

## Usage

``` r
create_PRISM_plot_list_sa(
  drug_name_vec,
  dt_metrics,
  dt_average = NULL,
  normalization_type_vec = "RV",
  metric = c("xc50", "x_mean", "x_max"),
  fit_source = "gDR",
  meta_data_path,
  feat_data_path,
  feature_sets,
  metadata_columns = NULL,
  clear_taxonomy_info = TRUE,
  with_decoding = FALSE
)
```

## Arguments

- drug_name_vec:

  character vector with drug names to be plotted (identifiers
  `DrugName`)

- dt_metrics:

  `data.table` representing data from the `Metrics` assay, outputted by
  `gDRutils::convert_se_assay_to_dt(se, "Metrics")` and single-agent
  `SummarizedExperiment`

- dt_average:

  `data.table` representing data from the `Averaged` assay, outputted by
  `gDRutils::convert_se_assay_to_dt(se, "Averaged")` and single-agent
  `SummarizedExperiment`

- normalization_type_vec:

  character vector with normalization types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability") or both

- metric:

  character vector with names of metric; chosen from: "xc50" ("GR50" or
  "IC50" - respectively depending on `normalization_type`), "x_max" ("GR
  Max" or "E Max") or "x_mean" ("GR Mean" or "RV Mean")

- fit_source:

  string source name for metrics

- meta_data_path:

  string path to metadata file describing all cancer models/cell lines
  which are referenced by a dataset contained within the DepMap portal.
  It is usually a file named `Model.csv` or `Model.csv.gz`.

- feat_data_path:

  string path to the directory containing the molecular feature set file
  to load from DepMap.

- feature_sets:

  character vector containing the names of the molecular feature sets to
  load from DepMap. These names should also correspond to the file names
  containing the feature data (without the extension, which is assumed
  to be `csv` or `csv.gz`)

- metadata_columns:

  character vector with the metadata columns to load for DepMap cell
  lines

- clear_taxonomy_info:

  logical flag whether to remove taxonomy information for gene names in
  table with the molecular feature sets from DepMap.

- with_decoding:

  logical whether the feature OmicsArmLevelCNA,
  OmicsSomaticMutationsMatrixHotspot and
  OmicsSomaticMutationsMatrixDamaging should be encoded into a 0-1
  scheme

## Value

A named list with elements:

- `ls_plot` nested list of plots for selected type of experiment

- `ls_assoc_data` nested list of table with association data

## Author

Janina Smoła <janina.smola@contractors.roche.com>
