# Create a nested list of plots for PRISM data with combo metrics

Create a nested list of plots for PRISM data with combo metrics

## Usage

``` r
create_PRISM_plot_list_combo(
  drug1_name_vec,
  drug2_name_vec,
  dt_metrics,
  dt_scores = NULL,
  normalization_type_vec = "RV",
  metric = c("xc50", "x_mean", "x_max"),
  metric_scores = c("hsa_score", "bliss_score"),
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

- drug1_name_vec:

  character vector with drug names to be plotted (identifiers
  `DrugName`)

- drug2_name_vec:

  character vector with co-drug names to be plotted (identifiers
  `DrugName_2`)

- dt_metrics:

  `data.table` representing data from the `Metrics` assay, outputted by
  `gDRutils::convert_se_assay_to_dt(se, "Metrics")` and combo
  `SummarizedExperiment`

- dt_scores:

  `data.table` representing data from the `scores` assay, outputted by
  `gDRutils::convert_se_assay_to_dt(se, "scores")` and combo
  `SummarizedExperiment`

- normalization_type_vec:

  character vector with normalization types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability") or both

- metric:

  character vector with names of metric; chosen from: "xc50" ("GR50" or
  "IC50" - respectively depending on `normalization_type`), "x_max" ("GR
  Max" or "E Max") or "x_mean" ("GR Mean" or "RV Mean")

- metric_scores:

  character vector with names of combo metric; chosen from:
  "hsa_score"("Bliss Excess GR" or "Bliss Excess RV" - respectively
  depending on `normalization_type`), "bliss_score" ("Bliss Score GR" or
  "Bliss Score RV")

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

Janina Smoła <janina.smola@external.roche.com>

## Examples

``` r
mae <- qs2::qs_read(system.file("testdata/finalMAE_combo_matrix_small.qs2",
                                 package = "gDRtestData"))
#> Loading required namespace: MultiAssayExperiment
se_combo <- mae[[gDRutils::get_supported_experiments("combo")]]
dt_metrics <- gDRutils::convert_se_assay_to_dt(se_combo, "Metrics")
#> Loading required namespace: BumpyMatrix
drug1 <- unique(dt_metrics[[gDRutils::get_env_identifiers("drug_name")]])[1]
drug2 <- unique(dt_metrics[[gDRutils::get_env_identifiers("drug_name2")]])[1]
create_PRISM_plot_list_combo(
  drug1_name_vec = drug1,
  drug2_name_vec = drug2,
  dt_metrics = dt_metrics,
  meta_data_path = system.file("depmap_data/Model.csv.gz", package = "gDRtestData"),
  feat_data_path = system.file("depmap_data", package = "gDRtestData"),
  feature_sets = c("CRISPRGeneEffect")
)
#> $ls_plot
#> $ls_plot$CRISPRGeneEffect
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.001_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.001_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.00316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.00316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.01_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.01_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.0316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.0316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.1_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.1_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_1_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_1_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_3.16_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_3.16_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.001_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.001_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.00316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.00316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.01_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.01_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.0316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.0316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.1_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.1_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_1_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_1_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_3.16_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_3.16_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.001_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.001_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.00316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.00316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.01_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.01_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.0316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.0316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.1_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.1_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_1_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_1_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_3.16_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_3.16_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.001_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.001_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.00316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.00316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.01_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.01_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.0316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.0316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.1_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.1_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_1_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_1_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_3.16_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_3.16_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.001_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.001_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.00316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.00316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.01_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.01_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.0316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.0316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.1_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.1_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_1_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_1_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_3.16_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_3.16_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.001_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.001_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.00316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.00316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.01_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.01_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.0316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.0316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.1_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.1_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_1_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_1_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_3.16_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_3.16_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.001_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.001_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.00316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.00316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.01_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.01_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.0316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.0316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.1_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.1_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_1_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_1_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_3.16_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_3.16_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.001_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.001_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.00316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.00316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.01_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.01_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.0316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.0316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.1_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.1_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_1_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_1_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_3.16_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_3.16_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.001_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.001_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.00316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.00316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.01_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.01_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.0316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.0316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.1_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.1_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.316_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.316_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_1_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_1_drug_2

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_3.16_drug_1

#> 
#> $ls_plot$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_3.16_drug_2

#> 
#> 
#> 
#> 
#> 
#> $ls_assoc_data
#> $ls_assoc_data$CRISPRGeneEffect
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.001_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.001_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.00316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.00316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.01_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.01_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.0316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.0316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.1_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.1_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_0.316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_1_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_1_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_3.16_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_zero_3.16_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.001_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.001_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.00316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.00316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.01_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.01_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.0316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.0316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.1_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.1_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_0.316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_1_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_1_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_3.16_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_zero_3.16_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.001_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.001_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.00316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.00316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.01_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.01_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.0316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.0316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.1_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.1_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_0.316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_1_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_1_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_3.16_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_zero_3.16_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.001_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.001_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.00316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.00316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.01_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.01_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.0316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.0316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.1_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.1_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_0.316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_1_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_1_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_3.16_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_3.16_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.001_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.001_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.00316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.00316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.01_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.01_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.0316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.0316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.1_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.1_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_0.316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_1_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_1_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_3.16_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_3.16_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.001_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.001_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.00316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.00316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.01_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.01_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.0316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.0316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.1_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.1_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_0.316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_1_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_1_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_3.16_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_3.16_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.001_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.001_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.00316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.00316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.01_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.01_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.0316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.0316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.1_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.1_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_0.316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_1_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_1_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_3.16_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_log10_xc50_cotrt_diff_3.16_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.001_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.001_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.00316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.00316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.01_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.01_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.0316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.0316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.1_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.1_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_0.316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_1_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_1_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_3.16_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_mean_cotrt_diff_3.16_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.001_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.001_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.00316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.00316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.01_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.01_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.0316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.0316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.1_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.1_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.316_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_0.316_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_1_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_1_drug_2
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_3.16_drug_1
#> NULL
#> 
#> $ls_assoc_data$CRISPRGeneEffect$`drug_004 x drug_021`$RV$RV_gDR_x_max_cotrt_diff_3.16_drug_2
#> NULL
#> 
#> 
#> 
#> 
#> 
```
