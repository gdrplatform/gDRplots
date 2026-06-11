# Plot pretty heatmap with annotations for single-agent data

Plot pretty heatmap with annotations for single-agent data

## Usage

``` r
pheatmap_with_anno_sa(
  dt_metrics,
  dt_metrics_capped = NULL,
  normalization_type = "GR",
  metric = "xc50",
  fit_source = "gDR",
  hm_title = NA,
  colors_vec = NULL,
  no_breaks = 50,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  distfun = compute_distances,
  annotation_row = NULL,
  annotation_col = NULL,
  annotation_colors = NULL,
  max_hm_lbl_length = gDRutils::get_settings_from_json("MAX_HM_LBL_LENGTH",
    system.file(package = "gDRplots", "settings.json"))
)
```

## Arguments

- dt_metrics:

  `data.table` representing data from the `"Metrics"` assay, outputted
  by
  [`gDRutils::convert_se_assay_to_dt`](https://gdrplatform.github.io/gDRstyle/reference/convert_se_assay_to_dt.html)
  and single-agent `SummarizedExperiment`

- dt_metrics_capped:

  `data.table` representing data from the `"Metrics"` assay, the same as
  `dt_metrics` but with capped values with
  [`gDRutils::cap_assay_infinities`](https://gdrplatform.github.io/gDRstyle/reference/cap_assay_infinities.html)

- normalization_type:

  string with normalization types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- metric:

  string name of the metric; one of: "xc50" ("GR50" or "IC50" -
  respectively depending on `normalization_type`), "x_max" ("GR Max" or
  "E Max") or "x_mean" ("GR Mean" or "RV Mean"); but the values from any
  numeric column can be displayed.

- fit_source:

  string source name for metrics

- hm_title:

  string plot title

- colors_vec:

  character vector of colors (valid name or hex) used in heatmap; note
  that the first color will be assigned to the min value, and the last
  one - to the max

- no_breaks:

  numeric number of breaks on scale used for mapping values to colors

- cluster_rows:

  logical flag whether rows should be clustered; the dendrogram will not
  be shown for the matrix with any dimension greater than 200.

- cluster_cols:

  logical flag whether columns should be clustered; the dendrogram will
  not be shown for the matrix with any dimension greater than 200.

- distfun:

  function used to compute the distance (dissimilarity) between both
  rows and columns; used for the dendrogram when `cluster_rows` or
  `cluster_cols` is set to TRUE the default is
  [`compute_distances`](compute_distances.md) using Spearman method.

- annotation_row:

  `data.table` that specifies the annotations shown on left side of the
  heatmap. Each row defines the features for a specific row. The rows in
  the data and in the annotation are matched using corresponding names
  from the required `DrugName` column. Note that color schemes takes
  into account if variable is continuous or discrete.

- annotation_col:

  `data.table` that specifies the annotations shown above the heatmap.
  Each row defines the features for a specific column. The columns in
  the data and in the annotation are matched using corresponding names
  from the required `CellLineName` column. Note that color schemes takes
  into account if variable is continuous or discrete.

- annotation_colors:

  named list for specifying `annotation_col` and `annotation_row` track
  colors manually; note list is named with annotation name (column names
  of `annotation_row` - without `DrugName` and column names of
  `annotation_col` - without `CellLineName`), each list item is named
  vector with valid color name for each value described in
  `annotation_row` and in `annotation_col` - respectively. Not described
  elements will be colored in default.

- max_hm_lbl_length:

  numeric value for the maximum number of characters in the label; if
  set to Inf, no trimming will be performed; for better readability, it
  is recommended to use the default number.

## Value

A named list with elements:

- `data` a list containing the information visualized in the heatmap:

  - `matrix` data shown in the heatmap for the selected metric.

  - `annotation_row` a table with row annotations (for `DrugName`), if
    provided.

  - `annotation_col` a table with column annotations (for
    `CellLineName`), if provided.

- `pheatmap` object containing the heatmap itself.

## See also

[`pheatmap::pheatmap`](https://rdrr.io/pkg/pheatmap/man/pheatmap.html)

## Author

Janina Smoła <janina.smola@contractors.roche.com>

## Examples

``` r
mae <- gDRutils::get_synthetic_data("combo_matrix")
se <- mae[[gDRutils::get_supported_experiments("sa")]]
dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
                                               assay_name = "Metrics")
dt_averaged <- gDRutils::convert_se_assay_to_dt(se = se,
                                                assay_name = "Averaged")
dt_metrics_capped <-
  gDRutils::cap_assay_infinities(
    conc_assay_dt = dt_averaged,
    assay_dt = dt_metrics,
    experiment_name = gDRutils::get_supported_experiments("sa"),
    col = "xc50",
    capping_fold = 5)

output <- pheatmap_with_anno_sa(dt_metrics = dt_metrics)
hm_0 <- output[["heatmap"]]
ggpubr::as_ggplot(hm_0[["gtable"]])


output <- pheatmap_with_anno_sa(dt_metrics = dt_metrics,
                                dt_metrics_capped = dt_metrics_capped)
hm_1 <- output[["heatmap"]]
ggpubr::as_ggplot(hm_1[["gtable"]])


annotation_manual_col <-
  unique(dt_metrics[, c("CellLineName", "Tissue"), with = FALSE])
annotation_manual_row <-
  unique(dt_metrics[, c("DrugName", "drug_moa"), with = FALSE])
annotation_map <-
  get_ann_color_map(unique(dt_metrics[, c("Tissue", "drug_moa"), with = FALSE]))

output <- pheatmap_with_anno_sa(dt_metrics = dt_metrics,
                                normalization_type = "RV",
                                metric = "x_mean",
                                colors_vec = c("darkblue", "grey90"),
                                annotation_row = annotation_manual_row,
                                annotation_col = annotation_manual_col,
                                annotation_colors = annotation_map)
hm_2 <- output[["heatmap"]]
ggpubr::as_ggplot(hm_2[["gtable"]])


annotation_manual <- data.table::data.table(
  CellLineName =
    c("cellline_AA", "cellline_EA", "cellline_IB", "cellline_MC", "cellline_BC"),
  mut_A = c(1, 1, 1, 0, 0),
  mut_B = c("yes", "yes", "no", "no", "no")
)
annotation_map <- list(
  mut_A = c("1" = "coral", "0" = "cadetblue"),
  mut_B = c("yes" = "black", "no" = "grey90")
)

output <- pheatmap_with_anno_sa(dt_metrics = dt_metrics,
                                annotation_col = annotation_manual,
                                annotation_colors = annotation_map,
                                hm_title = get_hm_title(
                                  normalization_type = "GR",
                                  metric = "hsa_score",
                                  dataset_name = "Combo Matrix - combo data"))
hm_3 <- output[["heatmap"]]
ggpubr::as_ggplot(hm_3[["gtable"]])

```
