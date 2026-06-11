# Plot box plots for metric for single-agent data grouped by selected variable

Plot box plots for metric for single-agent data grouped by selected
variable

## Usage

``` r
plot_boxplot_metric_sa(
  dt_metrics,
  group_var,
  normalization_type = "GR",
  metric = "xc50",
  fit_source = "gDR",
  grouped_flag = FALSE,
  colored_pts_flag = FALSE,
  colors_vec = NULL,
  with_inf = FALSE
)
```

## Arguments

- dt_metrics:

  data.table representing data from the `Metrics` assay, outputted by
  `gDRutils::convert_se_assay_to_dt(se, "Metrics")` and single-agent
  `SummarizedExperiment`

- group_var:

  string name of group variable; one of: `"CellLineName"` or
  `"DrugName"`

- normalization_type:

  string with normalization types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- metric:

  string name of the metric; one of: "xc50" ("GR50" or "IC50" -
  respectively depending on `normalization_type`), "x_max" ("GR Max" or
  "E Max") or "x_mean" ("GR Mean" or "RV Mean"), but the values from any
  numeric column can be displayed

- fit_source:

  string source name for metrics

- grouped_flag:

  a logical flag whether the boxplots should be grouped and colored by
  `Tissue` for `group_var` set as `"CellLineName"` and `drug_moa` - for
  `"DrugName"`

- colored_pts_flag:

  a logical flag whether the points should be colored by grouped
  variable - for `group_var` equal `"CellLineName"` points will be
  colored by `"DrugName"` and similarly vice versa

- colors_vec:

  character vector with colors (name or hex value) to color boxplots;
  for `grouped_flag` set as `FALSE` only first from vector will be used

- with_inf:

  a logical flag indicating whether infinite values should be shown on
  boxplots

## Value

`ggplot` object containing boxplots for selected single-agent metric
grouped by `group_var`

## Author

Janina Smoła <janina.smola@contractors.roche.com>

## Examples

``` r
mae <- gDRutils::get_synthetic_data("combo_matrix")
se <- mae[[gDRutils::get_supported_experiments("sa")]]

dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")

plot_boxplot_metric_sa(dt_metrics,
                       group_var = "CellLineName")


plot_boxplot_metric_sa(dt_metrics,
                       group_var = "DrugName")


plot_boxplot_metric_sa(dt_metrics,
                       group_var = "DrugName",
                       metric = "x_AOC_range",
                       colors_vec = "grey",
                       colored_pts_flag = TRUE)


plot_boxplot_metric_sa(dt_metrics,
                       group_var = "CellLineName",
                       normalization_type = "RV",
                       metric = "x_max",
                       colors_vec = c("gold", "darkorange", "darkcyan", "darkblue"),
                       grouped_flag = TRUE)

```
