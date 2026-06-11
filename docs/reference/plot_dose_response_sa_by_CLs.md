# Plot drug response curves for single-agent data for selected call lines and drugs

Plot drug response curves for single-agent data for selected call lines
and drugs

## Usage

``` r
plot_dose_response_sa_by_CLs(
  dt_metrics,
  dt_average,
  cellline_name_vec = NULL,
  drug_name_vec = NULL,
  normalization_type = "GR",
  colors_vec = NULL,
  plot_averaged_flag = TRUE,
  plot_fit_flag = TRUE,
  fit_source = "gDR"
)
```

## Arguments

- dt_metrics:

  data.table representing data from the `Metrics` assay, outputted by
  `gDRutils::convert_se_assay_to_dt(se, "Metrics")` and single-agent
  `SummarizedExperiment`

- dt_average:

  data.table representing data from the `Averaged` assay, outputted by
  `gDRutils::convert_se_assay_to_dt(se, "Averaged")` and single-agent
  `SummarizedExperiment`

- cellline_name_vec:

  character vector with cell line names to be plotted (Cell Line Name)

- drug_name_vec:

  character vector with drug names to be plotted (Drug Name)

- normalization_type:

  string with normalization_types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- colors_vec:

  character vector with colors for `group_names` - name or hex value

- plot_averaged_flag:

  logical flag whether plot points with average values

- plot_fit_flag:

  logical flag whether plot points with fitted values

- fit_source:

  string source name for metrics

## Value

list of `ggplot objects` containing plots of dose-response curves

## Examples

``` r
mae <- gDRutils::get_synthetic_data("small")
se <- mae[[gDRutils::get_supported_experiments("sa")]]
dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
cellline_name_vec <- unique(dt_metrics[["CellLineName"]])[2:5]
drug_name_vec <- unique(dt_metrics[["DrugName"]])[5:7]

plot_dose_response_sa_by_CLs(dt_metrics = dt_metrics,
                             dt_average = dt_average,
                             cellline_name_vec = cellline_name_vec,
                             drug_name_vec = drug_name_vec,
                             normalization_type = "RV",
                             colors_vec = c("#00008B", "#FF6347", "#4CBB17"))
#> $drug_006

#> 
#> $drug_007

#> 
#> $drug_008

#> 
```
