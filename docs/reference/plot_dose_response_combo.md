# Plot drug response curves for combo data

Plot drug response curves for combo data

## Usage

``` r
plot_dose_response_combo(
  dt_average,
  drug1_name,
  drug2_name,
  cl_name,
  normalization_type = "GR",
  colors_vec = NULL,
  split_by_conc = FALSE
)
```

## Arguments

- dt_average:

  data.table representing data from the `Averaged` assay, outputted by
  `gDRutils::convert_se_assay_to_dt(se, "Averaged")` and combo
  `SummarizedExperiment`

- drug1_name:

  string with drug name to be plotted (identifiers `DrugName`)

- drug2_name:

  string with co-drug name to be plotted (identifiers `DrugName_2`)

- cl_name:

  string with cell line name to be plotted (identifiers `CellLineName`)

- normalization_type:

  string with normalization_types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- colors_vec:

  character vector with colors for `group_names` - name or hex value
  note that the first color will be assigned to the min value of
  `Concentration_2`, and the last one - to the max of `Concentration_2`;
  the default is the orange-black palette

- split_by_conc:

  split_by_conc logical flag indicating whether curves for
  `Concentration_2` should be plotted on a single plot or separately

## Value

`ggplot` object containing plot with dose-response curves for combo data

## Examples

``` r
mae <- gDRutils::get_synthetic_data("combo_matrix")
se <- mae[[gDRutils::get_supported_experiments("combo")]]
dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")

cl_name <- "cellline_BC"
drug1_name <- "drug_011"
drug2_name <- "drug_021"

plot_dose_response_combo(dt_average = dt_average,
                         drug1_name = drug1_name,
                         drug2_name = drug2_name,
                         cl_name = cl_name)


plot_dose_response_combo(dt_average = dt_average,
                         drug1_name = drug1_name,
                         drug2_name = drug2_name,
                         cl_name = cl_name,
                         split_by_conc = TRUE)


plot_dose_response_combo(dt_average = dt_average,
                         drug1_name = drug1_name,
                         drug2_name = drug2_name,
                         cl_name = cl_name,
                         normalization_type = "RV",
                         colors_vec = c("lightblue", "darkblue"))

```
