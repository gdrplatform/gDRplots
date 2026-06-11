# Plot line plot of combination index

Plot line plot of combination index

## Usage

``` r
plot_combination_index(
  dt_excess,
  dt_isobolograms,
  drug1_name,
  drug2_name,
  cl_name,
  normalization_type = "GR",
  iso_levels = c("0.25", "0.5", "0.75"),
  colors_vec_iso = NULL
)
```

## Arguments

- dt_excess:

  data.table representing data from the `excess` assay, outputted by
  `gDRutils::convert_se_assay_to_dt(se, "excess")` and combo
  `SummarizedExperiment`

- dt_isobolograms:

  data.table representing data from the `isobolograms` assay, outputted
  by `gDRutils::convert_se_assay_to_dt(se, "isobolograms")` and combo
  `SummarizedExperiment`

- drug1_name:

  string with drug name to be plotted (identifiers `DrugName`)

- drug2_name:

  string with co-drug name to be plotted (identifiers `DrugName_2`)

- cl_name:

  string with cell line to be plotted (identifiers `CellLineName`)

- normalization_type:

  string with normalization_types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- iso_levels:

  character vector with isobologram levels to be selected

- colors_vec_iso:

  character vector of colors (valid name or hex) used for the isolines;
  the default is the dark red-orange palette

## Value

`ggplot` object containing combination index plot at different ratios of
the two drugs

## Examples

``` r
cl_name <- "cellline_BC"
drug1_name <- "drug_001"
drug2_name <- "drug_026"

mae <- gDRutils::get_synthetic_data("combo_matrix")
se <- mae[[gDRutils::get_supported_experiments("combo")]]
dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")

plot_combination_index(dt_excess,
                       dt_isobolograms,
                       drug1_name, drug2_name,
                       cl_name,
                       normalization_type = "GR")


plot_combination_index(dt_excess,
                       dt_isobolograms,
                       drug1_name, drug2_name,
                       cl_name,
                       normalization_type = "RV",
                       colors_vec_iso = c("darkblue", "darkcyan"))


cl_name <- "cellline_JE"
drug1_name <- "drug_011"
drug2_name <- "drug_026"

plot_combination_index(dt_excess,
                       dt_isobolograms,
                       drug1_name, drug2_name,
                       cl_name,
                       normalization_type = "RV",
                       iso_levels = "0.5")

```
