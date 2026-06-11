# Plot panel of heatmaps with fitted and reference data for isobolograms

This function is dedicated to cases in which given cell lines are
exposed to drugs of different concentrations and have almost no or no
common values.

## Usage

``` r
heatmap_combo_with_isoref_panel_independent(
  dt_excess,
  dt_isobolograms,
  drug1_name,
  drug2_name,
  cl_names,
  normalization_type = "GR",
  metric = "smooth",
  iso_levels = "0.5",
  colors_vec = NULL,
  no_breaks = 50,
  swap_axes = FALSE
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

- cl_names:

  character vector with cell line names to be plotted (Cell Line Name);
  if `NULL` - all available cell lines will be plotted

- normalization_type:

  string with normalization_types to be selected one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability")

- metric:

  string name of the combo metric; one of: "smooth" ("Smooth GR" or
  "Smooth RV" - respectively depending on `normalization_type`)
  "hsa_excess" ("Bliss Excess GR" or "Bliss Excess RV") or
  "bliss_excess" ("Bliss Excess GR" or "Bliss Excess RV")

- iso_levels:

  character vector with isobologram levels to be selected; when `NULL` -
  no isolines will be displayed

- colors_vec:

  character vector of colors (valid name or hex) used in heatmap; the
  default depends on `metric`: for "smooth" - the dark purple-light grey
  palette and for "hsa_excess" and "bliss_excess" - the blue-light
  grey-red color scale

- no_breaks:

  numeric number of breaks on scale

- swap_axes:

  logical flag whether to swap the axes with drugs of the heatmap

## Value

`ggplot` object containing panel with heatmaps for fitted values and
reference data for isobolograms for selected drug and co-drug by list of
cell lines

## Examples

``` r
cl_names <-
  c("cellline_AA", "cellline_EA", "cellline_IB",
  "cellline_MC", "cellline_BC", "cellline_FD")

drug1_name <- "drug_001"
drug2_name <- "drug_026"

mae <- gDRutils::get_synthetic_data("combo_matrix")
se <- mae[[gDRutils::get_supported_experiments("combo")]]
dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")

heatmap_combo_with_isoref_panel_independent(dt_excess,
                                            dt_isobolograms,
                                            drug1_name, drug2_name,
                                            cl_names)


heatmap_combo_with_isoref_panel_independent(dt_excess,
                                            dt_isobolograms,
                                            drug1_name, drug2_name,
                                            cl_names,
                                            iso_levels = c("-0.25", "0.25"))


heatmap_combo_with_isoref_panel_independent(
  dt_excess,
  dt_isobolograms,
  drug1_name, drug2_name,
  cl_names = c("cellline_FD", "cellline_MC", "cellline_AA", "cellline_EA"),
  iso_levels =  c("-0.25", "-0.05", "0.2", "0.65"))


heatmap_combo_with_isoref_panel_independent(dt_excess,
                                            dt_isobolograms,
                                            drug1_name, drug2_name,
                                            cl_names,
                                            normalization_type = "RV",
                                            iso_levels = NULL,
                                            colors_vec = c("darkcyan", "snow", "darkorange"),
                                            swap_axes = TRUE)

```
