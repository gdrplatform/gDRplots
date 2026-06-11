# Plot heatmap of mapping controls to treated for single-agent and combo data to control quality of the data

Plot heatmap of mapping controls to treated for single-agent and combo
data to control quality of the data

## Usage

``` r
heatmap_control_mapping_qc(dt_treat, dt_controls)
```

## Arguments

- dt_treat:

  data.table representation of the data in `RawTreated` assay, outputted
  by `gDRutils::convert_se_assay_to_dt(se, "RawTreated")` and
  `SummarizedExperiment` with chosen data type: single-agent or combo

- dt_controls:

  data.table representation of the data in `Controls` assay, outputted
  by `gDRutils::convert_se_assay_to_dt(se, "Controls")` and
  `SummarizedExperiment` with chosen data type: single-agent or combo

## Value

`pheatmap` object containing hetamap of mapping controls to treated

## Author

Bartosz Czech <czech.bartosz@external.gene.com>

## Examples

``` r
mae <- gDRutils::get_synthetic_data("small")
se <- mae[[gDRutils::get_supported_experiments("sa")]]

dt_treat <- gDRutils::convert_se_assay_to_dt(se, "RawTreated")
dt_controls <- gDRutils::convert_se_assay_to_dt(se, "Controls")

heatmap_control_mapping_qc(dt_treat = dt_treat,
                           dt_controls = dt_controls)


heatmap_control_mapping_qc(dt_treat = dt_treat[1:1350, ],
                           dt_controls = dt_controls)


dt_treat_NA <- dt_treat[-c(1:135, 270:405),]
heatmap_control_mapping_qc(dt_treat = dt_treat_NA,
                           dt_controls = dt_controls)


dt_controls_NA <- dt_controls[-c(1:305, 611:763, 1221:1750),]
heatmap_control_mapping_qc(dt_treat = dt_treat,
                           dt_controls = dt_controls_NA)


heatmap_control_mapping_qc(dt_treat = dt_treat,
                           dt_controls = dt_controls[1:3660, ])


```
