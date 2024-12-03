## gDRplots 0.0.50 - 2024-12-03
* revert color palette for dose-response curves for combo data

## gDRplots 0.0.49 - 2024-12-02
* add boxplots for single-agent and combo metrics

## gDRplots 0.0.48 - 2024-11-28
* add plate visualizations

## gDRplots 0.0.47 - 2024-11-26
* refactor `heatmap_combo_metrics_panel` function - split in single plot and panel
* add `plot_combination_index` function

## gDRplots 0.0.46 - 2024-11-25
* update `plot_volcano_assoc_panel`
* unify function name - rename `heatmap_combo_metrics` into `heatmap_combo_metrics_panel`

## gDRplots 0.0.45 - 2024-11-20
* add possibility for plotting isoline references for HSA Excess and Bliss Excess

## gDRplots 0.0.44 - 2024-11-18
* add possibility for switching axes in the combo heatmaps

## gDRplots 0.0.43 - 2024-11-14
* add symmetric limits on the heatmaps for Bliss and HSA excesses

## gDRplots 0.0.42 - 2024-11-04
* fix invalid title for dose-response curves for single-agent

## gDRplots 0.0.41 - 2024-10-30
* update handling special characters in rmd

## gDRplots 0.0.40 - 2024-10-23
* update settings.json with colors palette

## gDRplots 0.0.39 - 2024-10-22
* add support for visualizing combo data without isobolograms

## gDRplots 0.0.38 - 2024-10-04
* add helper functions for knitting Rmarkdown reports

## gDRplots 0.0.37 - 2024-10-02
* update vignette

## gDRplots 0.0.36 - 2024-09-25
* add plots and help functions for PRISM data 

## gDRplots 0.0.35 - 2024-09-23
* update clustering condition
* add `pheatmap_with_anno_cd` function
 
## gDRplots 0.0.34 - 2024-09-16
* add clustering in `pheatmap_with_anno_sa` and `pheatmap_with_anno_combo` functions

## gDRplots 0.0.33 - 2024-09-10
* increase the readability of QC dose-response curves for single-agent

## gDRplots 0.0.32 - 2024-09-09
* fix infinity values on the x-axis

## gDRplots 0.0.31 - 2024-09-03
* add vignette

## gDRplots 0.0.30 - 2024-08-28
* replace viridis with another set of colors

## gDRplots 0.0.29 - 2024-08-23
* fix bug in plot axis title
* add horizontal line in 1 for dose response curve

## gDRplots 0.0.28 - 2024-08-22
* update heatmap function to return plot and data - change file schema

## gDRplots 0.0.27 - 2024-08-20
* move functions related to interactive plots to the `gDRiPlots` package

## gDRplots 0.0.26 - 2024-08-20
* make plotting functions more versatile

## gDRplots 0.0.25 - 2024-08-13
* increase readability for reference isobolograms for combo data

## gDRplots 0.0.24 - 2024-08-12
* fix issue with coloring dose-response curve for combo data 

## gDRplots 0.0.23 - 2024-08-09
* add function for neutralizing spaces in file names
* add tests for rmarkdown helpers

## gDRplots 0.0.22 - 2024-08-08
* update heatmap function to return plot and data

## gDRplots 0.0.21 - 2024-08-05
* update plots for quality control of the distribution of normalized combo data

## gDRplots 0.0.20 - 2024-08-05
* extend the range of x-axis in dose-response curve

## gDRplots 0.0.19 - 2024-08-01
* prevent drugCombo plot from overwriting the source data

## gDRplots 0.0.18 - 2024-07-26
* improve the ratio of square plots
* unify text size across qc plots

## gDRplots 0.0.17 - 2024-07-25
* add a function for saving plots

## gDRplots 0.0.16 - 2024-07-24
* fix an issue with annotations coloring due to multiple occurrences
* fix the lack of color for dose-response curves for single-agent in a small set of data
* refactor naming for control_mapping heatmaps

## gDRplots 0.0.15 - 2024-07-16
* fix issue with wrong order of colorbar values for combo heatmap

## gDRplots 0.0.14 - 2024-07-15
* refactor code to properly visualize data with NAs

## gDRplots 0.0.13 - 2024-07-15
* fix issue with app crash when limits in `create_color_palette` are the same

## gDRplots 0.0.12 - 2024-07-12
* move `get_combo_col_settings` and `get_iso_colors` from `gDRutils` package

## gDRplots 0.0.11 - 2024-07-11
* add plots to quality control for combo data - combo metrics and reference isoline

## gDRplots 0.0.10 - 2024-07-10
* add plot to visualize fitting accuracy

## gDRplots 0.0.9 - 2024-07-03
* add plots to quality control of distribution of normalized data and metric stat

## gDRplots 0.0.8 - 2024-07-01
* add function for plot size estimation

## gDRplots 0.0.7 - 2024-06-21
* add plot to quality control of dose response

## gDRplots 0.0.6 - 2024-06-19
* update plot function for `processingReport` package

## gDRplots 0.0.5 - 2024-06-14
* move plot function from `metricClustering`, `responseGrid` and `metricDistribution` packages
* move `get_metrics_to_transform`, `convert_factor_to_character` and `paletteBrew` from `gDRcomponents` package
* rename function in snake case

## gDRplots 0.0.4 - 2024-06-07
* update plot functions for `processingReport` package

## gDRplots 0.0.3 - 2024-05-27
* add plot functions for `processingReport` package

## gDRplots 0.0.2 - 2024-05-20
* move plot functions from `gDRcomponents` package

## gDRplots 0.0.1 - 2024-05-15
* update package structure
* update SA plots and combo plots functions
* add basic tests
