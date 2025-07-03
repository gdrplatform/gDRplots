## gDRplots 0.0.86 - 2025-07-03
* refactor `create_PRISM_plot_list_sa` and `create_PRISM_plot_list_combo`

## gDRplots 0.0.85 - 2025-07-02
* add wrapper for `DT::datatable`

## gDRplots 0.0.84 - 2025-06-25
* test data for PRISM
* update and ad dmissing tests

## gDRplots 0.0.83 - 2025-06-24
* clean gene names in PRISM reports

## gDRplots 0.0.82 - 2025-06-06
* update `pheatmap_qc` with standardized concentrations

## gDRplots 0.0.81 - 2025-05-28
* update arg names in cgs functions

## gDRplots 0.0.80 - 2025-05-27
* update chunk naming to fulfill `knitr` requirements

## gDRplots 0.0.79 - 2025-05-05
* update `prep_dt_depmap_feat`

## gDRplots 0.0.78 - 2025-04-29
* add support for cgs analysis for difference between two cell lines

## gDRplots 0.0.77 - 2025-04-15
* add function to calculate linear associations

## gDRplots 0.0.76 - 2025-04-07
* update `prep_dt_depmap_meta`

## gDRplots 0.0.75 - 2025-03-27
* add link param in `prep_plot_chunk`
* add download feature in `prep_plot_chunk` and `prep_double_table_chunk`

## gDRplots 0.0.74 - 2025-03-26
* handle edge-case in `pheatmap_with_anno_combo`

## gDRplots 0.0.73 - 2025-03-25
* fix readability of number for PRISM boxplots

## gDRplots 0.0.72 - 2025-03-24
* change xc50 to log10 for PRISM analysis

## gDRplots 0.0.71 - 2025-03-12
* add coloring points in metric boxplots

## gDRplots 0.0.70 - 2025-03-10
* update xc50 capping logic for combo PRISM inputs

## gDRplots 0.0.69 - 2025-03-03
* protect against unexpected data aggregation in `data.table::dcast`

## gDRplots 0.0.68 - 2025-02-25
* move `is_color_dark` from gDRiPlots package
* add `.get_pheatmap_number_color` function

## gDRplots 0.0.67 - 2025-02-17
* add `prep_pheatmap_matrix` function
* add `dt_metrics_capped` param to `pheatmap_with_anno_sa` function
* update xc50 capping logic for PRISM inputs

## gDRplots 0.0.66 - 2024-02-15
* update `ls_selected_met` for PRISM data

## gDRplots 0.0.65 - 2024-02-12
* add support for chemical genomics analysis

## gDRplots 0.0.64 - 2025-02-12
* downgrade reqs for Bioc pkgs

## gDRplots 0.0.63 - 2025-02-10
* refactor `plot_dose_response_sa` to handle missing data consistently
* fix labeling in `pheatmap_qc`
* refactor panel shape for `heatmap_combo_metrics_panel`

## gDRplots 0.0.62 - 2025-02-06
* update `plot_boxplot_meta` to handle edge cases 
* add `plot_boxplot_num` and `plot_boxplot_num_panel` functions
* add new logic for selecting n-top based on q-val and rho

## gDRplots 0.0.61 - 2025-02-05
* add values to combo heatmaps and adjust position of plots in the panel

## gDRplots 0.0.60 - 2025-01-31
* fix y-axis range for dose-response curves

## gDRplots 0.0.59 - 2025-01-27
* update `.get_data_type` with check for relation number

## gDRplots 0.0.58 - 2025-01-24
* update handling special characters in rmd (double colon escape)

## gDRplots 0.0.57 - 2025-01-13
* move `compute_distances` from gDRiPlots package
* turn on clustering in the heatmaps with infinite and NA values

## gDRplots 0.0.56 - 2022-01-13
* move ggrepel from Suggests to Imports

## gDRplots 0.0.55 - 2025-01-03
* add `get_r_file_path` function

## gDRplots 0.0.54 - 2024-12-20
* refactor `prep_plot_chunk` to support simple nested lists of visualizations
* refactor plate visualizations to improve the performance

## gDRplots 0.0.53 - 2024-12-10
* refactor handling 0 in concentration column and -Inf and Inf for metrics

## gDRplots 0.0.52 - 2024-12-09
* order subplots in `plot_scatter_with_corr_panel` according to given list

## gDRplots 0.0.51 - 2024-12-04
* revert color palette for dose-response curves for combo data

## gDRplots 0.0.50 - 2024-12-03
* fix for edge-case in `pheatmap_with_anno_sa` (-Inf)

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
