# Package index

## Single-agent plots

Function to plot drug response curves for single-agent data

- [`plot_boxplot_metric_sa()`](plot_boxplot_metric_sa.md) : Plot box
  plots for metric for single-agent data grouped by selected variable
- [`plot_boxplot_metric_sa_by_CLs()`](plot_boxplot_metric_sa_by_CLs.md)
  : Plot box plots for metric for single-agent data grouped by cell line
  names
- [`plot_boxplot_metric_sa_by_drugs()`](plot_boxplot_metric_sa_by_drugs.md)
  : Plot box plots for metric for single-agent data grouped by drug
  names
- [`plot_boxplot_metric_sa_by_grp()`](plot_boxplot_metric_sa_by_grp.md)
  : Plot box plots for metric for single-agent data grouped by selected
  variable
- [`plot_dose_response_sa()`](plot_dose_response_sa.md) : Plot drug
  response curves for single-agent data
- [`plot_dose_response_sa_by_CLs()`](plot_dose_response_sa_by_CLs.md) :
  Plot drug response curves for single-agent data for selected call
  lines and drugs
- [`plot_dose_response_sa_by_drugs()`](plot_dose_response_sa_by_drugs.md)
  : Plot drug response curves for single-agent data for selected call
  lines and drugs

## Combo plots

Function to plot averaged values heatmaps for combo data

- [`heatmap_combo_metrics()`](heatmap_combo_metrics.md) : Plot heatmaps
  of fitted values for combination metrics data
- [`heatmap_combo_metrics_panel()`](heatmap_combo_metrics_panel.md) :
  Plot panel of heatmaps of fitted values for combination metrics data
- [`heatmap_combo_with_isoref()`](heatmap_combo_with_isoref.md) : Plot
  heatmaps of averaged values for combination data
- [`heatmap_combo_with_isoref_panel()`](heatmap_combo_with_isoref_panel.md)
  : Plot panel of heatmaps with fitted and reference data for
  isobolograms
- [`heatmap_combo_with_isoref_panel_common()`](heatmap_combo_with_isoref_panel_common.md)
  : \#' Plot panel of heatmaps with fitted and reference data for
  isobolograms
- [`heatmap_combo_with_isoref_panel_independent()`](heatmap_combo_with_isoref_panel_independent.md)
  : Plot panel of heatmaps with fitted and reference data for
  isobolograms
- [`plot_boxplot_metric_combo()`](plot_boxplot_metric_combo.md) : Plot
  box plots for metric for combo data grouped by selected variable
- [`plot_boxplot_metric_combo_by_CLs()`](plot_boxplot_metric_combo_by_CLs.md)
  : Plot box plots for metric for combo data grouped by cell line names
- [`plot_boxplot_metric_combo_by_drugs()`](plot_boxplot_metric_combo_by_drugs.md)
  : Plot box plots for metric for combo data grouped by drug names
- [`plot_boxplot_metric_combo_by_grp()`](plot_boxplot_metric_combo_by_grp.md)
  : Plot box plots for metric for combo data grouped by selected
  variable
- [`plot_combination_index()`](plot_combination_index.md) : Plot line
  plot of combination index
- [`plot_dose_response_combo()`](plot_dose_response_combo.md) : Plot
  drug response curves for combo data
- [`plot_dose_response_combo_panel()`](plot_dose_response_combo_panel.md)
  : Plot panel with drug response curves for single-agent data to
  control quality of the data

## Quality control plots for Processing Report

QC plots for Processing Report with all helpers

- [`heatmap_control_mapping_qc()`](heatmap_control_mapping_qc.md) : Plot
  heatmap of mapping controls to treated for single-agent and combo data
  to control quality of the data
- [`pheatmap_qc()`](pheatmap_qc.md) : Plot pretty heatmap for
  single-agent or combo data to control quality of the data
- [`plot_dose_response_sa_qc()`](plot_dose_response_sa_qc.md) : Plot
  drug response curves for single-agent data to control quality of the
  data
- [`plot_dose_response_sa_qc_panel()`](plot_dose_response_sa_qc_panel.md)
  : Plot panel with drug response curves for single-agent data to
  control quality of the data
- [`plot_fitting_acc()`](plot_fitting_acc.md) : Visualization for the
  quality control of the fitting for single-agent data
- [`plot_plate()`](plot_plate.md) : Plot data from a specific column
- [`plot_plate_stack_info()`](plot_plate_stack_info.md) : Plot plate
  data (Heatmap + Dose Ranks + Smart Legend + Explicit QC)
- [`plot_single_plate_stack_info()`](plot_single_plate_stack_info.md) :
  Plot a single plate's stack info
- [`plot_var_distribution_qc()`](plot_var_distribution_qc.md) : Plot a
  violin plot for normalized or averaged data (single-agent and combo)
  to control the quality of the data
- [`plot_var_stat_qc()`](plot_var_stat_qc.md) : Lollipop plot for metric
  single-agent data to control quality of the data

## Pretty heatmap with annotation

Functions to plot annotated heatmap with helpers

- [`get_hm_title()`](get_hm_title.md) : Get Legend Title
- [`pheatmap_with_anno_cd()`](pheatmap_with_anno_cd.md) : Plot pretty
  heatmap with annotations for co-dilution data
- [`pheatmap_with_anno_combo()`](pheatmap_with_anno_combo.md) : Plot
  pretty heatmap with annotations for combo data
- [`pheatmap_with_anno_combo_metrics()`](pheatmap_with_anno_combo_metrics.md)
  : Plot heatmap with annotations for combo metrics data
- [`pheatmap_with_anno_sa()`](pheatmap_with_anno_sa.md) : Plot pretty
  heatmap with annotations for single-agent data
- [`prep_pheatmap_matrix()`](prep_pheatmap_matrix.md) : Prep matrix with
  metric value based on the Metrics assay

## PRISM visualisation

Functions to plot PRISM visualisation with helpers

- [`create_PRISM_plot_list_combo()`](create_PRISM_plot_list_combo.md) :
  Create a nested list of plots for PRISM data with combo metrics
- [`create_PRISM_plot_list_sa()`](create_PRISM_plot_list_sa.md) : Create
  a nested list of plots for PRISM data with single-agent metrics
- [`create_PRISM_summary_list()`](create_PRISM_summary_list.md) : Create
  a list of PRISM association table
- [`.create_PRISM_plot_list()`](dot-create_PRISM_plot_list.md) : Create
  a list of PRISM plots for a selected type of experiment
- [`.get_data_type()`](dot-get_data_type.md) : Check data type
- [`.get_n_top_asssoc()`](dot-get_n_top_asssoc.md) : Get n-top linear
  associations
- [`plot_boxplot_meta()`](plot_boxplot_meta.md) : Plot boxplot for
  metric values grouped by metadata from DepMap
- [`plot_boxplot_num()`](plot_boxplot_num.md) : Plot boxplot for
  categorical features
- [`plot_boxplot_num_panel()`](plot_boxplot_num_panel.md) : Plot panel
  with boxplots
- [`plot_scatter_with_corr()`](plot_scatter_with_corr.md) : Plot scatter
  with correlation
- [`plot_scatter_with_corr_panel()`](plot_scatter_with_corr_panel.md) :
  Plot panel with scatter with correlation
- [`plot_volcano_assoc()`](plot_volcano_assoc.md) : Volcano plot with
  association
- [`plot_volcano_assoc_panel()`](plot_volcano_assoc_panel.md) : Plot
  panel with volcano plot and according to the data type - scatter plots
  or box plots
- [`prep_dt_assoc()`](prep_dt_assoc.md) : Prep table with calculated
  linear associations
- [`prep_dt_response_dose_sa()`](prep_dt_response_dose_sa.md) : Prep
  table with metric values by doses for single-agent experiment
- [`prep_dt_response_metric_diff()`](prep_dt_response_metric_diff.md) :
  Prep table with metric values for combination experiment
- [`prep_dt_response_metric_sa()`](prep_dt_response_metric_sa.md) : Prep
  table with metric values for single-agent experiment
- [`prep_dt_response_scores()`](prep_dt_response_scores.md) : Prep table
  with metric values for combination experiment

## Color helpers

Plot color helpers

- [`fill_ann_color_map()`](fill_ann_color_map.md) : Fill missing values
  in the color map for annotation
- [`get_ann_color_map()`](get_ann_color_map.md) : Create color map for
  annotation
- [`get_col_luminance()`](get_col_luminance.md) : Calculate the
  luminance of a color
- [`get_iso_colors()`](get_iso_colors.md) : get_iso_colors
- [`get_qual_colors()`](get_qual_colors.md) : Create list of qualitative
  colors
- [`is_color_dark()`](is_color_dark.md) : Determine whether or not a
  color is dark
- [`is_valid_color()`](is_valid_color.md) : Determine whether or not a
  color name is valid

## Utils

Utils function

- [`compute_distances()`](compute_distances.md) : Compute distance
  between rows of a matrix.
- [`create_log_seq()`](create_log_seq.md) : create a log-sequence
- [`round_to_unique_string()`](round_to_unique_string.md) : Round
  numbers to unique string

## Chemical genomics visualisation

Functions to plot chemical genomics visualisation with helpers

- [`analyze_cgs()`](analyze_cgs.md) : Analyze Chemical Genomics (CGS)
  Data and Perform GSEA
- [`plot_cgs_ranking()`](plot_cgs_ranking.md) : Plot Chemical Genomics
  Screen GSEA Results
