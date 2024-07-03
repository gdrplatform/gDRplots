#' Plot violin for normalized or averaged single-agent data to check quality of data
#'
#' @param dt_assay data.table representation of the data in assay
#'    output from \code{gDRutils::convert_se_assay_to_dt(se, <assay_name>)} 
#'    for assay_name like "Normalized" and "Averaged"
#' @param cl_name string cell line name to be plotted (Cell Line Name)
#' @param metric string with variable name to be plotted; it has to be in \code{dt_assay}
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#'
#' @return plot with violin for each drug
#'
#' @keywords QC_plot
#' @examples
#' mae <- gDRutils::get_synthetic_data("small")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#' 
#' dt_norm <- gDRutils::convert_se_assay_to_dt(se, "Normalized")
#' cl_name <- dt_norm[["CellLineName"]][1]
#' 
#' plot_var_distribution_qc(dt_assay = dt_norm,
#'                          cl_name = cl_name)
#' 
#' dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
#' plot_var_distribution_qc(dt_assay = dt_average,
#'                          cl_name = cl_name,
#'                          normalization_type = "RV")
#'                          
#' plot_var_distribution_qc(dt_assay = dt_average,
#'                          cl_name = cl_name,
#'                          metric = "x_std",
#'                          normalization_type = "RV")
#'                          
#' @export
plot_var_distribution_qc <- function(dt_assay,
                                     cl_name,
                                     metric = "x",
                                     normalization_type = "GR") {
  
  checkmate::expect_data_table(dt_assay)
  checkmate::expect_string(cl_name)
  checkmate::expect_choice(metric, choices = names(dt_assay))
  checkmate::expect_choice(normalization_type, choices = c("GR", "RV"))
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  
  # ----
  cl_clid <- unique(dt_assay[get(cellline_name) == cl_name, ][[clid]]) 
  tab_subplot <- dt_assay[normalization_type == metric_growth & get(cellline_name) == cl_name, ]
  
  plt_title <- sprintf("%s (%s)", cl_name, cl_clid)
  color_palette <- get_qual_colors(NROW(unique(tab_subplot[[drug_name]])))
  
  plt <- ggplot2::ggplot(tab_subplot, ggplot2::aes(x = get(drug_name), y = !!rlang::sym(metric))) +
    ggplot2::geom_hline(yintercept = c(0, 1), color = "#2c3e50", linetype = "dashed") +
    ggplot2::geom_violin(ggplot2::aes(fill = get(drug_name), color = get(drug_name)), 
                         alpha = 0.25, na.rm = TRUE, drop = FALSE) +
    ggplot2::geom_jitter(width = 0.2, height = 0, color = "#2c3e50") +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = color_palette) +
    ggplot2::scale_color_manual(values = color_palette) +
    ggplot2::labs(y = sprintf("%s for %s", metric, normalization_type), x = drug_name, title = plt_title) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))
  
  return(plt)
}

#' Lollipop plot for metric single-agent data to check quality of data
#'
#' @param dt_assay data.table representation of the data in assay
#'    output from \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#' @param cl_name string cell line name to be plotted (Cell Line Name)
#' @param metric string with variable name to be plotted; it has to be in \code{dt_assay}
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param with_table logical whether table with metric values should be shown next to the plot
#'
#' @return lollipop plot with stat value for each drug
#'
#' @keywords QC_plot
#' @examples
#' mae <- gDRutils::get_synthetic_data("small")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#' 
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
#' cl_name <- dt_metrics[["CellLineName"]][1]
#' 
#' plot_var_stat_qc(dt_assay = dt_metrics,
#'                  cl_name = cl_name)
#'                          
#' plot_var_stat_qc(dt_assay = dt_metrics,
#'                  cl_name = cl_name,
#'                  metric = "r2",
#'                  normalization_type = "RV")
#'                          
#' plot_var_stat_qc(dt_assay = dt_metrics,
#'                  cl_name = cl_name,
#'                  metric = "x_AOC",
#'                  normalization_type = "RV",
#'                  with_table = TRUE)
#'                          
#' @export
plot_var_stat_qc <- function(dt_assay,
                             cl_name,
                             metric = "r2",
                             normalization_type = "GR", 
                             with_table = FALSE) {
  
  checkmate::expect_data_table(dt_assay)
  checkmate::expect_string(cl_name)
  checkmate::expect_choice(metric, choices = names(dt_assay))
  checkmate::expect_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::expect_flag(with_table)
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  
  # --- 
  cl_clid <- unique(dt_assay[get(cellline_name) == cl_name, clid]) 
  tab_subplot <- dt_assay[normalization_type == metric_growth & get(cellline_name) == cl_name, ]
  
  plt_title <- sprintf("%s (%s)", cl_name, cl_clid)
  color_palette <- get_qual_colors(NROW(unique(tab_subplot[[drug_name]])))
  
  plt <- ggplot2::ggplot(tab_subplot, ggplot2::aes(x = get(drug_name), y = !!rlang::sym(metric))) +
    ggplot2::geom_hline(yintercept = c(0, 1), color = "#2c3e50", linetype = "dashed") +
    ggplot2::geom_segment(
      ggplot2::aes(x = get(drug_name), xend = get(drug_name), y = 0, yend = !!rlang::sym(metric))) +
    ggplot2::geom_point(ggplot2::aes(fill = get(drug_name), color = get(drug_name)), 
                        alpha = 0.75, size = 5, shape = 21, stroke = 1) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = color_palette) +
    ggplot2::scale_color_manual(values = color_palette) +
    ggplot2::labs(y = sprintf("%s for %s", metric, normalization_type), x = drug_name, title = plt_title) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))
  
  if (with_table) {
    tab_metric <- ggpubr::ggtexttable(
      tab_subplot[, .SD, .SDcols = c(drug_name, metric)][order(get(metric))], 
      rows = NULL, theme = ggpubr::ttheme("light")) 
    
    plt <- ggpubr::ggarrange(plt, tab_metric, nrow = 1, widths = c(2, 1))
  }
  
  return(plt)
}

#' Plot drug response curves for single-agent data to check quality of data
#'
#' @param dt_treat data.table representation of the data in \code{RawTreated} assay
#'    output from \code{gDRutils::convert_se_assay_to_dt(se, "RawTreated")}
#' @param dt_controls data.table representation of the data in \code{Controls} assay
#'    output from \code{gDRutils::convert_se_assay_to_dt(se, "Controls")}
#'
#' @return hetamap with check of mapping controls to treated
#'
#' @keywords QC_plot
#' @examples
#' mae <- gDRutils::get_synthetic_data("small")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#' 
#' dt_treat <- gDRutils::convert_se_assay_to_dt(se, "RawTreated")
#' dt_controls <- gDRutils::convert_se_assay_to_dt(se, "Controls")
#' 
#' heatmap_control_mapping_qc(dt_treat = dt_treat,
#'                            dt_controls = dt_controls)
#' 
#' heatmap_control_mapping_qc(dt_treat = dt_treat[1:1350, ],
#'                            dt_controls = dt_controls)
#' 
#' dt_treat_NA <- dt_treat[-c(1:135, 270:405),]
#' heatmap_control_mapping_qc(dt_treat = dt_treat_NA,
#'                            dt_controls = dt_controls)
#' 
#' dt_controls_NA <- dt_controls[-c(1:305, 611:763, 1221:1750),]
#' heatmap_control_mapping_qc(dt_treat = dt_treat,
#'                            dt_controls = dt_controls_NA)
#' 
#' heatmap_control_mapping_qc(dt_treat = dt_treat,
#'                            dt_controls = dt_controls[1:3660, ])
#' 
#' 
#' @export
heatmap_control_mapping_qc <- function(dt_treat,
                                       dt_controls) {
  checkmate::expect_data_table(dt_treat)
  checkmate::expect_data_table(dt_controls)
  
  # calculate the frequency of each (rID, cID) combination in Controls 
  frequency <- dt_controls[, .N, by = .(rId, cId)]
  # merge the frequency with the Treated data.table
  result <- merge(unique(dt_treat[, c("rId", "cId")]), frequency, by = c("rId", "cId"), all.x = TRUE)
  
  # Convert the result to a matrix format suitable for pheatmap
  result_matrix <- data.table::dcast(result, rId ~ cId, value.var = "N")
  rownames <- result_matrix$rId
  result_matrix <- as.matrix(result_matrix[, !("rId"), with = FALSE])
  rownames(result_matrix) <- rownames
  
  # Replace 0 with NA to use na_col for red color
  result_matrix[result_matrix == 0] <- NA
  
  # Generate the breaks for integers
  max_count <- max(result_matrix, na.rm = TRUE)
  breaks <- seq(1, max_count, by = 1)
  unique_values <- unique(stats::na.omit(as.vector(result_matrix)))
  
  # Generate the heatmap
  pheatmap::pheatmap(result_matrix,
                     color = grDevices::colorRampPalette(c("#CEEFC8", "#76d364"))(length(breaks) - 1),
                     breaks = breaks,
                     na_col = "red",
                     main = "Counts of mapped controls",
                     cluster_cols = FALSE,
                     cluster_rows = FALSE,
                     angle_col = 90,
                     legend_breaks = unique_values,
                     legend_labels = unique_values)
}
