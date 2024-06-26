#' Plot violin for normalized single-agent data to check quality of data
#'
#' @param dt_assay data.table representation of the data in assay
#'    output from \code{gDRutils::convert_se_assay_to_dt(se, <assya_name>)}
#' @param cl_name string cell line name to be plotted (Cell Line Name)
#' @param metric string with variable name to be plotted; it has to be in \code{dt_assay}
#' @param metric_growth string with normalization_types to be selected
#'                      one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
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
#'                          metric_growth = "RV")
#' 
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
#' plot_var_distribution_qc(dt_assay = dt_metrics,
#'                          cl_name = cl_name,
#'                          metric = "r2")
#'  
#' @export
plot_var_distribution_qc <- function(dt_assay,
                                     cl_name,
                                     metric = "x",
                                     metric_growth = "GR") {
  
  checkmate::expect_data_table(dt_assay)
  checkmate::expect_string(cl_name)
  checkmate::expect_choice(metric, choices = names(dt_assay))
  checkmate::expect_choice(metric_growth, choices = c("GR", "RV"))
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  
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
    ggplot2::labs(y = metric_growth, x = drug_name, title = plt_title) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))

  return(suppressWarnings(plt))
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
  
  # calculate frequency  
  frequency <- dt_controls[, .N, by = .(rId, cId)]
  result <- merge(unique(dt_treat[, c("rId", "cId")]), frequency, by = c("rId", "cId"), all.x = TRUE)
  data.table::setnames(result, "N", "Count")
  Count <- NULL # due to NSE notes in R CMD check
  
  # final plot
  plt <- ggplot2::ggplot(result, 
                         ggplot2::aes(x = cId, y = rId, fill = Count, colour = "")) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "#CEEFC8", high = "#76d364", na.value = "red") +
    ggplot2::scale_colour_manual(values = NA) +              
    ggplot2::guides(colour = ggplot2::guide_legend("No data", override.aes = list(fill = "red", colour = "black"))) +
    ggplot2::labs(title = "Mapping Counts Comparison between Treated and Control",
                  x = "col id",
                  y = "row id") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(hjust = 0)
    )
  
  return(plt)
}
