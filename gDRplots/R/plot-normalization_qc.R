#' Plot violin for normalized or averaged single-agent data to control quality of the data
#'
#' @param dt_assay data.table representing data from the assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, <assay_name>)}
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
  
  cl_clid <- unique(dt_assay[get(cellline_name) == cl_name, clid])
  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type, list(norm_type = normalization_type))
  dt_assay <- dt_assay[eval(filter_expr)]
  
  tab_subplot <- dt_assay[get(cellline_name) == cl_name, ]
  
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
