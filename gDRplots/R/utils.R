#' Help function to change column suffix
#'
#' @param dt_ data.table with data
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' mae <- gDRutils::get_synthetic_data("small")
#' se <- mae[[1]]
#' dt_ <- gDRutils::convert_se_assay_to_dt(se, assay_name = "Metrics", wide_structure = TRUE)
#' 
#' swapped_dt_ <- swap_drugs_1_2(dt_)
#' }
#' 
#' @return swapped data.table
#' 
swap_drugs_1_2 <- function(dt_) {
  
  checkmate::assert_data_table(dt_)
  
  conc <- gDRutils::get_env_identifiers("concentration")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug <- gDRutils::get_env_identifiers("drug")
  drug_moa <- gDRutils::get_env_identifiers("drug_moa")

  for (dr_var in intersect(c(drug_name, drug, drug_moa, conc), colnames(dt_))) {
    dt_[[paste0("temp_", dr_var)]] <- dt_[[paste0(dr_var, "_2")]]
    dt_[[paste0(dr_var, "_2")]] <- dt_[[dr_var]] 
    dt_[[dr_var]] <- dt_[[paste0("temp_", dr_var)]]
    dt_[[paste0("temp_", dr_var)]] <- NULL
  }
  dt_
}


#' Estimate the optimal plot size (either ggplot or pheatmap) for saving plots
#'
#' @param plot a ggplot or pheatmap object
#' @param base_width an integer with default base_width
#' @param base_height an integer with default base_height
#' @param scale_factor an integer with default scale_factor
#'
#' @return named vector with optimal width and height used in ggsave function
#' @export
estimate_plot_size <- function(plot,
                               base_width = 10,
                               base_height = 6,
                               scale_factor = 0.5) {
  
  # Assert inputs
  checkmate::assert_multi_class(plot, c("ggplot", "pheatmap"))
  checkmate::assert_numeric(base_width, lower = 0, finite = TRUE)
  checkmate::assert_numeric(base_height, lower = 0, finite = TRUE)
  checkmate::assert_numeric(scale_factor, lower = 0, finite = TRUE)
  
  
  if (inherits(plot, "ggplot")) {
    # For ggplot2 objects
    plot_data <- ggplot2::ggplot_build(plot)$data
    num_elements <- length(unique(plot_data[[1]]$group))
    estimated_width <- base_width + num_elements * scale_factor
    estimated_height <- base_height + num_elements * scale_factor
  } else if (inherits(plot, "pheatmap")) {
    # For pheatmap objects
    matrix_position <- which(plot$gtable$layout$name == "matrix")
    matrix_dim <- dim(plot$gtable$grobs[[matrix_position]]$children[[1]]$gp$fill)
    num_rows <- matrix_dim[[1]]
    num_cols <- matrix_dim[[2]]
    estimated_width <- base_width + num_cols * scale_factor
    estimated_height <- base_height + num_rows * scale_factor
  } else {
    stop("Unsupported plot type. Only ggplot2 and pheatmap objects are supported.")
  }
  return(c(width = estimated_width, height = estimated_height))
}

