
#' Get Legend Title
#' 
#' prepare legend title formatted as HTML code
#'
#' @param var character name of the variable to be shown in legend
#' @param has_codrug_data logical whether the data is combo type
#' @param default_var character name of the variable to be shown in legend when \code{var} is not selected
#' 
#' @examples
#' get_legend_title("Drug")
#' get_legend_title("Drug 2", has_codrug_data = TRUE, default_var = "Conc_2")
#'
#' @keywords utils_legend
#' @return single element list with character in HTML format
#' 
#' @export
get_legend_title <-
  function(var,
           has_codrug_data = FALSE,
           default_var = gDRutils::get_prettified_identifiers("concentration2", simplify = TRUE)) {
    checkmate::assert_string(var, null.ok = TRUE)
    checkmate::assert_flag(has_codrug_data)
    checkmate::assert_string(default_var)
    
    if (is.null(var) || tolower(var) == "none") {
      if (has_codrug_data) {
        list(text = sprintf("<b>%s</b>", default_var))
      } else {
        NULL
      }
    } else {
      if (has_codrug_data && var != default_var) {
        list(text = sprintf("<b>%s</b> && <b>%s</b>",
                            default_var,
                            var))
      } else {
        list(text = sprintf("<b>%s</b>", var))
      }
    }
  }


#' Helper function to evaluate if legend should be added to the figure
#'
#' Current logic is as follows: show legend only if:
#' (1) user colors by "Drug MOA" or "Primary Tissue" or
#' (2) there is a codrug data with at least two unique values
#' 
#' @param var_col name of column used for coloring dots
#' @param data \code{data.table} prepared by \code{prepare_data_metric_contrast}
#' 
#' @examples
#' dt <- data.table::data.table(
#'   x = seq_len(10),
#'   y = seq_len(10) + 0.5,
#'   `Tissue` = c("lung", "brain")
#' )
#' do_show_legend("none", dt)
#' do_show_legend("Tissue", dt)
#' 
#' @keywords utils_legend
#' @return logical value
#' 
#' @export
do_show_legend <- function(var_col, 
                           data) {
  
  checkmate::assert_string(var_col, null.ok = FALSE)
  checkmate::assert_data_table(data)
  
  pidfs <- gDRutils::get_prettified_identifiers(simplify = TRUE)
  
  var_col %in% c(pidfs[["drug_moa"]], pidfs[["cellline_tissue"]]) ||
    length(unique(data[[pidfs[["concentration2"]]]])) >= 2
}
