#' Plot QCS heatmap
#'
#' @param tab_response \code{data.table} containing drug response metrics
#'    output from \code{\link[gDRutils]{convert_se_assay_to_dt}}
#' @param metric_growth string with normalization_types to be selected
#'    one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param metric string name of metric
#'    one of: "xc50"("GR50" or "IC50" - respectively depending on \code{metric_growth}), 
#'    "x_max" ("GR Max" or "E Max") or x_mean" ("GR Mean" or "RV Mean")
#' @param fit_source fit source name for new GDS metrics
#' @param dataset_name string name of dataset to be shown in t
#' @param mapcolor character representing a valid colour - name or hex
#' @param no_breaks numeric number of breaks on scale
#' @param annotation_col data.table that specifies the annotations shown above the heatmap.
#'   Each row defines the features for a specific row. The rows in the data and in the annotation
#'   are matched using corresponding names from \code{CellLineName} column. 
#'   Note that color schemes takes into account if variable is continuous or discrete.
#' @param annotation_colors named list for specifying \code{annotation_col} track colors manually;
#'   note list is named wit annotation name (columne names of \code{annotation_col} - 
#'   without \code{CellLineName}), each list item is named vector with valid colour name for 
#'   each value describe in \code{annotation_col}). Not described elements will be colored in default.
#'   more detail see \code{\link[pheatmap]{pheatmap}}
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#' response_metrics <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "Metrics")
#' 
#' annotation_manual <- data.table::data.table(
#'   CellLineName = c("cellline_AA", "cellline_EA", "cellline_IB", "cellline_MC", "cellline_BC"),
#'   mut_A = c(1, 1, 1, 0, 0),
#'   mut_B = c("yes", "yes", "no", "no", "no")
#' )
#' 
#' heatmap_QCS(tab_response = response_metrics)
#' heatmap_QCS(tab_response = response_metrics, annotation_col = annotation_manual)
#' 
#' @return heatmap for selected metric wit annotation - if given
#' @export 
heatmap_QCS <- function(tab_response,
                        metric_growth = "GR",
                        metric = "xc50",
                        fit_source = "gDR",
                        dataset_name = NULL, # TODO change int plot_title
                        mapcolor = c("firebrick2", "white"),
                        no_breaks = 50,
                        annotation_col = NULL,
                        annotation_colors = NULL) {
  
  checkmate::assert_data_table(tab_response)
  checkmate::assert_choice(metric_growth, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, choices = c("xc50", "x_max", "x_mean"))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_string(dataset_name, null.ok = TRUE)
  checkmate::assert_character(mapcolor)
  stopifnot("Must be valid color name" = all(vapply(mapcolor, gDRplots::isValidColor, logical(1))))
  checkmate::assert_int(no_breaks, lower = 2)
  checkmate::assert_data_table(annotation_col, null.ok = TRUE)
  checkmate::assert_list(annotation_colors, null.ok = TRUE)
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  
  # prep data
  tab_response <- tab_response[normalization_type == metric_growth,]
  
  if (fit_source %in% names(tab_response)) {
    data.table::setkeyv(tab_response, "fit_source")
    tab_response <- tab_response[fit_source]
    data.table::setkey(tab_response, NULL)
  }
  
  qmfun <- switch(metric,
                  "xc50" = log10, 
                  "x_max" = identity, 
                  "x_mean" = identity)
  
  # select data for normalization type
  tab_plot <- data.table::dcast(
    data = tab_response,
    formula = get(cellline_name) ~  get(drug_name),
    # TODO add DrugNamePlot
    value.var = metric)
  data.table::setnames(tab_plot, "cellline_name", cellline_name)
  
  # prep matrix
  mat_cvd <- as.matrix(tab_plot[, .SD, .SDcols = -cellline_name])
  rownames(mat_cvd) <- tab_plot[[cellline_name]]
  rm_col <- vapply(colnames(mat_cvd), function(i) !all(is.na(mat_cvd[,i])), logical(1))
  rm_row <- vapply(seq_along(rownames(mat_cvd)), function(i) !all(is.na(mat_cvd[i,])), logical(1))
  if (!all(rm_col)) mat_cvd <- mat_cvd[,rm_col]
  if (!all(rm_row)) mat_cvd <- mat_cvd[rm_row,]
  mat_cvd[] <- vapply(mat_cvd, function(x) qmfun(x), numeric(1))
  
  # check completeness of annotation - TODO wrap in separate function
  if (!is.null(annotation_col)) {
    if (!all(rownames(mat_cvd) %in% annotation_col[[cellline_name]])) {
      tab_missing_ann <- data.table::data.table(
        missing = rownames(mat_cvd)[!rownames(mat_cvd) %in% annotation_col[[cellline_name]]]
      )
      data.table::setnames(tab_missing_ann, "missing", cellline_name)
      
      annotation_col <- data.table::rbindlist(list(annotation_col, tab_missing_ann), fill = TRUE)
      # data.table::nafill does not support character
      cols <- names(annotation_col)[names(annotation_col) != cellline_name]
      annotation_col[ , (cols) := lapply(.SD, change_NA_into_char), .SDcols = cols]
    }
  }
  
  if (!is.null(annotation_col) && !is.null(annotation_colors)) {
    # TODO
  }
  
  # flip 
  t_mat_cvd <- t(mat_cvd)
  
  # prep hmcolors
  breaks <- seq(from = min(na.omit(mat_cvd)), to = 1.0, length.out = no_breaks)
  hmcol <- grDevices::colorRampPalette(mapcolor)(no_breaks + 1)
  
  hm <- pheatmap::pheatmap(t_mat_cvd,
                           scale = "none",
                           display_numbers = TRUE, 
                           fontsize_number = 6,
                           number_color = "black", 
                           color = hmcol,
                           breaks = breaks, 
                           angle_col = 90, 
                           fontsize = 6,
                           treeheight_row = 30, 
                           treeheight_col = 30,
                           main = paste(metric, dataset_name),
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           # manual annotation
                           annotation_col = annotation_col,
                           annotation_colors = annotation_colors
  )
  return(hm)
}


# helpers ----
#' Change NA into given string
#'
#' @param x vector with items suspected of being NA
#' @param lbl_NA string - replacement for NA - as default "NA"
#'
#' @return character (for NA -> given string)
#' @keywords internal
#' @md

change_NA_into_char <- function(x,
                                lbl_NA = "NA") {
  
  checkmate::assert_string(lbl_NA)
  
  ifelse(is.na(x), lbl_NA, as.character(x))
  
}

