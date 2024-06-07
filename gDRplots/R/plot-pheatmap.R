#' Plot pretty heatmap with annotationsfor single-agent data
#'
#' @param tab_response \code{data.table} containing drug response metrics
#'    output from \code{\link[gDRutils]{convert_se_assay_to_dt}} for assay "Metrics" 
#'    and single-agent \code{SummarizedExperiment}
#' @param metric_growth string with normalization types to be selected
#'    one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param metric string name of metric;
#'    one of: "xc50"("GR50" or "IC50" - respectively depending on \code{metric_growth}), 
#'    "x_max" ("GR Max" or "E Max") or x_mean" ("GR Mean" or "RV Mean")
#' @param fit_source string source name for metrics
#' @param hm_title string plot title
#' @param colors_vec character vector of colors (valid name or hex) used in heatmap
#' @param no_breaks numeric number of breaks on scale
#' @param annotation_col \code{data.table} that specifies the annotations shown above the heatmap.
#'   Each row defines the features for a specific row. The rows in the data and in the annotation
#'   are matched using corresponding names from \code{CellLineName} column. 
#'   Note that color schemes takes into account if variable is continuous or discrete.
#' @param annotation_colors named list for specifying \code{annotation_col} track colors manually;
#'   note list is named wit annotation name (column names of \code{annotation_col} - 
#'   without \code{CellLineName}), each list item is named vector with valid colour name for 
#'   each value describe in \code{annotation_col}). Not described elements will be colored in default.
#' 
#' @seealso \code{\link[pheatmap]{pheatmap}}
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
#' annotation_map <- list(
#'   mut_A = c("1" = "coral", "0" = "cadetblue"),
#'   mut_B = c("yes" = "black", "no" = "grey90")
#' )
#' 
#' pheatmap_with_anno_sa(tab_response = response_metrics)
#' pheatmap_with_anno_sa(tab_response = response_metrics, 
#'                       metric_growth = "RV",
#'                       metric = "x_mean",
#'                       colors_vec = c("darkblue", "grey90"),
#'                       annotation_col = annotation_manual)
#' pheatmap_with_anno_sa(tab_response = response_metrics, 
#'                       annotation_col = annotation_manual,
#'                       annotation_colors = annotation_map,
#'                       hm_title = get_hm_title("Combo Matrix - single-agent data"))
#' 
#' @keywords QCS_plot
#' 
#' @return heatmap for selected metric with annotation - if given
#' @export 
pheatmap_with_anno_sa <- function(
    tab_response,
    metric_growth = "GR",
    metric = "xc50",
    fit_source = "gDR",
    hm_title = NA,
    colors_vec = c("firebrick2", "white"),
    no_breaks = 50,
    annotation_col = NULL,
    annotation_colors = NULL) {
  
  checkmate::assert_data_table(tab_response)
  checkmate::assert_choice(metric_growth, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, choices = c("x", "xc50", "x_max", "x_mean"))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_string(hm_title, na.ok = TRUE)
  checkmate::assert_character(colors_vec)
  stopifnot("Must be valid color name" = all(vapply(colors_vec, gDRplots::isValidColor, logical(1))))
  checkmate::assert_int(no_breaks, lower = 2)
  checkmate::assert_data_table(annotation_col, null.ok = TRUE)
  checkmate::assert_list(annotation_colors, null.ok = TRUE)
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  
  # prep data
  tab_response <- tab_response[normalization_type == metric_growth, ]
  
  if (fit_source %in% names(tab_response)) {
    data.table::setkeyv(tab_response, "fit_source")
    tab_response <- tab_response[fit_source]
    data.table::setkey(tab_response, NULL)
  }
  
  qmfun <- switch(metric,
                  "x" = identity,
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
  rm_col <- vapply(colnames(mat_cvd), function(i) !all(is.na(mat_cvd[, i])), logical(1))
  rm_row <- vapply(seq_along(rownames(mat_cvd)), function(i) !all(is.na(mat_cvd[i, ])), logical(1))
  if (!all(rm_col)) mat_cvd <- mat_cvd[, rm_col]
  if (!all(rm_row)) mat_cvd <- mat_cvd[rm_row, ]
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
      annotation_col[, (cols) := lapply(.SD, change_NA_into_char, "NA"), .SDcols = cols]
    }
    # select annotation acc to matrix
    annotation_col <- annotation_col[get(cellline_name) %in% rownames(mat_cvd), ]
    rownames(annotation_col) <- annotation_col[[cellline_name]] # required by pheatmap::pheatmap
    annotation_col <- annotation_col[, .SD, .SDcol = -cellline_name]
  }
  
  if (!is.null(annotation_col) && !is.null(annotation_colors)) {
    ls_ann_with_colors <- names(annotation_col)[names(annotation_col) %in% names(annotation_colors)]
    for (ann in ls_ann_with_colors) {
      reqired_lvl <- unique(annotation_col[[ann]])
      avaialable_lvl <- names(annotation_colors[[ann]])
      missing_lvl <- reqired_lvl[!reqired_lvl %in% avaialable_lvl]
      if (NROW(missing_lvl) == 1 && missing_lvl == "NA") {
        annotation_colors[[ann]] <- c(annotation_colors[[ann]], "NA" = "darkred")
      } else if (NROW(missing_lvl) > 0) {
        annotation_colors[[ann]] <- NULL # allow default colouring
      }
    }
  }
  
  # flip 
  t_mat_cvd <- t(mat_cvd)
  
  # prep hm color palette
  breaks <- seq(from = min(stats::na.omit(mat_cvd)), to = 1.0, length.out = no_breaks)
  hm_color_palette <- grDevices::colorRampPalette(colors_vec)(no_breaks + 1)
  
  hm <- pheatmap::pheatmap(t_mat_cvd,
                           scale = "none",
                           display_numbers = TRUE, 
                           number_color = "black", 
                           color = hm_color_palette,
                           breaks = breaks, 
                           angle_col = 90, 
                           main = hm_title,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           # manual annotation
                           annotation_col = annotation_col,
                           annotation_colors = annotation_colors
  )
  return(hm)
}


#' Plot pretty heatmap with annotations for combo data
#' 
#' @param tab_response \code{data.table} containing drug response metrics
#'    output from \code{\link[gDRutils]{convert_se_assay_to_dt}} for assay "scores" 
#'    and combo \code{SummarizedExperiment}
#' @param metric string name of combo metric;
#'    one of: "hsa_score"("Bliss Excess GR" or "Bliss Excess RV" - respectively depending on \code{metric_growth}), 
#'    "bliss_score" ("Bliss Score GR" or "Bliss Score RV")
#' @inheritParams pheatmap_with_anno_sa
#' 
#' @seealso \code{\link[pheatmap]{pheatmap}}
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' response_metrics <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "scores")
#' 
#' annotation_manual <- data.table::data.table(
#'   CellLineName = c("cellline_AA", "cellline_EA", "cellline_IB", "cellline_MC", "cellline_BC"),
#'   mut_A = c(1, 1, 1, 0, 0),
#'   mut_B = c("yes", "yes", "no", "no", "no")
#' )
#' 
#' annotation_map <- list(
#'   mut_A = c("1" = "coral", "0" = "cadetblue"),
#'   mut_B = c("yes" = "black", "no" = "grey90")
#' )
#' 
#' 
#' pheatmap_with_anno_combo(tab_response = response_metrics)
#' pheatmap_with_anno_combo(tab_response = response_metrics, 
#'                          metric_growth = "RV",
#'                          metric = "bliss_score",
#'                          colors_vec = c("darkblue", "grey90", "darkred"),
#'                          annotation_col = annotation_manual)
#' pheatmap_with_anno_combo(tab_response = response_metrics, 
#'                          annotation_col = annotation_manual,
#'                          annotation_colors = annotation_map,
#'                          hm_title = get_hm_title("Combo Matrix - combo data"))
#'             
#' @keywords QCS_plot
#' 
#' @return heatmap for selected metric with annotation - if given
#' @export 
pheatmap_with_anno_combo <- function(
    tab_response,
    metric_growth = "GR",
    metric = "hsa_score",
    fit_source = "gDR",
    hm_title = NA,
    colors_vec = c("royalblue3", "royalblue1", "grey95", "grey95", "firebrick1", "firebrick3"),
    no_breaks = 50,
    annotation_col = NULL,
    annotation_colors = NULL) {
  
  checkmate::assert_data_table(tab_response)
  checkmate::assert_choice(metric_growth, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, choices = c("hsa_score", "bliss_score"))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_string(hm_title, na.ok = TRUE)
  checkmate::assert_character(colors_vec)
  stopifnot("Must be valid color name" = all(vapply(colors_vec, gDRplots::isValidColor, logical(1))))
  checkmate::assert_int(no_breaks, lower = 2)
  checkmate::assert_data_table(annotation_col, null.ok = TRUE)
  checkmate::assert_list(annotation_colors, null.ok = TRUE)
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  
  # prep data
  tab_response <- tab_response[normalization_type == metric_growth, ]
  
  if (fit_source %in% names(tab_response)) {
    data.table::setkeyv(tab_response, "fit_source")
    tab_response <- tab_response[fit_source]
    data.table::setkey(tab_response, NULL)
  }
  
  # select data for normalization type
  tab_plot <- data.table::dcast(
    data = tab_response,
    formula = get(cellline_name) ~ paste(get(drug_name), "x", get(drug_name_2)),
    # TODO add DrugNamePlot
    value.var = metric)
  data.table::setnames(tab_plot, "cellline_name", cellline_name)
  
  # prep matrix
  mat_cvd <- as.matrix(tab_plot[, .SD, .SDcols = -cellline_name])
  rownames(mat_cvd) <- tab_plot[[cellline_name]]
  rm_col <- vapply(colnames(mat_cvd), function(i) !all(is.na(mat_cvd[, i])), logical(1))
  rm_row <- vapply(seq_along(rownames(mat_cvd)), function(i) !all(is.na(mat_cvd[i, ])), logical(1))
  if (!all(rm_col)) mat_cvd <- mat_cvd[, rm_col]
  if (!all(rm_row)) mat_cvd <- mat_cvd[rm_row, ]
  
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
      annotation_col[, (cols) := lapply(.SD, change_NA_into_char, "NA"), .SDcols = cols]
    }
    # select annotation acc to matrix
    annotation_col <- annotation_col[get(cellline_name) %in% rownames(mat_cvd), ]
    rownames(annotation_col) <- annotation_col[[cellline_name]] # required by pheatmap::pheatmap
    annotation_col <- annotation_col[, .SD, .SDcol = -cellline_name]
  }
  
  if (!is.null(annotation_col) && !is.null(annotation_colors)) {
    ls_ann_with_colors <- names(annotation_col)[names(annotation_col) %in% names(annotation_colors)]
    for (ann in ls_ann_with_colors) {
      reqired_lvl <- unique(annotation_col[[ann]])
      avaialable_lvl <- names(annotation_colors[[ann]])
      missing_lvl <- reqired_lvl[!reqired_lvl %in% avaialable_lvl]
      if (NROW(missing_lvl) == 1 && missing_lvl == "NA") {
        annotation_colors[[ann]] <- c(annotation_colors[[ann]], "NA" = "darkred")
      } else if (NROW(missing_lvl) > 0) {
        annotation_colors[[ann]] <- NULL # allow default colouring
      }
    }
  }
  
  # flip 
  t_mat_cvd <- t(mat_cvd)
  
  # prep hm color palette
  breaks <- seq(from = -0.7, to = 0.7, length.out = no_breaks)
  hm_color_palette <- grDevices::colorRampPalette(colors_vec)(no_breaks + 1)
  
  hm <- pheatmap::pheatmap(t_mat_cvd,
                           scale = "none",
                           display_numbers = TRUE,
                           number_color = "black", 
                           color = hm_color_palette,
                           breaks = breaks, 
                           angle_col = 90,
                           main = hm_title,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           # manual annotation
                           annotation_col = annotation_col,
                           annotation_colors = annotation_colors
  )
  return(hm)
}

# helpers ----
#' Get Legend Title
#' 
#' @param dataset_name string name of dataset
#' @param metric_growth string with normalization types to be selected
#'    one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param metric string name of metric
#'    one of: "xc50"("GR50" or "IC50" - respectively depending on \code{metric_growth}), 
#'    "x_max" ("GR Max" or "E Max") or x_mean" ("GR Mean" or "RV Mean")
#'
#' @examples
#' get_hm_title(dataset_name = "Dateset DX123",
#'              metric = "x_mean", 
#'              metric_growth = "GR")
#' 
#' @keywords QCS_plot
#' 
#' @return character title for heatmap
#' @export 
get_hm_title <- function(dataset_name,
                         metric = "xc50", 
                         metric_growth = "GR") {
  
  checkmate::assert_string(dataset_name)
  checkmate::assert_choice(metric_growth, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, 
                           choices = c("xc50", "x_max", "x_mean", "hsa_score", "bliss_score"))
  
  title_metric <- 
    gDRutils::prettify_flat_metrics(sprintf("%s_%s", metric, metric_growth), human_readable = TRUE)
  
  sprintf("%s (%s)", dataset_name, title_metric)
}

#' Change NA into given string
#'
#' @param x vector with items suspected of being NA
#' @param lbl_NA string - replacement for NA - as default "NA"
#'
#' @return character (for NA -> given string)
#' @keywords internal
#' 
change_NA_into_char <- function(x,
                                lbl_NA = "NA") {
  
  checkmate::assert_string(lbl_NA)
  
  ifelse(is.na(x), lbl_NA, as.character(x))
  
}
