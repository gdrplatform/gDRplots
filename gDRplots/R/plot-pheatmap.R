#' Plot pretty heatmap for single-agent or combo data to check quality of data
#'
#' @param tab_response \code{data.table} containing drug response metrics
#'    output from \code{\link[gDRutils]{convert_se_assay_to_dt}} for assay "Averaged" 
#'    and \code{SummarizedExperiment} with chosen data type: single-agent or combo
#' @param metric_growth string with normalization types to be selected
#'    one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param metric string name of metric;
#'    one of: "x" (value of "GR" or "RV" itself - respectively depending on \code{metric_growth}), 
#'    or "x_std" (standard deviation)
#' @param fit_source string source name for metrics
#' @param hm_title string plot title
#' @param colors_vec character vector of colors (valid name or hex) used in heatmap 
#'   (first colour for min value, last colour - for max value)
#' @param no_breaks numeric number of breaks on scale
#' @param cluster_rows logical flag whether ows should be clustered
#' @param lbl_by_CellLineName logical flag whether heatmap should be described by CellLineNames instead of clid
#' @param lbl_by_DrugName logical flag whether heatmap should be described by DrugName instead of Gnumber
#' 
#' @seealso \code{\link[pheatmap]{pheatmap}}
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]][2:5, ]
#' response_metrics <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "Averaged")
#' 
#' pheatmap_qc(tab_response = response_metrics)
#' pheatmap_qc(tab_response = response_metrics, 
#'             metric_growth = "RV",
#'             colors_vec = c("darkblue", "grey90"),
#'             lbl_by_CellLineName = TRUE,
#'             lbl_by_DrugName = TRUE)
#'              
#'              
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' response_metrics <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "Averaged")
#' pheatmap_qc(tab_response = response_metrics,
#'            cluster_rows = FALSE)
#'              
#' pheatmap_qc(tab_response = response_metrics,
#'             metric = "x_std",
#'             cluster_rows = FALSE)
#' 
#' @keywords QC_plot
#' 
#' @return heatmap for selected metric with annotation - if given
#' @export 
pheatmap_qc <- function(
    tab_response,
    metric_growth = "GR",
    metric = "x",
    fit_source = "gDR",
    hm_title = NA,
    colors_vec = c("black", "grey99"),
    no_breaks = 50,
    cluster_rows = TRUE,
    lbl_by_CellLineName = FALSE,
    lbl_by_DrugName = FALSE) {
  
  checkmate::assert_data_table(tab_response)
  checkmate::assert_choice(metric_growth, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, choices = c("x", "x_std"))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_string(hm_title, na.ok = TRUE)
  checkmate::assert_character(colors_vec)
  stopifnot("Must be valid color name" = all(vapply(colors_vec, gDRplots::is_valid_color, logical(1))))
  checkmate::assert_int(no_breaks, lower = 2)
  checkmate::assert_flag(cluster_rows)
  checkmate::assert_flag(lbl_by_CellLineName)
  checkmate::assert_flag(lbl_by_DrugName)
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  drug_moa <- gDRutils::get_env_identifiers("drug_moa")
  conc <- gDRutils::get_env_identifiers("concentration")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  gnumber_2 <- gDRutils::get_env_identifiers("drug2")
  conc_2 <- gDRutils::get_env_identifiers("concentration2")
  drug_moa_2 <- gDRutils::get_env_identifiers("drug_moa2")
  
  # select data for normalization type
  tab_response <- tab_response[normalization_type == metric_growth, ]
  
  if (fit_source %in% names(tab_response)) {
    data.table::setkeyv(tab_response, "fit_source")
    tab_response <- tab_response[fit_source]
    data.table::setkey(tab_response, NULL)
  }
  
  # fill column
  if (conc_2 %in% names(tab_response)) {
    tab_response <- tab_response[(!is.na(get(conc)) & !is.na(get(conc_2))), ]
  } else {
    tab_response <- tab_response[!is.na(get(conc)), ]
  }
  required_column <- c(gnumber_2, conc_2, drug_name_2, drug_moa_2)
  fill_val <- c(gnumber_2 = "untreated",
                conc_2 = 0.0,
                drug_name_2 = "untreated",
                drug_moa_2 = "untreated")
  names(fill_val) <- required_column
  need_to_be_filled <- required_column[!required_column %in% names(tab_response)]
  
  if (NROW(need_to_be_filled) > 0) {
    for (nm in need_to_be_filled) {
      tab_response[[nm]] <- fill_val[nm]
      if (nm == conc_2) tab_response[[conc_2]] <- as.numeric(tab_response[[conc_2]])
    }
  }
  
  # prep pivot data
  tab_response <- data.table::setorderv(tab_response, c(gnumber, gnumber_2, conc, conc_2))
  tab_response$col_pivot_name <- sprintf("%s_%s_%s_%s",
                                         tab_response[[gnumber]],
                                         tab_response[[gnumber_2]],
                                         tab_response[[conc]],
                                         tab_response[[conc_2]])
  col_pivot_name <- "col_pivot_name"
  tab_plot <- data.table::dcast(
    data = tab_response,
    formula = clid ~ col_pivot_name,
    value.var = metric
  )
  data.table::setcolorder(tab_plot, c(clid, unique(tab_response$col_pivot_name)))
  
  # prep matrix
  mat_cvd <- as.matrix(tab_plot[, .SD, .SDcols = -clid])
  rownames(mat_cvd) <- tab_plot[[clid]]
  rm_col <- vapply(colnames(mat_cvd), function(i) !all(is.na(mat_cvd[, i])), logical(1))
  rm_row <- vapply(seq_along(rownames(mat_cvd)), function(i) !all(is.na(mat_cvd[i, ])), logical(1))
  if (!all(rm_col)) mat_cvd <- mat_cvd[, rm_col]
  if (!all(rm_row)) mat_cvd <- mat_cvd[rm_row, ]
  if (metric == "x_std") mat_cvd <- mat_cvd ^ 2
  
  # annotation
  info_drug <- unique(tab_response[, .SD, .SDcols = c("col_pivot_name", gnumber, conc)])
  info_drug_2 <- unique(tab_response[, .SD, .SDcols = c("col_pivot_name", gnumber_2, conc_2)])
  data.table::setnames(info_drug_2, old = c(gnumber_2, conc_2), new = c(gnumber, conc))
  info_drug <- rbind(info_drug, info_drug_2)[get(gnumber) != "untreated"]
  drug_annotation <- data.table::dcast(
    data = info_drug,
    formula = col_pivot_name ~ get(gnumber), 
    value.var = conc
  )
  rownames(drug_annotation) <- drug_annotation$col_pivot_name # required by pheatmap::pheatmap
  drug_annotation <- drug_annotation[, .SD, .SDcol = -col_pivot_name]
  drug_annotation <- log10(drug_annotation)
  drug_annotation[drug_annotation == -Inf] <- NA # Q: when conc = 0
  
  # annotation colouring  
  drug_to_colored <- names(drug_annotation)
  ls_col <- get_qual_colors(NROW(drug_to_colored))
  drug_annotation_colors <- 
    lapply(seq_along(drug_to_colored), 
           function(i) { 
             c(colorspace::lighten(ls_col[i], 0.8, space = "HLS"), 
               colorspace::darken(ls_col[i], 0.1, space = "HLS"))
           }) 
  names(drug_annotation_colors) <- drug_to_colored
  
  # dendogram
  if (cluster_rows) {
    cluster_rows <- stats::hclust(stats::dist(mat_cvd))
  }
  
  # heatmap labels
  if (lbl_by_CellLineName) {
    row_lbls <- tab_response[, unique(.SD), .SDcols = c(cellline_name, clid)][order(rownames(mat_cvd))]
    # re-label
    rownames(mat_cvd) <- row_lbls[[cellline_name]]
  }
  
  if (lbl_by_DrugName) {
    col_lbls <- tab_response[, unique(.SD), .SDcols = c(drug_name, gnumber)]
    col_lbls_2 <- tab_response[, unique(.SD), .SDcols = c(drug_name_2, gnumber_2)]
    data.table::setnames(col_lbls_2, old = c(drug_name_2, gnumber_2), new = c(drug_name, gnumber))
    col_lbls <- rbind(col_lbls, col_lbls_2)
    # re-label
    colnames(drug_annotation) <- 
      col_lbls[get(gnumber) %in% colnames(drug_annotation), ][order(colnames(drug_annotation))][[drug_name]]
    names(drug_annotation_colors) <- 
      col_lbls[get(gnumber) %in% names(drug_annotation_colors), ][order(names(drug_annotation_colors))][[drug_name]]
  }
  
  annotation_legend_flag <- NROW(drug_to_colored) <= 3 # TODO Find better solution
  
  # prep hm color palette
  maxval <- switch(metric, "x" = 1.1, "x_std" = 0.5)
  minval <- min(c(0, round(min(stats::na.omit(mat_cvd)), digits = 2)))
  
  breaks <- seq(from = minval, to = maxval, length.out = no_breaks)
  hm_color_palette <- grDevices::colorRampPalette(colors_vec)(no_breaks + 1)
  if (metric == "x_std") hm_color_palette <- rev(hm_color_palette)
  
  hm <- pheatmap::pheatmap(
    mat = mat_cvd, 
    scale = "none", 
    display_numbers = FALSE, 
    number_color = "black",
    fontsize_number = 16, 
    color = hm_color_palette, 
    breaks = breaks, 
    angle_col = 45, 
    fontsize = 10, 
    show_colnames = FALSE,
    main = hm_title,
    na_col = "red",
    annotation_legend = annotation_legend_flag,
    # dendogram
    treeheight_row = 70, 
    treeheight_col = 70, 
    cluster_cols = FALSE, 
    cluster_rows = cluster_rows,
    # manual annotation
    annotation_col = drug_annotation, 
    annotation_colors = drug_annotation_colors 
  )
  return(hm)
}


#' Plot pretty heatmap with annotationsfor single-agent data
#'
#' @param tab_response \code{data.table} containing drug response metrics
#'    output from \code{\link[gDRutils]{convert_se_assay_to_dt}} for assay "Metrics" 
#'    and single-agent \code{SummarizedExperiment}
#' @param metric_growth string with normalization types to be selected
#'    one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param metric string name of metric;
#'    one of: "xc50" ("GR50" or "IC50" - respectively depending on \code{metric_growth}), 
#'    "x_max" ("GR Max" or "E Max") or "x_mean" ("GR Mean" or "RV Mean")
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
#'                       hm_title = get_hm_title(metric_growth = "GR",
#'                                               metric = "hsa_score",
#'                                               dataset_name = "Combo Matrix - combo data"))
#' 
#' @keywords QC_plot
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
  stopifnot("Must be valid color name" = all(vapply(colors_vec, gDRplots::is_valid_color, logical(1))))
  checkmate::assert_int(no_breaks, lower = 2)
  checkmate::assert_data_table(annotation_col, null.ok = TRUE)
  checkmate::assert_list(annotation_colors, null.ok = TRUE)
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  
  # select data for normalization type
  tab_response <- tab_response[normalization_type == metric_growth, ]
  
  if (fit_source %in% names(tab_response)) {
    data.table::setkeyv(tab_response, "fit_source")
    tab_response <- tab_response[fit_source]
    data.table::setkey(tab_response, NULL)
  }
  
  qmfun <- switch(metric,
                  "xc50" = log10, 
                  "x_max" = identity, 
                  "x_mean" = identity)
  
  # prep data
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
  
  hm <- pheatmap::pheatmap(mat = t_mat_cvd,
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
#'                          hm_title = get_hm_title(metric_growth = "GR",
#'                                                  metric = "hsa_score",
#'                                                  dataset_name = "Combo Matrix - combo data"))
#'             
#' @keywords QC_plot
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
  stopifnot("Must be valid color name" = all(vapply(colors_vec, gDRplots::is_valid_color, logical(1))))
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
#' @param metric_growth string with normalization types to be selected
#'    one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param metric string name of metric
#'    one of: "xc50"("GR50" or "IC50" - respectively depending on \code{metric_growth}), 
#'    "x_max" ("GR Max" or "E Max") or x_mean" ("GR Mean" or "RV Mean")
#' @param dataset_name string name of dataset
#'
#' @examples
#' get_hm_title(dataset_name = "Dateset DX123",
#'              metric = "x_mean", 
#'              metric_growth = "GR")
#' 
#' get_hm_title(metric = "xc50", 
#'              metric_growth = "GR")
#'              
#' @keywords QC_plot
#' 
#' @return character title for heatmap
#' @export 
get_hm_title <- function(metric = "xc50", 
                         metric_growth = "GR",
                         dataset_name = NULL) {
  
  checkmate::assert_string(dataset_name, null.ok = TRUE)
  checkmate::assert_choice(metric_growth, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, 
                           choices = c("xc50", "x_max", "x_mean", "hsa_score", "bliss_score"))
  
  title_metric <- 
    gDRutils::prettify_flat_metrics(sprintf("%s_%s", metric, metric_growth), human_readable = TRUE)
  
  if (metric == "xc50") title_metric <- sprintf("log10(%s)", title_metric)
  
  hm_title <- if (!is.null(dataset_name))  {
    sprintf("%s (%s)", dataset_name, title_metric)
  } else {
    title_metric
  }
  return(hm_title)
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


#' Create list of qualitative colors
#'
#' @param n number of required colors
#'
#' @return vector with hex colors from qualitative palettes
#' 
#' @examples
#' \dontrun{
#' get_qual_colors()
#' get_qual_colors(0)
#' get_qual_colors(5)
#' get_qual_colors(35) 
#' }
#' 
#' @keywords internal
#' 
get_qual_colors <- function(n = NULL) {
  checkmate::assert_int(n, null.ok = TRUE)
  
  if (!is.null(n) && n == 0) return("#000000")
  
  # list of colors: qualitative and friendly for user with color vision deficiency
  qual_col_pals <- RColorBrewer::brewer.pal.info[
    RColorBrewer::brewer.pal.info$category == "qual" & 
      RColorBrewer::brewer.pal.info$colorblind == TRUE, ]
  all_colors <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  if (is.null(n)) return(all_colors)
  
  if (n > length(all_colors)) {
    ls_light <- colorspace::lighten(all_colors, 0.3)
    ls_dark <-  colorspace::darken(all_colors, 0.3) # darker
    all_colors <- append(all_colors, values = c(ls_light, ls_dark))
  }
  
  rep(all_colors, length.out = n)
}
