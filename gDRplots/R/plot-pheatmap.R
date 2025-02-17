#' Plot pretty heatmap for single-agent or combo data to control quality of the data
#'    
#' @param dt_average  \code{data.table} representing data from the \code{Averaged} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Averaged")}
#'    and \code{SummarizedExperiment} with chosen data type: single-agent or combo
#' @param normalization_type string with normalization types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param metric string name of the metric;
#'    one of: "x" (value of "GR" or "RV" itself - respectively depending on \code{normalization_type}),
#'    or "x_std" (standard deviation)
#' @param fit_source string source name for metrics
#' @param hm_title string plot title
#' @param colors_vec character vector of colors (valid name or hex) used in heatmap
#'    note that for \code{metric} "x" the first color will be assigned to the min value, and 
#'    the last one - to the max; for "x_std" - that will be reversed
#' @param no_breaks numeric number of breaks on scale
#' @param cluster_rows logical flag whether rows should be clustered;
#'   the dendrogram will not be shown for the matrix with any dimension greater than 200.
#' @param distfun function used to compute the distance (dissimilarity) between rows;
#'   used for the dendrogram when \code{cluster_rows} is set to TRUE; 
#'   the default is \code{\link[gDRplots]{compute_distances}} using Spearman method.
#' @param lbl_by_CellLineName logical flag whether heatmap should be described by CellLineNames instead of clid
#' @param lbl_by_DrugName logical flag whether heatmap should be described by DrugName instead of Gnumber
#' 
#' @seealso \code{\link[pheatmap]{pheatmap}}
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]][2:5, ]
#' dt_average <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                                assay_name = "Averaged")
#' 
#' hm_1 <- pheatmap_qc(dt_average = dt_average)
#' 
#' hm_2 <- pheatmap_qc(dt_average = dt_average,
#'                     normalization_type = "RV",
#'                     colors_vec = c("darkblue", "grey90"),
#'                     lbl_by_CellLineName = TRUE,
#'                     lbl_by_DrugName = TRUE)
#' 
#' ggpubr::as_ggplot(hm_1[["gtable"]])
#' ggpubr::as_ggplot(hm_2[["gtable"]])
#' 
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' dt_average <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                                assay_name = "Averaged")
#' hm_3 <- pheatmap_qc(dt_average = dt_average,
#'                     cluster_rows = FALSE)
#' hm_4 <- pheatmap_qc(dt_average = dt_average,
#'                     metric = "x_std",
#'                     cluster_rows = FALSE)
#' 
#' ggpubr::as_ggplot(hm_3[["gtable"]])
#' ggpubr::as_ggplot(hm_4[["gtable"]])
#' 
#' @keywords QC_plot
#' 
#' @return \code{pheatmap} object containing heatmap for selected metric with annotation - if given
#' @export
pheatmap_qc <- function(
    dt_average,
    normalization_type = "GR",
    metric = "x",
    fit_source = "gDR",
    hm_title = NA,
    colors_vec = c("black", "grey99"),
    no_breaks = 50,
    cluster_rows = TRUE,
    distfun = compute_distances,
    lbl_by_CellLineName = FALSE,
    lbl_by_DrugName = FALSE) {
  
  checkmate::assert_data_table(dt_average)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, choices = c("x", "x_std"))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_string(hm_title, na.ok = TRUE)
  checkmate::assert_character(colors_vec)
  stopifnot("Must be a valid color name" = all(vapply(colors_vec, is_valid_color, logical(1))))
  checkmate::assert_int(no_breaks, lower = 2)
  checkmate::assert_flag(cluster_rows)
  checkmate::assert_function(distfun)
  checkmate::assert_flag(lbl_by_CellLineName)
  checkmate::assert_flag(lbl_by_DrugName)
  zero_conc_scaling_factor <- 
    gDRutils::get_settings_from_json("ZERO_CONC_SCALING_FACTOR",
                                     system.file(package = "gDRplots", "settings.json"))
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  conc <- gDRutils::get_env_identifiers("concentration")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  gnumber_2 <- gDRutils::get_env_identifiers("drug2")
  conc_2 <- gDRutils::get_env_identifiers("concentration2")
  drug_moa_2 <- gDRutils::get_env_identifiers("drug_moa2")
  
  # select data for normalization type
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  tab_response <- dt_average[eval(filter_expr)]
  
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
  if (!all(rm_col)) mat_cvd <- mat_cvd[, rm_col, drop = FALSE]
  if (!all(rm_row)) mat_cvd <- mat_cvd[rm_row, , drop = FALSE]
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
  # handle conc = 0
  min_val <-
    min(unlist(drug_annotation)[!is.na(unlist(drug_annotation)) & unlist(drug_annotation) != 0])
  drug_annotation[drug_annotation == 0] <- min_val / zero_conc_scaling_factor
  # log 10 (conc)
  drug_annotation <- log10(drug_annotation)
  
  # annotation coloring 
  drug_to_colored <- names(drug_annotation)
  ls_col <- get_qual_colors(NROW(drug_to_colored))
  drug_annotation_colors <-
    lapply(seq_along(drug_to_colored),
           function(i) {
             c(colorspace::lighten(ls_col[i], 0.8, space = "HLS"),
               colorspace::darken(ls_col[i], 0.1, space = "HLS"))
           })
  names(drug_annotation_colors) <- drug_to_colored
  
  # dendrogram
  max_dim_matrix_cluster <- 
    gDRutils::get_settings_from_json("MAX_DIM_MATRIX_CLUSTER",
                                     system.file(package = "gDRplots", "settings.json"))
  gDR_cluster_condition <- any(dim(mat_cvd) < max_dim_matrix_cluster)  # gDR standard
  if (cluster_rows) {
    cluster_rows <- .get_pheatmap_cluster_param(mat_to_cluster = mat_cvd,
                                                distfun = distfun,
                                                additional_condition = gDR_cluster_condition)
  }
  
  # heatmap labels
  if (lbl_by_CellLineName) {
    row_lbls <- tab_response[match(rownames(mat_cvd), clid), unique(.SD), .SDcols = c(cellline_name, clid)]
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
      col_lbls[get(gnumber) %in% colnames(drug_annotation), ][
        match(colnames(drug_annotation), get(gnumber)), ][[drug_name]]
    names(drug_annotation_colors) <-
      col_lbls[get(gnumber) %in% names(drug_annotation_colors), ][
        match(names(drug_annotation_colors), get(gnumber)), ][[drug_name]]
  }
  
  annotation_legend_flag <- NROW(drug_to_colored) <= 3 # TODO Find better solution
  
  # prep hm color palette
  maxval <- switch(metric, "x" = 1.1, "x_std" = 0.5)
  minval <- min(c(0, round(min(mat_cvd, na.rm = TRUE), digits = 2)))
  
  breaks <- seq(from = minval, to = maxval, length.out = no_breaks)
  hm_color_palette <- grDevices::colorRampPalette(colors_vec)(no_breaks + 1)
  if (metric == "x_std") hm_color_palette <- rev(hm_color_palette)
  
  hm <- 
    pheatmap::pheatmap(mat = mat_cvd,
                       scale = "none",
                       display_numbers = FALSE,
                       number_color = "black",
                       fontsize_number = 1.2 * 8,
                       color = hm_color_palette,
                       breaks = breaks,
                       angle_col = 45,
                       show_colnames = FALSE,
                       main = hm_title,
                       fontsize = 8,
                       fontsize_row = ifelse(NROW(mat_cvd) > 40, 0.6 * 8, 8),
                       na_col = "red",
                       annotation_legend = annotation_legend_flag,
                       # dendrogram
                       treeheight_row = 70,
                       treeheight_col = 70,
                       cluster_cols = FALSE,
                       cluster_rows = cluster_rows,
                       # manual annotation
                       annotation_col = drug_annotation,
                       annotation_colors = drug_annotation_colors,
                       silent = TRUE
    )
  return(hm)
}


#' Plot pretty heatmap with annotations for single-agent data
#'
#' @param dt_metrics \code{data.table} representing data from the \code{Metrics} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'  and single-agent \code{SummarizedExperiment}
#' @param dt_metrics_capped \code{data.table} representing data from the \code{Metrics} assay,
#'  the same as \code{dt_metrics} but with capped values with \code{\link[gDRutils]{cap_assay_infinities}}
#' @param normalization_type string with normalization types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param metric string name of the metric;
#'  one of: "xc50" ("GR50" or "IC50" - respectively depending on \code{normalization_type}), 
#'  "x_max" ("GR Max" or "E Max") or "x_mean" ("GR Mean" or "RV Mean"); 
#'  but the values from any numeric column can be displayed.
#' @param fit_source string source name for metrics
#' @param hm_title string plot title
#' @param colors_vec character vector of colors (valid name or hex) used in heatmap;
#'   note that the first color will be assigned to the min value, and the last one - to the max
#' @param no_breaks numeric number of breaks on scale
#' @param annotation_row \code{data.table} that specifies the annotations shown on left side of the heatmap.
#'   Each row defines the features for a specific row. The rows in the data and in the annotation
#'   are matched using corresponding names from the required  \code{DrugName} column.
#'   Note that color schemes takes into account if variable is continuous or discrete.
#' @param cluster_rows logical flag whether rows should be clustered;
#'   the dendrogram will not be shown for the matrix with any dimension greater than 200.
#' @param cluster_cols logical flag whether columns should be clustered;
#'   the dendrogram will not be shown for the matrix with any dimension greater than 200.
#' @param distfun function used to compute the distance (dissimilarity) between both rows and columns;
#'   used for the dendrogram when \code{cluster_rows} or \code{cluster_cols} is set to TRUE
#'   the default is \code{\link[gDRplots]{compute_distances}} using Spearman method.
#' @param annotation_col \code{data.table} that specifies the annotations shown above the heatmap.
#'   Each row defines the features for a specific column. The columns in the data and in the annotation
#'   are matched using corresponding names from the required  \code{CellLineName} column.
#'   Note that color schemes takes into account if variable is continuous or discrete.
#' @param annotation_colors named list for specifying \code{annotation_col} and \code{annotation_row} 
#'   track colors manually; note list is named with annotation name (column names of \code{annotation_row} - 
#'   without \code{DrugName} and column names of \code{annotation_col} - 
#'   without \code{CellLineName}), each list item is named vector with valid color name for 
#'   each value described in \code{annotation_row} and in \code{annotation_col} - respectively.
#'   Not described elements will be colored in default.
#' 
#' @seealso \code{\link[pheatmap]{pheatmap}}
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                                assay_name = "Metrics")
#' 
#' output <- pheatmap_with_anno_sa(dt_metrics = dt_metrics)
#' hm_1 <- output[["heatmap"]]
#' ggpubr::as_ggplot(hm_1[["gtable"]])
#' 
#' annotation_manual_col <-
#'   unique(dt_metrics[, c("CellLineName", "Tissue"), with = FALSE])
#' annotation_manual_row <-
#'   unique(dt_metrics[, c("DrugName", "drug_moa"), with = FALSE])
#' annotation_map <-
#'   get_ann_color_map(unique(dt_metrics[, c("Tissue", "drug_moa"), with = FALSE]))
#' 
#' output <- pheatmap_with_anno_sa(dt_metrics = dt_metrics,
#'                                 normalization_type = "RV",
#'                                 metric = "x_mean",
#'                                 colors_vec = c("darkblue", "grey90"),
#'                                 annotation_row = annotation_manual_row,
#'                                 annotation_col = annotation_manual_col,
#'                                 annotation_colors = annotation_map)
#' hm_2 <- output[["heatmap"]]
#' ggpubr::as_ggplot(hm_2[["gtable"]])
#' 
#' annotation_manual <- data.table::data.table(
#'   CellLineName =
#'     c("cellline_AA", "cellline_EA", "cellline_IB", "cellline_MC", "cellline_BC"),
#'   mut_A = c(1, 1, 1, 0, 0),
#'   mut_B = c("yes", "yes", "no", "no", "no")
#' )
#' annotation_map <- list(
#'   mut_A = c("1" = "coral", "0" = "cadetblue"),
#'   mut_B = c("yes" = "black", "no" = "grey90")
#' )
#' 
#' output <- pheatmap_with_anno_sa(dt_metrics = dt_metrics,
#'                                 annotation_col = annotation_manual,
#'                                 annotation_colors = annotation_map,
#'                                 hm_title = get_hm_title(
#'                                   normalization_type = "GR",
#'                                   metric = "hsa_score",
#'                                   dataset_name = "Combo Matrix - combo data"))
#' hm_3 <- output[["heatmap"]]
#' ggpubr::as_ggplot(hm_3[["gtable"]])
#' 
#' @keywords pheat_ann
#' 
#' @return A named list with elements:
#' \itemize{
#'   \item \code{data} a list containing the information visualized in the heatmap:
#'     \itemize{
#'       \item \code{matrix} data shown in the heatmap for the selected metric.
#'       \item \code{annotation_row} a table with row annotations (for \code{DrugName}), if provided.
#'       \item \code{annotation_col} a table with column annotations (for \code{CellLineName}), if provided.
#'     }
#'   \item \code{pheatmap} object containing the heatmap itself.
#' }
#' 
#' @export
pheatmap_with_anno_sa <- function(
    dt_metrics,
    dt_metrics_capped = NULL,
    normalization_type = "GR",
    metric = "xc50",
    fit_source = "gDR",
    hm_title = NA,
    colors_vec = NULL,
    no_breaks = 50,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    distfun = compute_distances,
    annotation_row = NULL,
    annotation_col = NULL,
    annotation_colors = NULL) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  
  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_data_table(dt_metrics_capped, null.ok = TRUE)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  numeric_columns <- names(dt_metrics)[vapply(dt_metrics, is.numeric, logical(1))]
  checkmate::assert_choice(metric, choices = numeric_columns)
  if (!is.null(dt_metrics_capped)) {
    numeric_columns <- names(dt_metrics_capped)[vapply(dt_metrics_capped, is.numeric, logical(1))]
    checkmate::assert_choice(metric, choices = numeric_columns)
  }
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_string(hm_title, na.ok = TRUE)
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  checkmate::assert_int(no_breaks, lower = 2)
  checkmate::assert_flag(cluster_rows)
  checkmate::assert_flag(cluster_cols)
  checkmate::assert_function(distfun)
  checkmate::assert_data_table(annotation_row, null.ok = TRUE)
  if (!is.null(annotation_row)) {
    checkmate::assert_names(names(annotation_row), must.include = drug_name)
  }
  checkmate::assert_data_table(annotation_col, null.ok = TRUE)
  if (!is.null(annotation_col)) {
    checkmate::assert_names(names(annotation_col), must.include = cellline_name)
  }
  checkmate::assert_list(annotation_colors, null.ok = TRUE)
  
  # output
  ls_output <- list(data = list(matrix = NULL,
                                annotation_col = NULL,
                                annotation_row = NULL),
                    heatmap = NULL)
  
  qmfun <- switch(metric,
                  "xc50" = log10,
                  identity)
  
  # prep matrix
  mat_cvd_raw <- prep_pheatmap_matrix(dt_response = dt_metrics,
                                      normalization_type = normalization_type,
                                      metric = metric,
                                      fit_source = fit_source,
                                      experiment_type = gDRutils::get_supported_experiments("sa"))
  
  mat_cvd <- if (is.null(dt_metrics_capped)) {
    mat_cvd_raw
  } else {
    prep_pheatmap_matrix(dt_response = dt_metrics_capped,
                         normalization_type = normalization_type,
                         metric = metric,
                         fit_source = fit_source,
                         experiment_type = gDRutils::get_supported_experiments("sa"))
  }
  
  # check completeness of annotation - TODO wrap in separate function
  if (!is.null(annotation_col)) {
    annotation_col <- .fill_pheatmap_annotation(dt_anno = annotation_col,
                                                mat_with_metric = mat_cvd,
                                                anno_var = cellline_name)
    ls_output[["data"]][["annotation_col"]] <- annotation_col
    
    rownames(annotation_col) <- annotation_col[[cellline_name]] # required by pheatmap::pheatmap
    annotation_col <- annotation_col[, .SD, .SDcol = -cellline_name]
    # order matrix
    mat_cvd <- mat_cvd[rownames(annotation_col), , drop = FALSE]
    if (!is.null(dt_metrics_capped)) {
      mat_cvd_raw <- mat_cvd_raw[rownames(annotation_col), , drop = FALSE]
    }
  }
  
  if (!is.null(annotation_row)) {
    annotation_row <- .fill_pheatmap_annotation(dt_anno = annotation_row,
                                                mat_with_metric = mat_cvd,
                                                anno_var = drug_name)
    ls_output[["data"]][["annotation_row"]] <- annotation_row
    
    rownames(annotation_row) <- annotation_row[[drug_name]] # required by pheatmap::pheatmap
    annotation_row <- annotation_row[, .SD, .SDcol = -drug_name]
    # order matrix
    mat_cvd <- mat_cvd[, rownames(annotation_row), drop = FALSE]
    if (!is.null(dt_metrics_capped)) {
      mat_cvd_raw <- mat_cvd_raw[, rownames(annotation_row), drop = FALSE]
    }
  }
  
  # filling missing values
  if (!is.null(annotation_col) && !is.null(annotation_colors)) {
    annotation_colors <- fill_ann_color_map(annotation_col, annotation_colors)
  }
  if (!is.null(annotation_row) && !is.null(annotation_colors)) {
    annotation_colors <- fill_ann_color_map(annotation_row, annotation_colors)
  }
  
  ls_output[["data"]][["matrix"]] <- data.table::as.data.table(mat_cvd, keep.rownames = cellline_name)
  # flip
  t_mat_cvd <- t(mat_cvd)
  t_mat_cvd[] <- vapply(t_mat_cvd, function(x) qmfun(x), numeric(1))
  
  # dendrogram
  max_dim_matrix_cluster <- 
    gDRutils::get_settings_from_json("MAX_DIM_MATRIX_CLUSTER",
                                     system.file(package = "gDRplots", "settings.json"))
  gDR_cluster_condition <- any(dim(t_mat_cvd) < max_dim_matrix_cluster)  # gDR standard
  if (cluster_rows) {
    cluster_rows <- .get_pheatmap_cluster_param(mat_to_cluster = t_mat_cvd,
                                                distfun = distfun,
                                                additional_condition = gDR_cluster_condition)
  }
  if (cluster_cols) {
    cluster_cols <- .get_pheatmap_cluster_param(mat_to_cluster = t(t_mat_cvd),
                                                distfun = distfun,
                                                additional_condition = gDR_cluster_condition)
  }
  
  # prep hm color palette
  min_val <- ifelse(
    all(is.infinite(t_mat_cvd) | is.na(t_mat_cvd)),
    -1,
    min(t_mat_cvd[!is.infinite(t_mat_cvd)], na.rm = TRUE)
  )
  max_val <- ifelse(metric %in% c("x", "xc50", "x_max", "x_mean"), 1.0, max(t_mat_cvd, na.rm = TRUE))
  
  if (min_val == max_val) {
    min_val <- min_val - 1
  }
  
  breaks <- seq(from = min_val, to = max_val, length.out = no_breaks)
  hm_color_palette <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
    .get_smooth_palette(no_breaks)
  } else {
    grDevices::colorRampPalette(colors_vec)(no_breaks + 1)
  }
  
  # display numbers - for readability, turn it off for matrices larger than 15x15
  display_numbers_flag <- !any(dim(t_mat_cvd) > c(15, 15))
  
  ls_output[["heatmap"]] <- 
    pheatmap::pheatmap(mat = t_mat_cvd,
                       scale = "none",
                       display_numbers = display_numbers_flag,
                       number_color = "black",
                       fontsize_number = 1.2 * 8,
                       color = hm_color_palette,
                       breaks = breaks,
                       angle_col = 90,
                       main = hm_title,
                       fontsize = 8,
                       fontsize_col = ifelse(NCOL(t_mat_cvd) > 40, 0.6 * 8, 8),
                       na_col = "darkgray",
                       # dendrogram
                       cluster_rows = cluster_rows,
                       cluster_cols = cluster_cols,
                       treeheight_row = 25,
                       treeheight_col = 25,
                       # manual annotation
                       annotation_row = annotation_row,
                       annotation_col = annotation_col,
                       annotation_colors = annotation_colors,
                       silent = TRUE)
  
  return(ls_output)
}

#' Plot pretty heatmap with annotations for co-dilution data
#'
#' @param dt_metrics \code{data.table} representing data from the \code{Metrics} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'  and co-dilution \code{SummarizedExperiment}
#' @param normalization_type string with normalization types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param metric string name of the metric;
#'  one of: "xc50" ("GR50" or "IC50" - respectively depending on \code{normalization_type}), 
#'  "x_max" ("GR Max" or "E Max") or "x_mean" ("GR Mean" or "RV Mean")
#' @param fit_source string source name for metrics
#' @param hm_title string plot title
#' @param colors_vec character vector of colors (valid name or hex) used in heatmap;
#'   note that the first color will be assigned to the min value, and the last one - to the max
#' @param no_breaks numeric number of breaks on scale
#' @param annotation_row \code{data.table} that specifies the annotations shown on the left side 
#'   of the heatmap.
#'   Each row defines the features for a specific row. The rows in the data and in the annotation
#'   are matched using corresponding names from the required  \code{DrugName}, 
#'   \code{DrugName_2} and \code{Concentration_2} columns.
#'   Note that color schemes takes into account if the variable is continuous or discrete.
#' @param cluster_rows logical flag whether indicating rows should be clustered;
#'   the dendrogram will not be shown for the matrix with any dimension greater than 200.
#' @param cluster_cols logical flag indicating whether columns should be clustered;
#'   the dendrogram will not be shown for the matrix with any dimension greater than 200.
#' @param distfun function used to compute the distance (dissimilarity) between both rows and columns;
#'   used for the dendrogram when \code{cluster_rows} or \code{cluster_cols} is set to TRUE
#'   the default is \code{\link[gDRplots]{compute_distances}} using Spearman method.
#' @param annotation_col \code{data.table} that specifies the annotations shown above the heatmap.
#'   Each row defines the features for a specific column. The columns in the data and in the annotation
#'   are matched using corresponding names from the required  \code{CellLineName} column.
#'   Note that color schemes takes into account if the variable is continuous or discrete.
#' @param annotation_colors named list for specifying \code{annotation_col} and \code{annotation_row} 
#'   track colors manually; note list is named with annotation name (column names of \code{annotation_row} - 
#'   without \code{DrugName} and column names of \code{annotation_col} - 
#'   without \code{CellLineName}), each list item is named vector with valid color name for 
#'   each value described in \code{annotation_row} and in \code{annotation_col}, respectively.
#'   Not described elements will be colored by default.
#' 
#' @seealso \code{\link[pheatmap]{pheatmap}}
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_codilution_small")
#' se <- mae[[gDRutils::get_supported_experiments("cd")]]
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                                assay_name = "Metrics")
#' 
#' output <- pheatmap_with_anno_cd(dt_metrics = dt_metrics)
#' hm_1 <- output[["heatmap"]]
#' ggpubr::as_ggplot(hm_1[["gtable"]])
#' 
#' annotation_manual_col <-
#'   unique(dt_metrics[, c("CellLineName", "Tissue"), with = FALSE])
#' annotation_manual_row <-
#'   unique(dt_metrics[, c("DrugName", "DrugName_2", "Concentration_2",
#'                         "drug_moa", "drug_moa_2"),
#'                     with = FALSE])
#' annotation_map <-
#'   get_ann_color_map(unique(dt_metrics[, c("Tissue", "drug_moa", "drug_moa_2"), with = FALSE]))
#' 
#' output <- pheatmap_with_anno_cd(dt_metrics = dt_metrics,
#'                                 normalization_type = "RV",
#'                                 metric = "x_mean",
#'                                 colors_vec = c("darkblue", "grey90"),
#'                                 annotation_row = annotation_manual_row,
#'                                 annotation_col = annotation_manual_col,
#'                                 annotation_colors = annotation_map)
#' hm_2 <- output[["heatmap"]]
#' ggpubr::as_ggplot(hm_2[["gtable"]])
#' 
#' annotation_manual <- data.table::data.table(
#'   CellLineName =
#'     c("cellline_AA", "cellline_EA", "cellline_IB", "cellline_MC", "cellline_BC"),
#'   mut_A = c(1, 1, 1, 0, 0),
#'   mut_B = c("yes", "yes", "no", "no", "no")
#' )
#' annotation_map <- list(
#'   mut_A = c("1" = "coral", "0" = "cadetblue"),
#'   mut_B = c("yes" = "black", "no" = "grey90")
#' )
#' 
#' output <- pheatmap_with_anno_cd(dt_metrics = dt_metrics,
#'                                 annotation_col = annotation_manual,
#'                                 annotation_colors = annotation_map,
#'                                 hm_title = get_hm_title(
#'                                   normalization_type = "GR",
#'                                   metric = "hsa_score",
#'                                   dataset_name = "Co-dilution data"))
#' hm_3 <- output[["heatmap"]]
#' ggpubr::as_ggplot(hm_3[["gtable"]])
#' 
#' @keywords pheat_ann
#' 
#' @return A named list with elements:
#' \itemize{
#'   \item \code{data} a list containing the information visualized in the heatmap:
#'     \itemize{
#'       \item \code{matrix} data shown in the heatmap for the selected metric.
#'       \item \code{annotation_row} a table with row annotations (for \code{DrugName}), if provided.
#'       \item \code{annotation_col} a table with column annotations (for \code{CellLineName}), if provided.
#'     }
#'   \item \code{pheatmap} object containing the heatmap itself.
#' }
#' 
#' @export
pheatmap_with_anno_cd <- function(
    dt_metrics,
    normalization_type = "GR",
    metric = "xc50",
    fit_source = "gDR",
    hm_title = NA,
    colors_vec = NULL,
    no_breaks = 50,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    distfun = compute_distances,
    annotation_row = NULL,
    annotation_col = NULL,
    annotation_colors = NULL) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  conc_2 <- gDRutils::get_env_identifiers("concentration2")
  
  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  numeric_columns <- names(dt_metrics)[vapply(dt_metrics, is.numeric, logical(1))]
  checkmate::assert_choice(metric, choices = numeric_columns)
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_string(hm_title, na.ok = TRUE)
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  checkmate::assert_int(no_breaks, lower = 2)
  checkmate::assert_flag(cluster_rows)
  checkmate::assert_flag(cluster_cols)
  checkmate::assert_function(distfun)
  checkmate::assert_data_table(annotation_row, null.ok = TRUE)
  if (!is.null(annotation_row)) {
    checkmate::assert_names(names(annotation_row), must.include = c(drug_name, drug_name_2, conc_2))
  }
  checkmate::assert_data_table(annotation_col, null.ok = TRUE)
  if (!is.null(annotation_col)) {
    checkmate::assert_names(names(annotation_col), must.include = cellline_name)
  }
  checkmate::assert_list(annotation_colors, null.ok = TRUE)
  
  # output
  ls_output <- list(data = list(matrix = NULL,
                                annotation_col = NULL,
                                annotation_row = NULL),
                    heatmap = NULL)
  
  qmfun <- switch(metric,
                  "xc50" = log10,
                  identity)
  
  # prep matrix
  mat_cvd <- prep_pheatmap_matrix(dt_response = dt_metrics,
                                  normalization_type = normalization_type,
                                  metric = metric,
                                  fit_source = fit_source,
                                  experiment_type = gDRutils::get_supported_experiments("cd"))
  
  # check completeness of annotation - TODO wrap in separate function
  if (!is.null(annotation_col)) {
    if (!all(rownames(mat_cvd) %in% annotation_col[[cellline_name]])) {
      tab_missing_ann <- data.table::data.table(
        missing = rownames(mat_cvd)[!rownames(mat_cvd) %in% annotation_col[[cellline_name]]]
      )
      data.table::setnames(tab_missing_ann, "missing", cellline_name)
      
      annotation_col <- data.table::rbindlist(list(annotation_col, tab_missing_ann), fill = TRUE)
    }
    # data.table::nafill does not support character
    cols <- names(annotation_col)[names(annotation_col) != cellline_name]
    data.table::setorderv(annotation_col, cols = cols, na.last = TRUE)
    annotation_col[, (cols) := lapply(.SD, change_NA_into_char, "NA"), .SDcols = cols]
    # select annotation acc to matrix
    annotation_col <- annotation_col[get(cellline_name) %in% rownames(mat_cvd), ]
    ls_output[["data"]][["annotation_col"]] <- annotation_col
    
    rownames(annotation_col) <- annotation_col[[cellline_name]] # required by pheatmap::pheatmap
    annotation_col <- annotation_col[, .SD, .SDcol = -cellline_name]
    # order matrix
    mat_cvd <- mat_cvd[rownames(annotation_col), , drop = FALSE]
  }
  
  if (!is.null(annotation_row)) {
    annotation_row$DrugCombination <-
      paste(annotation_row[[drug_name]], "x", paste0(annotation_row[[drug_name_2]], "__", annotation_row[[conc_2]]))
    
    if (!all(colnames(mat_cvd) %in% annotation_row[["DrugCombination"]])) {
      tab_missing_ann <- data.table::data.table(
        missing = colnames(mat_cvd)[!colnames(mat_cvd) %in% annotation_row[["DrugCombination"]]]
      )
      data.table::setnames(tab_missing_ann, "missing", "DrugCombination")
      tab_missing_ann[, c(drug_name, drug_name_2) := data.table::tstrsplit(DrugCombination, " x ", fixed = TRUE)]
      tab_missing_ann[, c(drug_name_2, conc_2) := data.table::tstrsplit(get(drug_name_2), "__", fixed = TRUE)]
      tab_missing_ann[[conc_2]] <- as.numeric(tab_missing_ann[[conc_2]])
      
      annotation_row <- data.table::rbindlist(list(annotation_row, tab_missing_ann), fill = TRUE)
    }
    # data.table::nafill does not support character
    cols <- names(annotation_row)[!names(annotation_row) %in% c(drug_name, drug_name_2, conc_2, "DrugCombination")]
    data.table::setorderv(annotation_row, cols = cols, na.last = TRUE)
    annotation_row[, (cols) := lapply(.SD, change_NA_into_char, "NA"), .SDcols = cols, drop = FALSE]
    # select annotation acc to matrix
    annotation_row <- annotation_row[DrugCombination %in% colnames(mat_cvd), ]
    ls_output[["data"]][["annotation_row"]] <- annotation_row[, !c("DrugCombination"), with = FALSE]
    
    rownames(annotation_row) <- annotation_row[["DrugCombination"]] # required by pheatmap::pheatmap
    annotation_row <- 
      annotation_row[, .SD, .SDcols = -c(drug_name, drug_name_2, conc_2, "DrugCombination")]
    # order matrix
    mat_cvd <- mat_cvd[, rownames(annotation_row), drop = FALSE]
  }
  
  # filling missing values
  if (!is.null(annotation_col) && !is.null(annotation_colors)) {
    annotation_colors <- fill_ann_color_map(annotation_col, annotation_colors)
  }
  if (!is.null(annotation_row) && !is.null(annotation_colors)) {
    annotation_colors <- fill_ann_color_map(annotation_row, annotation_colors)
  }
  
  ls_output[["data"]][["matrix"]] <- data.table::as.data.table(mat_cvd, keep.rownames = cellline_name)
  # flip
  t_mat_cvd <- t(mat_cvd)
  t_mat_cvd[] <- vapply(t_mat_cvd, function(x) purrr::quietly(qmfun)(x)$result, numeric(1))
  
  # dendrogram
  max_dim_matrix_cluster <- 
    gDRutils::get_settings_from_json("MAX_DIM_MATRIX_CLUSTER",
                                     system.file(package = "gDRplots", "settings.json"))
  gDR_cluster_condition <- any(dim(t_mat_cvd) < max_dim_matrix_cluster)  # gDR standard
  if (cluster_rows) {
    cluster_rows <- .get_pheatmap_cluster_param(mat_to_cluster = t_mat_cvd,
                                                distfun = distfun,
                                                additional_condition = gDR_cluster_condition)
  }
  if (cluster_cols) {
    cluster_cols <- .get_pheatmap_cluster_param(mat_to_cluster = t(t_mat_cvd),
                                                distfun = distfun,
                                                additional_condition = gDR_cluster_condition)
  }
  
  # prep hm color palette
  min_val <- min(t_mat_cvd, na.rm = TRUE)
  max_val <- ifelse(metric %in% c("x", "xc50", "x_max", "x_mean"), 1.0, max(t_mat_cvd, na.rm = TRUE))
  
  if (min_val == max_val) {
    min_val <- min_val - 1
  } else if (is.infinite(min_val)) {
    min_val <- max_val / 100
  }
  
  breaks <- seq(from = min_val, to = max_val, length.out = no_breaks)
  hm_color_palette <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
    .get_smooth_palette(no_breaks)
  } else {
    grDevices::colorRampPalette(colors_vec)(no_breaks + 1)
  }
  
  # display numbers - for readability, turn it off for matrices larger than 15x15
  display_numbers_flag <- !any(dim(t_mat_cvd) > c(15, 15))
  
  ls_output[["heatmap"]] <- 
    pheatmap::pheatmap(mat = t_mat_cvd,
                       scale = "none",
                       display_numbers = display_numbers_flag,
                       number_color = "black",
                       fontsize_number = 1.2 * 8,
                       color = hm_color_palette,
                       breaks = breaks,
                       angle_col = 90,
                       main = hm_title,
                       fontsize = 8,
                       fontsize_col = ifelse(NCOL(t_mat_cvd) > 40, 0.6 * 8, 8),
                       na_col = "darkgray",
                       # dendrogram
                       cluster_rows = cluster_rows,
                       cluster_cols = cluster_cols,
                       treeheight_row = 25,
                       treeheight_col = 25,
                       # manual annotation
                       annotation_row = annotation_row,
                       annotation_col = annotation_col,
                       annotation_colors = annotation_colors,
                       silent = TRUE)
  
  return(ls_output)
}

#' Plot pretty heatmap with annotations for combo data
#' 
#' @param dt_scores \code{data.table} representing data from the \code{scores} assay,
#'   outputted by \code{gDRutils::convert_se_assay_to_dt(se, "scores")}
#'   and combo \code{SummarizedExperiment}
#' @param metric string name of the combo metric;
#'   one of: "hsa_score"("Bliss Excess GR" or "Bliss Excess RV" - respectively 
#'   depending on \code{normalization_type}), "bliss_score" ("Bliss Score GR" or "Bliss Score RV")
#' @param annotation_row \code{data.table} that specifies the annotations shown on left side of the heatmap.
#'   Each row defines the features for a specific row. The rows in the data and in the annotation
#'   are matched using corresponding combination of names from \code{DrugName} and \code{DrugName_2} 
#'   columns. Both columns are required.
#'   Note that color schemes takes into account if variable is continuous or discrete.
#' @param annotation_colors named list for specifying \code{annotation_col} and \code{annotation_row} 
#'   track colors manually; note list is named with annotation name (column names of \code{annotation_row} - 
#'   without \code{DrugName} and \code{DrugName_2}, and column names of \code{annotation_col} - 
#'   without \code{CellLineName}), each list item is named vector with a valid color name for 
#'   each value described in \code{annotation_row}) and in \code{annotation_col}) - respectively.
#'   Not described elements will be colored in default.
#' @inheritParams pheatmap_with_anno_sa
#' 
#' @seealso \code{\link[pheatmap]{pheatmap}}
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' dt_scores <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                               assay_name = "scores")
#' annotation_manual_col <-
#'   unique(dt_scores[, c("CellLineName", "Tissue"), with = FALSE])
#' annotation_manual_row <-
#'   unique(dt_scores[, c("DrugName", "DrugName_2", "drug_moa", "drug_moa_2"), with = FALSE])
#' annotation_map <-
#'   get_ann_color_map(unique(dt_scores[, c("Tissue", "drug_moa", "drug_moa_2"), with = FALSE]))
#' 
#' output <- pheatmap_with_anno_combo(dt_scores = dt_scores,
#'                                    normalization_type = "RV",
#'                                    metric = "bliss_score",
#'                                    colors_vec = c("darkblue", "grey90", "darkred"),
#'                                    annotation_row = annotation_manual_row,
#'                                    annotation_col = annotation_manual_col,
#'                                    annotation_colors = annotation_map)
#' hm_1 <- output[["heatmap"]]
#' ggpubr::as_ggplot(hm_1[["gtable"]])
#' 
#' annotation_manual <- data.table::data.table(
#'   CellLineName =
#'     c("cellline_AA", "cellline_EA", "cellline_IB", "cellline_MC", "cellline_BC"),
#'   mut_A = c(1, 1, 1, 0, 0),
#'   mut_B = c("yes", "yes", "no", "no", "no")
#' )
#' 
#' annotation_map <- list(
#'   mut_A = c("1" = "coral", "0" = "cadetblue"),
#'   mut_B = c("yes" = "black", "no" = "grey90")
#' )
#' 
#' output <- pheatmap_with_anno_combo(dt_scores = dt_scores,
#'                                    cluster_cols = FALSE,
#'                                    annotation_col = annotation_manual,
#'                                    annotation_colors = annotation_map,
#'                                    hm_title = get_hm_title(
#'                                      normalization_type = "GR",
#'                                      metric = "hsa_score",
#'                                      dataset_name = "Combo Matrix - combo data"))
#' hm_2 <- output[["heatmap"]]
#' ggpubr::as_ggplot(hm_2[["gtable"]])
#'             
#' @keywords pheat_ann
#' 
#' @return A named list with elements:
#' \itemize{
#'   \item \code{data} a list containing the information visualized in the heatmap:
#'     \itemize{
#'       \item \code{matrix} data shown in the heatmap for the selected metric.
#'       \item \code{annotation_row} a table with row annotations (for \code{DrugName}), if provided.
#'       \item \code{annotation_col} a table with column annotations (for \code{CellLineName}), if provided.
#'     }
#'   \item \code{pheatmap} object containing the heatmap itself.
#' }
#' 
#' @export
pheatmap_with_anno_combo <- function(
    dt_scores,
    normalization_type = "GR",
    metric = "hsa_score",
    fit_source = "gDR",
    hm_title = NA,
    colors_vec = NULL,
    no_breaks = 50,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    distfun = compute_distances,
    annotation_row = NULL,
    annotation_col = NULL,
    annotation_colors = NULL) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  
  checkmate::assert_data_table(dt_scores)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, choices = c("hsa_score", "bliss_score"))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_string(hm_title, na.ok = TRUE)
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  checkmate::assert_int(no_breaks, lower = 2)
  checkmate::assert_flag(cluster_rows)
  checkmate::assert_flag(cluster_cols)
  checkmate::assert_function(distfun)
  checkmate::assert_data_table(annotation_row, null.ok = TRUE)
  if (!is.null(annotation_row)) {
    checkmate::assert_names(names(annotation_row), must.include = c(drug_name, drug_name_2))
  }
  checkmate::assert_data_table(annotation_col, null.ok = TRUE)
  if (!is.null(annotation_col)) {
    checkmate::assert_names(names(annotation_col), must.include = cellline_name)
  }
  checkmate::assert_list(annotation_colors, null.ok = TRUE)
  
  # output
  ls_output <- list(data = list(matrix = NULL,
                                annotation_col = NULL,
                                annotation_row = NULL),
                    heatmap = NULL)
  
  # prep matrix
  mat_cvd <- prep_pheatmap_matrix(dt_response = dt_scores,
                                  normalization_type = normalization_type,
                                  metric = metric,
                                  fit_source = fit_source,
                                  experiment_type = gDRutils::get_supported_experiments("combo"))
  
  # check completeness of annotation - TODO wrap in separate function
  if (!is.null(annotation_col)) {
    if (!all(rownames(mat_cvd) %in% annotation_col[[cellline_name]])) {
      tab_missing_ann <- data.table::data.table(
        missing = rownames(mat_cvd)[!rownames(mat_cvd) %in% annotation_col[[cellline_name]]]
      )
      data.table::setnames(tab_missing_ann, "missing", cellline_name)
      
      annotation_col <- data.table::rbindlist(list(annotation_col, tab_missing_ann), fill = TRUE)
    }
    # data.table::nafill does not support character
    cols <- names(annotation_col)[names(annotation_col) != cellline_name]
    data.table::setorderv(annotation_col, cols = cols, na.last = TRUE)
    annotation_col[, (cols) := lapply(.SD, change_NA_into_char, "NA"), .SDcols = cols]
    # select annotation acc to matrix
    annotation_col <- annotation_col[get(cellline_name) %in% rownames(mat_cvd), ]
    ls_output[["data"]][["annotation_col"]] <- annotation_col
    
    rownames(annotation_col) <- annotation_col[[cellline_name]] # required by pheatmap::pheatmap
    annotation_col <- annotation_col[, .SD, .SDcol = -cellline_name]
    # order matrix
    mat_cvd <- mat_cvd[rownames(annotation_col), , drop = FALSE]
  }
  
  if (!is.null(annotation_row)) {
    annotation_row$DrugCombination <-
      paste(annotation_row[[drug_name]], "x", annotation_row[[drug_name_2]])
    
    if (!all(colnames(mat_cvd) %in% annotation_row[["DrugCombination"]])) {
      tab_missing_ann <- data.table::data.table(
        missing = colnames(mat_cvd)[!colnames(mat_cvd) %in% annotation_row[["DrugCombination"]]]
      )
      data.table::setnames(tab_missing_ann, "missing", "DrugCombination")
      tab_missing_ann[, c(drug_name, drug_name_2) := data.table::tstrsplit(DrugCombination, " x ", fixed = TRUE)]
      
      annotation_row <- data.table::rbindlist(list(annotation_row, tab_missing_ann), fill = TRUE)
    }
    # data.table::nafill does not support character
    cols <- names(annotation_row)[!names(annotation_row) %in% c(drug_name, drug_name_2, "DrugCombination")]
    data.table::setorderv(annotation_row, cols = cols, na.last = TRUE)
    annotation_row[, (cols) := lapply(.SD, change_NA_into_char, "NA"), .SDcols = cols, drop = FALSE]
    # select annotation acc to matrix
    annotation_row <- annotation_row[DrugCombination %in% colnames(mat_cvd), ]
    ls_output[["data"]][["annotation_row"]] <- annotation_row[, !c("DrugCombination"), with = FALSE]
    
    rownames(annotation_row) <- annotation_row[["DrugCombination"]] # required by pheatmap::pheatmap
    annotation_row <- annotation_row[, .SD, .SDcol = -c(drug_name, drug_name_2, "DrugCombination")]
    # order matrix
    mat_cvd <- mat_cvd[, rownames(annotation_row), drop = FALSE]
  }
  
  # filling missing values
  if (!is.null(annotation_col) && !is.null(annotation_colors)) {
    annotation_colors <- fill_ann_color_map(annotation_col, annotation_colors)
  }
  if (!is.null(annotation_row) && !is.null(annotation_colors)) {
    annotation_colors <- fill_ann_color_map(annotation_row, annotation_colors)
  }
  
  ls_output[["data"]][["matrix"]] <- data.table::as.data.table(mat_cvd, keep.rownames = cellline_name)
  # flip
  t_mat_cvd <- t(mat_cvd)
  
  # dendrogram
  max_dim_matrix_cluster <- 
    gDRutils::get_settings_from_json("MAX_DIM_MATRIX_CLUSTER",
                                     system.file(package = "gDRplots", "settings.json"))
  gDR_cluster_condition <- any(dim(t_mat_cvd) < max_dim_matrix_cluster)  # gDR standard
  if (cluster_rows) {
    cluster_rows <- .get_pheatmap_cluster_param(mat_to_cluster = t_mat_cvd,
                                                distfun = distfun,
                                                additional_condition = gDR_cluster_condition)
  }
  if (cluster_cols) {
    cluster_cols <- .get_pheatmap_cluster_param(mat_to_cluster = t(t_mat_cvd),
                                                distfun = distfun,
                                                additional_condition = gDR_cluster_condition)
  }
  
  # prep hm color palette
  breaks <- seq(from = -0.7, to = 0.7, length.out = no_breaks)
  hm_color_palette <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
    .get_excess_palette(no_breaks)
  } else {
    grDevices::colorRampPalette(colors_vec)(no_breaks + 1)
  }
  
  # display numbers - for readability, turn it off for matrices larger than 15x15
  display_numbers_flag <- !any(dim(t_mat_cvd) > c(15, 15))
  
  ls_output[["heatmap"]] <- 
    pheatmap::pheatmap(t_mat_cvd,
                       scale = "none",
                       display_numbers = display_numbers_flag,
                       number_color = "black",
                       fontsize_number = 1.2 * 8,
                       color = hm_color_palette,
                       breaks = breaks,
                       angle_col = 90,
                       main = hm_title,
                       fontsize = 8,
                       fontsize_col = ifelse(NCOL(t_mat_cvd) > 40, 0.6 * 8, 8),
                       na_col = "darkgray",
                       # dendrogram
                       cluster_rows = cluster_rows,
                       cluster_cols = cluster_cols,
                       treeheight_row = 25,
                       treeheight_col = 25,
                       # manual annotation
                       annotation_row = annotation_row,
                       annotation_col = annotation_col,
                       annotation_colors = annotation_colors,
                       silent = TRUE)
  return(ls_output)
}

# helpers ----
#' Get Legend Title
#' 
#' @param normalization_type string with normalization types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param metric string name of the metric
#'    one of: "xc50"("GR50" or "IC50" - respectively depending on \code{normalization_type}),
#'    "x_max" ("GR Max" or "E Max") or x_mean" ("GR Mean" or "RV Mean")
#' @param dataset_name string name of dataset
#'
#' @examples
#' get_hm_title(dataset_name = "Dateset DX123",
#'              metric = "x_mean",
#'              normalization_type = "GR")
#' 
#' get_hm_title(metric = "xc50",
#'              normalization_type = "GR")
#'              
#' @keywords pheat_ann
#' 
#' @return character title for heatmap
#' @export
get_hm_title <- function(metric = "xc50",
                         normalization_type = "GR",
                         dataset_name = NULL) {
  
  checkmate::assert_string(dataset_name, null.ok = TRUE)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  # Allow for using any column (even custom), e.g. xc50_sd
  checkmate::assert_character(metric)
  title_metric <-
    gDRutils::prettify_flat_metrics(sprintf("%s_%s", metric, normalization_type), human_readable = TRUE)
  
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
change_NA_into_char <- function(x,
                                lbl_NA = "NA") {
  
  checkmate::assert_string(lbl_NA)
  
  ifelse(is.na(x), lbl_NA, as.character(x))
}


#' Create color map for annotation
#'
#' @param dt_ann \code{data.table} with the annotations
#'
#' @return list with color mapping for the annotations
#' 
#' @seealso \code{\link{pheatmap_qc}}
#' 
#' @examples
#' mae <- gDRutils::get_synthetic_data("small")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]][2:5, ]
#' dt_average <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "Averaged")
#' dt_ann <- dt_average[,.SD, .SDcols = c("Tissue", "ReferenceDivisionTime")]
#' 
#' get_ann_color_map(dt_ann)
#' 
#' @keywords utils_color
#' @export 
get_ann_color_map <- function(dt_ann) {
  checkmate::assert_data_table(dt_ann)
  
  n_unique_ann <- sum(
    vapply(names(dt_ann), function(nm) NROW(unique(dt_ann[[nm]])), FUN.VALUE = numeric(1)))
  ls_colors <- get_qual_colors(n_unique_ann)
  
  annotation_colors <- list()
  for (ann in names(dt_ann)) {
    lvl <- as.character(unique(dt_ann[[ann]]))
    col_map <- ls_colors[seq_along(lvl)]
    names(col_map) <- lvl
    
    ls_colors <- ls_colors[-seq_along(lvl)]
    annotation_colors[[ann]] <- col_map
  }
  annotation_colors
}


#' Fill missing values in the color map for annotation
#'
#' @param dt_ann \code{data.table} with the annotations
#' @param map_ann \code{list} with the annotations
#'
#' @return list with color mapping for the annotations with missing items filled in
#' 
#' @seealso \code{\link{pheatmap_with_anno_sa}} \code{\link{pheatmap_with_anno_combo}}
#' 
#' @examples
#' annotation_manual <- data.table::data.table(
#'   CellLineName = c("cellline_AA", "cellline_EA", "cellline_IB", "cellline_MC", "cellline_BC"),
#'   mut_A = c(0, 0, 1, 2, 3),
#'   mut_B = c("yes", "yes", "no", NA, NA),
#'   mut_C = c("AA", "AA", "AB", NA, "B")
#' )
#' 
#' annotation_map <- list(
#'   mut_A = c("1" = "coral", "0" = "cadetblue"),
#'   mut_B = c("yes" = "black", "no" = "grey90", "not_checked" = "lightblue"),
#'   mut_C = c("AA" = "red", "AB" = "orange", "B" = "yellow")
#' )
#' 
#' fill_ann_color_map(dt_ann = annotation_manual, map_ann = annotation_map)
#' 
#' @keywords utils_color
#' @export 
fill_ann_color_map <- function(dt_ann,
                               map_ann) {
  checkmate::assert_data_table(dt_ann)
  checkmate::assert_list(map_ann)
  
  dt_ann <- dt_ann[, lapply(.SD, change_NA_into_char)] # annotation has to be character type without NA
  
  ls_ann_with_colors <- names(dt_ann)[names(dt_ann) %in% names(map_ann)]
  
  if (NROW(ls_ann_with_colors) > 0) {
    for (ann in ls_ann_with_colors) {
      required_lvl <- unique(dt_ann[[ann]])
      available_lvl <- names(map_ann[[ann]])
      missing_lvl <- required_lvl[!required_lvl %in% available_lvl]
      
      if (any(required_lvl == "NA")) {
        required_lvl <- c(required_lvl[required_lvl != "NA"], "NA")
      }
      
      if (NROW(missing_lvl) == 1 && missing_lvl == "NA") {
        col_na <- ifelse(any(map_ann[[ann]] %in% c("black", "#000000")), "darkred", "black")
        map_ann[[ann]] <- c(map_ann[[ann]], "NA" = col_na)
      } else if (NROW(missing_lvl) > 0) {
        map_ann[[ann]] <- NULL # allow default coloring
      }
      map_ann[[ann]] <- map_ann[[ann]][required_lvl]
    }
  }
  return(map_ann)
}

#' Compute value of param cluster_rows or cluster_cols in pheatmap::pheatmap
#' 
#' The \code{cluster_rows} and \code{cluster_cols} parameters  pheatmap::pheatmap may take values:
#' - boolean determining if rows/columns should be clustered
#' - \code{hclust} object,
#' this function allows to calculate proper value taking into account matrix, additional condition and 
#' selected function used to compute the distance in \code{hclust} object
#' 
#' To compute the correct value when clustering columns - use the transposed source matrix as \code{mat_to_cluster}
#'
#' @param mat_to_cluster numeric matrix to be clustered; cluster dimension must be named
#' @param distfun function used to compute the distance (dissimilarity) between rows;
#'   defaults to \code{\link[stats]{dist}} using euclidean euclidean.
#' @param additional_condition additional logical flag whether rows/columns should be clustered 
#'
#' @return logical flag determining if rows should be clustered or \code{hclust} object.
#' 
#' @seealso \code{\link{compute_distances}}
#' 
#' @examples
#' \dontrun{
#' mat <- matrix(1:24, nrow = 4)
#' rownames(mat) <- sprintf("row_%s", 1:4)
#' colnames(mat) <- sprintf("col_%s", 1:6)
#' .get_pheatmap_cluster_param(mat)
#' .get_pheatmap_cluster_param(t(mat))
#' .get_pheatmap_cluster_param(t(mat), distfun = compute_distances)
#' 
#' mat[2,2] <- NA
#' mat[2,1] <- Inf
#' .get_pheatmap_cluster_param(mat)
#' .get_pheatmap_cluster_param(mat, distfun = compute_distances)
#' .get_pheatmap_cluster_param(t(mat), distfun = compute_distances)
#' add_cond <- NCOL(mat) > 10
#' .get_pheatmap_cluster_param(mat, distfun = compute_distances, additional_condition = add_cond)
#' }
#' 
#' @keywords internal
.get_pheatmap_cluster_param <- function(mat_to_cluster,
                                        distfun = stats::dist,
                                        additional_condition = TRUE) {
  
  checkmate::assert_matrix(mat_to_cluster, mode = "numeric", row.names = "unique")
  checkmate::assert_function(distfun)
  checkmate::assert_flag(additional_condition)
  
  if (additional_condition && NROW(mat_to_cluster) >= 2) {
    tryCatch(stats::hclust(distfun(mat_to_cluster)),
             error = function(e) return(FALSE)) # if any problem with distfun -> no clustering
  } else {
    FALSE
  }
}


#' Prep matrix with metric value based on the Metrics assay
#'
#' @param dt_response \code{data.table} representing data from the \code{Metrics} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'  and single-agent \code{SummarizedExperiment} or 
#'  \code{data.table} representing data from the \code{scores} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "scores")}
#'  and combo \code{SummarizedExperiment}
#' @param normalization_type string with normalization types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param metric string name of the metric;
#'  one of: "xc50" ("GR50" or "IC50" - respectively depending on \code{normalization_type}), 
#'  "x_max" ("GR Max" or "E Max") or "x_mean" ("GR Mean" or "RV Mean"); 
#'  but the values from any numeric column can be displayed.
#' @param fit_source string source name for metrics
#' @param experiment_type string with experiment name;
#'                        one of: "single-agent", "combination" or "co-dilution"
#' 
#' @keywords pheat_ann
#' 
#' @return matrix with values for selected metric with \code{CellLinName} in the rows
#'    and  \code{DrugName} (or combination of \code{DrugName} and \code{DrugName_2}) in the columns
#' 
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                                assay_name = "Metrics")
#'
#' mat <- prep_pheatmap_matrix(dt_response = dt_metrics,
#'                             experiment_type = gDRutils::get_supported_experiments("sa"))
#' 
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' dt_scores <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                               assay_name = "scores")
#'
#' mat <- prep_pheatmap_matrix(dt_response = dt_scores,
#'                             metric = "hsa_score",
#'                             normalization_type = "RV",
#'                             experiment_type = gDRutils::get_supported_experiments("combo"))
#'                             
#' @export                           
prep_pheatmap_matrix <- function(dt_response,
                                 normalization_type = "GR",
                                 metric = "xc50",
                                 fit_source = "gDR",
                                 experiment_type = gDRutils::get_supported_experiments("sa")) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  conc_2 <- gDRutils::get_env_identifiers("concentration2")
  
  checkmate::assert_data_table(dt_response)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  numeric_columns <- names(dt_response)[vapply(dt_response, is.numeric, logical(1))]
  checkmate::assert_choice(metric, choices = numeric_columns)
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_choice(experiment_type, choices = gDRutils::get_supported_experiments())
  
  # select data for normalization type
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  tab_response <- dt_response[eval(filter_expr)]
  
  # prep data
  tab_dcast <- if (experiment_type == gDRutils::get_supported_experiments("sa")) {
    data.table::dcast(
      data = tab_response,
      formula = get(cellline_name) ~ paste(get(drug_name)),
      value.var = metric)
  } else if (experiment_type == gDRutils::get_supported_experiments("cd")) {
    data.table::dcast(
      data = tab_response,
      formula = get(cellline_name) ~ paste(get(drug_name), "x", paste0(get(drug_name_2), "__", get(conc_2))),
      value.var = metric)
  } else {
    data.table::dcast(
      data = tab_response,
      formula = get(cellline_name) ~ paste(get(drug_name), "x", get(drug_name_2)),
      value.var = metric)
  }
  data.table::setkey(tab_dcast, NULL)
  data.table::setnames(tab_dcast, "cellline_name", cellline_name)
  
  # prep matrix cellline vs drugs
  mat_cvd <- as.matrix(tab_dcast[, .SD, .SDcols = -cellline_name])
  rownames(mat_cvd) <- tab_dcast[[cellline_name]]
  # remove all-NA-rows and all-NA-columns
  rm_col <- vapply(colnames(mat_cvd), function(i) !all(is.na(mat_cvd[, i])), logical(1))
  rm_row <- vapply(seq_along(rownames(mat_cvd)), function(i) !all(is.na(mat_cvd[i, ])), logical(1))
  if (!all(rm_col)) mat_cvd <- mat_cvd[, rm_col, drop = FALSE]
  if (!all(rm_row)) mat_cvd <- mat_cvd[rm_row, , drop = FALSE]
  
  return(mat_cvd)
}


#' Prep annotation data.table acc to metric matrix for pheatmap::pheatmat
#' 
#' @param dt_anno \code{data.table} that specifies the annotations shown on left side of the heatmap 
#'   or shown above the heatmap - depending on the \code{anno_var}.
#'   Each row defines the features for a specific row. The rows in the data and in the annotation
#'   are matched using corresponding names from the required \code{anno_var} column.
#' @param mat_to_cluster numeric matrix to be clustered; cluster dimension must be named
#' @param anno_var string with variable describing annotation dimension:
#'   one of: \code{CellLineName} for rows or \code{DrugName} for column.
#'
#' @return \code{data.table} with 
#' 
#' @keywords internal
.fill_pheatmap_annotation <- function(dt_anno,
                                      mat_with_metric,
                                      anno_var = gDRutils::get_env_identifiers("cellline_name")) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  
  checkmate::assert_data_table(dt_anno)
  checkmate::assert_matrix(mat_with_metric, mode = "numeric", row.names = "unique", col.names = "unique")
  checkmate::assert_choice(anno_var, choices = c(cellline_name, drug_name))
  
  fun_names <- if (anno_var == cellline_name) rownames else colnames
  
  # check completeness of annotation
  if (!all(fun_names(mat_with_metric) %in% dt_anno[[anno_var]])) {
    tab_missing_ann <- data.table::data.table(
      missing = fun_names(mat_with_metric)[!fun_names(mat_with_metric) %in% dt_anno[[anno_var]]]
    )
    data.table::setnames(tab_missing_ann, "missing", anno_var)
    
    dt_anno <- data.table::rbindlist(list(dt_anno, tab_missing_ann), fill = TRUE)
  }
  
  # foll missing values (note: data.table::nafill does not support character)
  cols <- names(dt_anno)[names(dt_anno) != anno_var]
  data.table::setorderv(dt_anno, cols = cols, na.last = TRUE)
  dt_anno[, (cols) := lapply(.SD, change_NA_into_char, "NA"), .SDcols = cols]
  # select annotation acc to matrix
  dt_anno <- dt_anno[get(anno_var) %in% fun_names(mat_with_metric), ]
  
  return(dt_anno)
}