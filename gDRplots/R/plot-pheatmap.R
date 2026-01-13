#' Plot pretty heatmap for single-agent or combo data to control quality of the data
#'    
#' @param dt_average  \code{data.table} representing data from the \code{Averaged} assay,
#'    outputted by \code{\link[gDRutils:convert_se_assay_to_dt]{gDRutils::convert_se_assay_to_dt}}
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
#' @param no_breaks numeric number of breaks on scale used for mapping values to colors
#' @param cluster_rows logical flag whether rows should be clustered;
#'   the dendrogram will not be shown for the matrix with any dimension greater than 200.
#' @param distfun function used to compute the distance (dissimilarity) between rows;
#'   used for the dendrogram when \code{cluster_rows} is set to TRUE; 
#'   the default is \code{\link{compute_distances}} using Spearman method.
#' @param lbl_by_CellLineName logical flag whether heatmap should be described by CellLineNames instead of clid
#' @param lbl_by_DrugName logical flag whether heatmap should be described by DrugName instead of Gnumber
#' 
#' @seealso \code{\link[pheatmap:pheatmap]{pheatmap::pheatmap}}
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
  # standardization of concentration
  conc_map <- gDRutils::map_conc_to_standardized_conc(tab_response[[conc]], 
                                                      tab_response[[conc_2]])
  tab_response <- merge(tab_response, conc_map, by.x = conc, by.y = "concs")
  tab_response <- merge(tab_response, conc_map, by.x = conc_2, by.y = "concs", suffixes = c("", "_2"))
  
  # prep pivot data
  tab_response <- data.table::setorderv(tab_response, c(gnumber, gnumber_2, "rconcs", "rconcs_2"))
  tab_response$col_pivot_name <- sprintf("%s_%s_%s_%s",
                                         tab_response[[gnumber]],
                                         tab_response[[gnumber_2]],
                                         trimws(format(tab_response[["rconcs"]], scientific = FALSE)),
                                         trimws(format(tab_response[["rconcs_2"]], scientific = FALSE)))
  col_pivot_name <- "col_pivot_name"
  tryCatch({
    fm_string <- "clid ~ col_pivot_name"
    tab_plot <- data.table::dcast(
      data = tab_response,
      formula = stats::as.formula(fm_string),
      value.var = metric
    )
    # data.table  >= 1.17.0
  }, warning = function(w) {
    .stop_on_aggregation("pheatmap_qc", fm_string)
    # data.table  < 1.17.0
  }, message = function(m) {
    .stop_on_aggregation("pheatmap_qc", fm_string)
  })
  data.table::setcolorder(tab_plot, c(clid, unique(tab_response$col_pivot_name)))
  
  # prep matrix
  mat_cvd <- as.matrix(tab_plot[, .SD, .SDcols = -clid])
  rownames(mat_cvd) <- tab_plot[[clid]]
  rm_col <- vapply(colnames(mat_cvd), function(i) !all(is.na(mat_cvd[, i])), logical(1))
  rm_row <- vapply(seq_along(rownames(mat_cvd)), function(i) !all(is.na(mat_cvd[i, ])), logical(1))
  if (!all(sum(rm_col) == 0, sum(rm_col) == 0)) {
    if (!all(rm_col)) mat_cvd <- mat_cvd[, rm_col, drop = FALSE]
    if (!all(rm_row)) mat_cvd <- mat_cvd[rm_row, , drop = FALSE]
  }
  if (metric == "x_std") mat_cvd <- mat_cvd ^ 2
  
  # annotation
  info_drug <- unique(tab_response[col_pivot_name %in% colnames(mat_cvd),
                                   .SD, .SDcols = c("col_pivot_name", gnumber, "rconcs")])
  info_drug_2 <- unique(tab_response[col_pivot_name %in% colnames(mat_cvd),
                                     .SD, .SDcols = c("col_pivot_name", gnumber_2, "rconcs_2")])
  data.table::setnames(info_drug_2, old = c(gnumber_2, "rconcs_2"), new = c(gnumber, "rconcs"))
  info_drug <- rbind(info_drug, info_drug_2)[get(gnumber) != "untreated"]
  tryCatch({
    fm_string <- "col_pivot_name ~ get(gnumber)"
    drug_annotation <- data.table::dcast(
      data = info_drug,
      formula = stats::as.formula(fm_string),
      value.var = "rconcs"
    )
  }, warning = function(w) {
    .stop_on_aggregation("pheatmap_qc", fm_string)
  }, message = function(m) {
    .stop_on_aggregation("pheatmap_qc", fm_string)
  })
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
  min_mat <- if (all(is.na(mat_cvd))) {
    0
  } else {
    round(min(mat_cvd, na.rm = TRUE), digits = 2)
  }
  minval <- min(c(0, min_mat))
  
  breaks <- seq(from = minval, to = maxval, length.out = no_breaks + 1)
  hm_color_palette <-   grDevices::colorRampPalette(colors_vec)(no_breaks)
  if (metric == "x_std") hm_color_palette <- rev(hm_color_palette)
  
  fontsize_row <- .get_pheatmap_fontsize(mat_cvd, "row")
  cellheight <- if (NROW(mat_cvd) == 1) {
    20
  } else {
    NA
  }
  
  hm <- 
    pheatmap::pheatmap(mat = mat_cvd,
                       scale = "none",
                       display_numbers = FALSE,
                       color = hm_color_palette,
                       breaks = breaks,
                       angle_col = 45,
                       show_colnames = FALSE,
                       main = hm_title,
                       fontsize = 8,
                       fontsize_row = fontsize_row,
                       cellheight = cellheight,
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
#' @param dt_metrics \code{data.table} representing data from the \code{"Metrics"} assay,
#'  outputted by \code{\link[gDRutils:convert_se_assay_to_dt]{gDRutils::convert_se_assay_to_dt}}
#'  and single-agent \code{SummarizedExperiment}
#' @param dt_metrics_capped \code{data.table} representing data from the \code{"Metrics"} assay,
#'  the same as \code{dt_metrics} but with capped values with 
#'  \code{\link[gDRutils:cap_assay_infinities]{gDRutils::cap_assay_infinities}}
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
#' @param no_breaks numeric number of breaks on scale used for mapping values to colors
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
#'   the default is \code{\link{compute_distances}} using Spearman method.
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
#' @param max_hm_lbl_length numeric value for the maximum number of characters in the label;
#'   if set to Inf, no trimming will be performed; for better readability, it is recommended to use 
#'   the default number.
#' 
#' @seealso \code{\link[pheatmap:pheatmap]{pheatmap::pheatmap}}
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                                assay_name = "Metrics")
#' dt_averaged <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                                 assay_name = "Averaged")
#' dt_metrics_capped <-
#'   gDRutils::cap_assay_infinities(
#'     conc_assay_dt = dt_averaged,
#'     assay_dt = dt_metrics,
#'     experiment_name = gDRutils::get_supported_experiments("sa"),
#'     col = "xc50",
#'     capping_fold = 5)
#'
#' output <- pheatmap_with_anno_sa(dt_metrics = dt_metrics)
#' hm_0 <- output[["heatmap"]]
#' ggpubr::as_ggplot(hm_0[["gtable"]])
#' 
#' output <- pheatmap_with_anno_sa(dt_metrics = dt_metrics,
#'                                 dt_metrics_capped = dt_metrics_capped)
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
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
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
    annotation_colors = NULL,
    max_hm_lbl_length = 
      gDRutils::get_settings_from_json("MAX_HM_LBL_LENGTH",
                                       system.file(package = "gDRplots", "settings.json"))
) {
  
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
  if (!is.null(annotation_colors)) checkmate::assert_named(annotation_colors)
  checkmate::assert_number(max_hm_lbl_length, lower = 5)
  
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
  
  if (is.null(dt_metrics_capped)) {
    mat_cvd <- mat_cvd_raw
  } else {
    mat_cvd <- prep_pheatmap_matrix(dt_response = dt_metrics_capped,
                                    normalization_type = normalization_type,
                                    metric = metric,
                                    fit_source = fit_source,
                                    experiment_type = gDRutils::get_supported_experiments("sa"))
  }
  
  # check completeness of annotation
  if (!is.null(annotation_col)) {
    annotation_col <- .fill_pheatmap_annotation(dt_anno = annotation_col,
                                                mat_with_metric = mat_cvd,
                                                anno_var = cellline_name)
    ls_output[["data"]][["annotation_col"]] <- annotation_col
    
    rownames(annotation_col) <- annotation_col[[cellline_name]] # required by pheatmap::pheatmap
    # order matrix
    mat_cvd <- mat_cvd[rownames(annotation_col), , drop = FALSE]
    mat_cvd_raw <- mat_cvd_raw[rownames(annotation_col), , drop = FALSE]
    
    # trim annotation
    if (is.finite(max_hm_lbl_length)) {
      ls_too_long_lbl <- 
        vapply(names(annotation_col), function(nm) any(nchar(annotation_col[[nm]]) > max_hm_lbl_length), logical(1))
      if (any(ls_too_long_lbl)) {
        ls_col <- names(annotation_col)[ls_too_long_lbl]
        annotation_col[, (ls_col) := lapply(.SD, .trim_labels, max_lbl_length = max_hm_lbl_length), 
                       .SDcols = ls_col]
        
        rownames(annotation_col) <- annotation_col[[cellline_name]] # update
        rownames(mat_cvd) <- annotation_col[[cellline_name]] # update
      }
    }
    annotation_col <- annotation_col[, .SD, .SDcol = -cellline_name] # required by pheatmap::pheatmap
  }
  
  if (!is.null(annotation_row)) {
    annotation_row <- .fill_pheatmap_annotation(dt_anno = annotation_row,
                                                mat_with_metric = mat_cvd,
                                                anno_var = drug_name)
    ls_output[["data"]][["annotation_row"]] <- annotation_row
    
    rownames(annotation_row) <- annotation_row[[drug_name]] # required by pheatmap::pheatmap
    # order matrix
    mat_cvd <- mat_cvd[, rownames(annotation_row), drop = FALSE]
    mat_cvd_raw <- mat_cvd_raw[, rownames(annotation_row), drop = FALSE]
    
    # trim annotation
    if (is.finite(max_hm_lbl_length)) {
      ls_too_long_lbl <- 
        vapply(names(annotation_row), function(nm) any(nchar(annotation_row[[nm]]) > max_hm_lbl_length), logical(1))
      if (any(ls_too_long_lbl)) {
        ls_col <- names(annotation_row)[ls_too_long_lbl]
        annotation_row[, (ls_col) := lapply(.SD, .trim_labels, max_lbl_length = max_hm_lbl_length), 
                       .SDcols = ls_col]
        
        rownames(annotation_row) <- annotation_row[[drug_name]] # update
        colnames(mat_cvd) <- annotation_row[[drug_name]] # update
      }
    }
    
    annotation_row <- annotation_row[, .SD, .SDcol = -drug_name] # required by pheatmap::pheatmap
  }
  
  if (!is.null(annotation_colors)) {
    # trim annotation
    if (is.finite(max_hm_lbl_length)) {
      ls_too_long_lbl <- 
        vapply(names(annotation_colors), function(nm) {
          any(nchar(names(annotation_colors[[nm]])) > max_hm_lbl_length) }, logical(1))
      if (any(ls_too_long_lbl) && is.finite(max_hm_lbl_length)) {
        for (nm in names(annotation_colors)[ls_too_long_lbl]) {
          names(annotation_colors[[nm]]) <- 
            .trim_labels(lbls_vec = names(annotation_colors[[nm]]), 
                         max_lbl_length = max_hm_lbl_length)
        }
      }
    }
    # filling missing values
    if (!is.null(annotation_col)) {
      annotation_colors <- fill_ann_color_map(dt_ann = annotation_col,
                                              map_ann = annotation_colors)
    }
    if (!is.null(annotation_row)) {
      annotation_colors <- fill_ann_color_map(dt_ann = annotation_row, 
                                              map_ann = annotation_colors)
    }
  }
  
  ls_output[["data"]][["matrix"]] <- 
    data.table::as.data.table(mat_cvd_raw, keep.rownames = cellline_name)
  
  # trim colnames & rownames for matrix
  if (is.finite(max_hm_lbl_length)) {
    if (any(nchar(colnames(mat_cvd)) > max_hm_lbl_length)) {
      colnames(mat_cvd) <- .trim_labels(lbls_vec = colnames(mat_cvd), 
                                        max_lbl_length = max_hm_lbl_length)
    }
    if (any(nchar(rownames(mat_cvd)) > max_hm_lbl_length)) {
      rownames(mat_cvd) <- .trim_labels(lbls_vec = rownames(mat_cvd), 
                                        max_lbl_length = max_hm_lbl_length)
    }
  }
  
  # flip
  t_mat_cvd <- t(mat_cvd)
  t_mat_cvd[] <- vapply(t_mat_cvd, function(x) qmfun(x), numeric(1))
  
  # dendrogram
  max_dim_matrix_cluster <- 
    gDRutils::get_settings_from_json("MAX_DIM_MATRIX_CLUSTER",
                                     system.file(package = "gDRplots", "settings.json"))
  gDR_cluster_condition <- any(dim(t_mat_cvd) < max_dim_matrix_cluster)  # gDR standard
  if (cluster_rows) {
    cluster_rows <- 
      .get_pheatmap_cluster_param(mat_to_cluster = t_mat_cvd,
                                  distfun = distfun,
                                  additional_condition = gDR_cluster_condition)
  }
  if (cluster_cols) {
    cluster_cols <- 
      .get_pheatmap_cluster_param(mat_to_cluster = t(t_mat_cvd),
                                  distfun = distfun,
                                  additional_condition = gDR_cluster_condition)
  }
  
  # prep hm color palette
  min_val <- if (all(is.infinite(t_mat_cvd) | is.na(t_mat_cvd))) {
    -1
  } else {
    min(t_mat_cvd[!is.infinite(t_mat_cvd)], na.rm = TRUE)
  }
  max_val <- if (metric %in% c("x", "xc50", "x_max", "x_mean")) {
    1.0
  } else {
    max(t_mat_cvd, na.rm = TRUE)
  }  
  if (min_val == max_val) {
    min_val <- min_val - 1
  }
  
  breaks <- seq(from = min_val, to = max_val, length.out = no_breaks + 1)
  hm_color_palette <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
    .get_smooth_palette(no_breaks)
  } else {
    grDevices::colorRampPalette(colors_vec)(no_breaks)
  }
  
  # display numbers - for readability, turn it off for matrices larger than 15x15
  display_numbers_flag <- !any(dim(t_mat_cvd) > c(15, 15))
  number_color <- 
    if (any(vapply(hm_color_palette, is_color_dark, logical(1))) && display_numbers_flag) {
      mat_ <- .get_pheatmap_number_color(mat_with_metric = t_mat_cvd,
                                         colors_vec = hm_color_palette,
                                         breaks = breaks)
      if (inherits(cluster_rows, "hclust")) mat_ <- mat_[cluster_rows$order, ]
      if (inherits(cluster_cols, "hclust")) mat_ <- mat_[, cluster_cols$order]
      mat_
    } else {
      "black"
    }
  
  fontsize_col <- .get_pheatmap_fontsize(t_mat_cvd, "col")
  
  ls_output[["heatmap"]] <- 
    pheatmap::pheatmap(mat = t_mat_cvd,
                       scale = "none",
                       display_numbers = display_numbers_flag,
                       number_color = number_color,
                       fontsize_number = 8,
                       color = hm_color_palette,
                       breaks = breaks,
                       border_color = "lightgrey",
                       angle_col = 90,
                       main = hm_title,
                       fontsize = 8,
                       fontsize_col = fontsize_col,
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
#'  outputted by  \code{\link[gDRutils:convert_se_assay_to_dt]{gDRutils::convert_se_assay_to_dt}}
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
#' @param no_breaks numeric number of breaks on scale used for mapping values to colors
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
#'   the default is \code{\link{compute_distances}} using Spearman method.
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
#' @seealso \code{\link[pheatmap:pheatmap]{pheatmap::pheatmap}}
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
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
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
  
  # check completeness of annotation
  if (!is.null(annotation_col)) {
    annotation_col <- .fill_pheatmap_annotation(dt_anno = annotation_col,
                                                mat_with_metric = mat_cvd,
                                                anno_var = cellline_name)
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
    annotation_colors <- fill_ann_color_map(dt_ann = annotation_col, 
                                            map_ann = annotation_colors)
  }
  if (!is.null(annotation_row) && !is.null(annotation_colors)) {
    annotation_colors <- fill_ann_color_map(dt_ann = annotation_row, 
                                            map_ann = annotation_colors)
  }
  
  ls_output[["data"]][["matrix"]] <- 
    data.table::as.data.table(mat_cvd, keep.rownames = cellline_name)
  
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
  max_val <- if (metric %in% c("x", "xc50", "x_max", "x_mean")) {
    1.0
  } else {
    max(t_mat_cvd, na.rm = TRUE)
  }  
  if (min_val == max_val) {
    min_val <- min_val - 1
  } else if (is.infinite(min_val)) {
    min_val <- max_val / 100
  }
  
  breaks <- seq(from = min_val, to = max_val, length.out = no_breaks + 1)
  hm_color_palette <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
    .get_smooth_palette(no_breaks)
  } else {
    grDevices::colorRampPalette(colors_vec)(no_breaks)
  }
  
  # display numbers - for readability, turn it off for matrices larger than 15x15
  display_numbers_flag <- !any(dim(t_mat_cvd) > c(15, 15))
  number_color <- 
    if (any(vapply(hm_color_palette, is_color_dark, logical(1))) && display_numbers_flag) {
      mat_ <- .get_pheatmap_number_color(mat_with_metric = t_mat_cvd,
                                         colors_vec = hm_color_palette,
                                         breaks = breaks)
      if (inherits(cluster_rows, "hclust")) mat_ <- mat_[cluster_rows$order, ]
      if (inherits(cluster_cols, "hclust")) mat_ <- mat_[, cluster_cols$order]
      mat_
    } else {
      "black"
    }
  
  fontsize_col <- .get_pheatmap_fontsize(t_mat_cvd, "col")
  
  ls_output[["heatmap"]] <- 
    pheatmap::pheatmap(mat = t_mat_cvd,
                       scale = "none",
                       display_numbers = display_numbers_flag,
                       number_color = number_color,
                       fontsize_number = 8,
                       color = hm_color_palette,
                       breaks = breaks,
                       border_color = "lightgrey",
                       angle_col = 90,
                       main = hm_title,
                       fontsize = 8,
                       fontsize_col = fontsize_col,
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
#'   outputted by \code{\link[gDRutils:convert_se_assay_to_dt]{gDRutils::convert_se_assay_to_dt}}
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
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
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
    annotation_colors = NULL,
    max_hm_lbl_length = 
      gDRutils::get_settings_from_json("MAX_HM_LBL_LENGTH",
                                       system.file(package = "gDRplots", "settings.json"))
) {
  
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
  checkmate::assert_number(max_hm_lbl_length, lower = 5)
  
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
  mat_cvd_raw <- mat_cvd
  
  # edge-case (no valid data in the matrix, usually matrix with NAs only)  
  if (NROW(mat_cvd) == 0) {
    return(ls_output)
  }
  
  # check completeness of annotation
  if (!is.null(annotation_col)) {
    annotation_col <- .fill_pheatmap_annotation(dt_anno = annotation_col,
                                                mat_with_metric = mat_cvd,
                                                anno_var = cellline_name)
    ls_output[["data"]][["annotation_col"]] <- annotation_col
    
    rownames(annotation_col) <- annotation_col[[cellline_name]] # required by pheatmap::pheatmap
    # order matrix
    mat_cvd <- mat_cvd[rownames(annotation_col), , drop = FALSE]
    mat_cvd_raw <- mat_cvd_raw[rownames(annotation_col), , drop = FALSE]
    
    # trim annotation
    if (is.finite(max_hm_lbl_length)) {
      ls_too_long_lbl <- 
        vapply(names(annotation_col), function(nm) any(nchar(annotation_col[[nm]]) > max_hm_lbl_length), logical(1))
      if (any(ls_too_long_lbl)) {
        ls_col <- names(annotation_col)[ls_too_long_lbl]
        annotation_col[, (ls_col) := lapply(.SD, .trim_labels, max_lbl_length = max_hm_lbl_length), 
                       .SDcols = ls_col]
        
        rownames(annotation_col) <- annotation_col[[cellline_name]] # update
        rownames(mat_cvd) <- annotation_col[[cellline_name]] # update
      }
    }
    
    annotation_col <- annotation_col[, .SD, .SDcol = -cellline_name]
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
    # order matrix
    mat_cvd <- mat_cvd[, rownames(annotation_row), drop = FALSE]
    mat_cvd_raw <- mat_cvd_raw[, rownames(annotation_row), drop = FALSE]
    
    # trim annotation
    if (is.finite(max_hm_lbl_length)) {
      ls_too_long_lbl <- 
        vapply(names(annotation_row), function(nm) any(nchar(annotation_row[[nm]]) > max_hm_lbl_length), logical(1))
      if (any(ls_too_long_lbl)) {
        ls_col <- names(annotation_row)[ls_too_long_lbl]
        annotation_row[, (ls_col) := lapply(.SD, .trim_labels, max_lbl_length = max_hm_lbl_length), 
                       .SDcols = ls_col]
        
        if (any(ls_col %in% c(drug_name, drug_name_2))) {
          annotation_row$DrugCombination <-
            paste(annotation_row[[drug_name]], "x", annotation_row[[drug_name_2]])
        }
        
        rownames(annotation_row) <- annotation_row[["DrugCombination"]] # update
        colnames(mat_cvd) <- annotation_row[["DrugCombination"]] # update
      }
    }
    
    annotation_row <- annotation_row[, .SD, .SDcol = -c(drug_name, drug_name_2, "DrugCombination")]
  }
  
  if (!is.null(annotation_colors)) {
    # trim annotation
    if (is.finite(max_hm_lbl_length)) {
      ls_too_long_lbl <- 
        vapply(names(annotation_colors), function(nm) {
          any(nchar(names(annotation_colors[[nm]])) > max_hm_lbl_length) }, logical(1))
      if (any(ls_too_long_lbl) && is.finite(max_hm_lbl_length)) {
        for (nm in names(annotation_colors)[ls_too_long_lbl]) {
          names(annotation_colors[[nm]]) <- 
            .trim_labels(lbls_vec = names(annotation_colors[[nm]]), 
                         max_lbl_length = max_hm_lbl_length)
        }
      }
    }
    
    # filling missing values
    if (!is.null(annotation_col)) {
      annotation_colors <- fill_ann_color_map(dt_ann = annotation_col,
                                              map_ann = annotation_colors)
    }
    if (!is.null(annotation_row)) {
      annotation_colors <- fill_ann_color_map(dt_ann = annotation_row, 
                                              map_ann = annotation_colors)
    }
  }
  
  ls_output[["data"]][["matrix"]] <- 
    data.table::as.data.table(mat_cvd_raw, keep.rownames = cellline_name)
  # trim colnames & rownames for matrix
  if (is.finite(max_hm_lbl_length)) {
    if (any(nchar(colnames(mat_cvd)) > max_hm_lbl_length)) {
      dt_lbls <- data.table::data.table(src = colnames(mat_cvd))
      drug_1 <- drug_2 <- trimmed <- NULL
      dt_lbls[, c("drug_1", "drug_2") := data.table::tstrsplit(src, " x ", fixed = TRUE)]
      dt_lbls[, c("drug_1", "drug_2") := lapply(.SD, .trim_labels, max_lbl_length = max_hm_lbl_length), 
              .SDcols = c("drug_1", "drug_2")]
      dt_lbls[, trimmed := paste(drug_1, drug_2, sep = " x ")]
      colnames(mat_cvd) <- dt_lbls[["trimmed"]]
    }
    if (any(nchar(rownames(mat_cvd)) > max_hm_lbl_length)) {
      rownames(mat_cvd) <- .trim_labels(lbls_vec = rownames(mat_cvd),
                                        max_lbl_length = max_hm_lbl_length)
    }
  }
  
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
  breaks <- seq(from = -0.7, to = 0.7, length.out = no_breaks + 1)
  hm_color_palette <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
    .get_excess_palette(no_breaks)
  } else {
    grDevices::colorRampPalette(colors_vec)(no_breaks)
  }
  
  # display numbers - for readability, turn it off for matrices larger than 15x15
  display_numbers_flag <- !any(dim(t_mat_cvd) > c(15, 15))
  number_color <- 
    if (any(vapply(hm_color_palette, is_color_dark, logical(1))) && display_numbers_flag) {
      mat_ <- .get_pheatmap_number_color(mat_with_metric = t_mat_cvd,
                                         colors_vec = hm_color_palette,
                                         breaks = breaks)
      if (inherits(cluster_rows, "hclust")) mat_ <- mat_[cluster_rows$order, ]
      if (inherits(cluster_cols, "hclust")) mat_ <- mat_[, cluster_cols$order]
      mat_
    } else {
      "black"
    }
  
  fontsize_col <- .get_pheatmap_fontsize(t_mat_cvd, "col")
  
  ls_output[["heatmap"]] <- 
    pheatmap::pheatmap(t_mat_cvd,
                       scale = "none",
                       display_numbers = display_numbers_flag,
                       number_color = "black",
                       fontsize_number = 8,
                       color = hm_color_palette,
                       breaks = breaks,
                       border_color = "lightgrey",
                       angle_col = 90,
                       main = hm_title,
                       fontsize = 8,
                       fontsize_col = fontsize_col,
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

#' Plot heatmap with annotations for combo metrics data
#' 
#' Plots a heatmap for metrics data (e.g., "xc50", "x_inf") covering both Single Agents and Combinations.
#' Rows are sorted alphabetically by the primary drug, prioritizing Single Agent conditions first, 
#' followed by combinations sorted by co-treatment name and concentration.
#'
#' @param dt_metrics \code{data.table} representing data from the \code{metrics} assay, 
#'    outputted by \code{\link[gDRutils:convert_se_assay_to_dt]{gDRutils::convert_se_assay_to_dt}}.
#'    Must contain \code{DrugName}, \code{DrugName_2}, \code{cotrt_value}, and \code{CellLineName}.
#' @param dt_metrics_capped \code{data.table} (optional) representing capped data from the \code{metrics} assay;
#'    if provided, it is used instead of \code{dt_metrics}.
#' @param normalization_type string representing the normalization type;
#'    one of: "GR", "RV".
#' @param metric string name of the metric to plot;
#'    e.g., "xc50", "x_inf", "ec50", "r2", "h_RV", "h_GR".
#' @param fit_source string representing the source of the fit; 
#'    defaults to "gDR".
#' @param annotation_row \code{data.table} (optional) that specifies the annotations shown on
#'    the left side of the heatmap.
#'    \emph{Note: In this specific function, row annotations are largely auto-generated based on the 
#'    co-treatment (Fixed Drug) names and concentrations to ensure consistent formatting.}
#' @param annotation_col \code{data.table} that specifies the annotations shown on the top of the heatmap.
#'    The rows in the annotation are matched using the \code{CellLineName} column.
#' @param annotation_colors named list for specifying \code{annotation_col} and \code{annotation_row} 
#'    track colors manually; note list is named with annotation name (column names of \code{annotation_row} 
#'    and \code{annotation_col}), and each list item is a named vector with a valid color name for 
#'    each value described in the annotations. Not described elements will be colored by default.
#' @inheritParams pheatmap_with_anno_sa
#' @export
pheatmap_with_anno_combo_metrics <- function(
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
    annotation_colors = NULL,
    max_hm_lbl_length = gDRutils::get_settings_from_json("MAX_HM_LBL_LENGTH",
                                                         system.file(package = "gDRplots",
                                                                     "settings.json"))
) {
  
  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")[1]
  drug_name_id <- gDRutils::get_env_identifiers("drug_name")
  drug_name_id2 <- gDRutils::get_env_identifiers("drug_name2")
  cellline_name_id <- gDRutils::get_env_identifiers("cellline_name")
  
  conc3 <- gDRutils::get_env_identifiers("concentration3")
  drug3 <- gDRutils::get_env_identifiers("drug3")
  drug_name_id3 <- gDRutils::get_env_identifiers("drug_name3")
  
  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_data_table(dt_metrics_capped, null.ok = TRUE)
  
  if (!is.null(dt_metrics_capped)) {
    dt_to_use <- dt_metrics_capped
  } else {
    dt_to_use <- dt_metrics
  }
  
  req_cols <- c(drug_name_id, drug_name_id2, "cotrt_value", cellline_name_id, 
                "normalization_type", "fit_source")
  
  checkmate::assert_names(names(dt_to_use), must.include = req_cols)
  
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_sub <- dt_to_use[eval(filter_expr)]
  
  checkmate::assert_data_table(dt_sub, min.rows = 1, .var.name = "dt_sub (filtered data)")
  
  dt_sub <- data.table::copy(dt_sub)
  
  pattern <- paste0("\\s*\\([^)]*(", drug3, "|", conc3, ")[^)]*\\)")
  dt_sub[[drug_name_id]] <- trimws(gsub(pattern, "", dt_sub[[drug_name_id]]))
  
  if ("dilution_drug" %in% names(dt_sub)) {
    dt_sub <- dt_sub[dilution_drug != "codilution"]
  }
  
  dt_sub[is.na(cotrt_value), cotrt_value := 0]
  
  dt_sub[, `:=`(
    Row_Display_Name = get(drug_name_id),
    Fixed_Name_1     = get(drug_name_id2),
    Fixed_Conc_1     = as.numeric(cotrt_value),
    Fixed_Name_2     = untreated_tag,
    Fixed_Conc_2     = 0
  )]
  
  dt_sub[Fixed_Conc_1 == 0, Fixed_Name_1 := untreated_tag]
  
  if (drug_name_id3 %in% names(dt_sub) && conc3 %in% names(dt_sub)) {
    dt_sub[is.na(get(conc3)), (conc3) := 0]
    
    dt_sub[get(conc3) > 0, `:=`(
      Fixed_Name_2 = get(drug_name_id3),
      Fixed_Conc_2 = as.numeric(get(conc3))
    )]
  }
  
  dt_sub[, Fixed_Conc_1 := .round_to_unique_string(Fixed_Conc_1)]
  dt_sub[, Fixed_Conc_2 := .round_to_unique_string(Fixed_Conc_2)]
  
  dt_sub[, Fixed_Label_1 := ifelse(Fixed_Name_1 == untreated_tag, untreated_tag, 
                                   paste0(Fixed_Name_1, " (", Fixed_Conc_1, ")"))]
  
  dt_sub[, Fixed_Label_2 := ifelse(Fixed_Name_2 == untreated_tag, untreated_tag, 
                                   paste0(Fixed_Name_2, " (", Fixed_Conc_2, ")"))]
  
  dt_sub[, Treatment_Key := paste(Row_Display_Name, Fixed_Label_1, Fixed_Label_2, sep = "__")]
  
  fm_string <- paste("Treatment_Key ~", cellline_name_id)
  
  tryCatch({
    tab_dcast <- data.table::dcast(dt_sub, 
                                   formula = stats::as.formula(fm_string), 
                                   value.var = metric, 
                                   fun.aggregate = mean, 
                                   na.rm = TRUE)
  }, error = function(e) {
    .stop_on_aggregation(e)
  })
  
  mat_cvd <- as.matrix(tab_dcast[, -1, with = FALSE])
  rownames(mat_cvd) <- tab_dcast[[1]]
  
  rm_col <- vapply(colnames(mat_cvd), function(i) !all(is.na(mat_cvd[, i])), logical(1))
  rm_row <- vapply(seq_along(rownames(mat_cvd)), function(i) !all(is.na(mat_cvd[i, ])), logical(1))
  
  if (!all(rm_col)) {
    mat_cvd <- mat_cvd[, rm_col, drop = FALSE]
  }
  if (!all(rm_row)) {
    mat_cvd <- mat_cvd[rm_row, , drop = FALSE]
  }
  
  mat_cvd_raw <- mat_cvd
  
  qmfun <- switch(metric, "xc50" = log10, identity)
  mat_cvd[] <- vapply(mat_cvd, function(x) qmfun(x), numeric(1))
  
  ls_output <- list(data = list(matrix = NULL, annotation_col = NULL, annotation_row = NULL), heatmap = NULL)
  
  if (!is.null(annotation_col)) {
    annotation_col <- .fill_pheatmap_annotation(annotation_col, t(mat_cvd), cellline_name_id)
    ls_output[["data"]][["annotation_col"]] <- annotation_col
    rownames(annotation_col) <- annotation_col[[cellline_name_id]]
    annotation_col <- annotation_col[, .SD, .SDcol = -cellline_name_id]
    
    mat_cvd <- mat_cvd[, rownames(annotation_col), drop = FALSE]
    mat_cvd_raw <- mat_cvd_raw[, rownames(annotation_col), drop = FALSE]
  }

  row_meta <- unique(dt_sub[Treatment_Key %in% rownames(mat_cvd), 
                            .(Treatment_Key, Row_Display_Name, 
                              Fixed_Label_1, Fixed_Label_2,
                              Fixed_Name_1, Fixed_Conc_1, 
                              Fixed_Name_2, Fixed_Conc_2)])
  
  ord <- order(
    row_meta$Row_Display_Name,
    row_meta$Fixed_Name_1 != untreated_tag,
    row_meta$Fixed_Name_1,
    as.numeric(row_meta$Fixed_Conc_1),
    row_meta$Fixed_Name_2 != untreated_tag, 
    row_meta$Fixed_Name_2,
    as.numeric(row_meta$Fixed_Conc_2)
  )
  
  row_meta <- row_meta[ord, ]
  
  mat_cvd <- mat_cvd[row_meta$Treatment_Key, , drop = FALSE]
  mat_cvd_raw <- mat_cvd_raw[row_meta$Treatment_Key, , drop = FALSE]
  
  anno_df <- data.frame(
    `Fixed_Drug` = row_meta$Fixed_Label_1,
    stringsAsFactors = FALSE
  )
  
  if (any(row_meta$Fixed_Label_2 != untreated_tag)) {
    anno_df$`Fixed_Drug_2` <- row_meta$Fixed_Label_2
  }
  
  rownames(anno_df) <- row_meta$Treatment_Key
  ls_output[["data"]][["annotation_row"]] <- anno_df
  
  if (is.null(annotation_colors)) {
    annotation_colors <- list()
  }
  if (!is.null(annotation_col)) {
    annotation_colors <- fill_ann_color_map(annotation_col, annotation_colors)
  }
  
  for (nm in names(anno_df)) {
    levels_vec <- unique(anno_df[[nm]])
    
    if (length(levels_vec) > 0) {
      new_colors <- get_qual_colors(length(levels_vec))
      names(new_colors) <- levels_vec
      
      if (untreated_tag %in% names(new_colors)) {
        new_colors[untreated_tag] <- "#E0E0E0"
      }
      annotation_colors[[nm]] <- new_colors
    }
  }
  
  row_labels_display <- row_meta$Row_Display_Name
  
  if (is.finite(max_hm_lbl_length)) {
    if (any(nchar(colnames(mat_cvd)) > max_hm_lbl_length)) {
      colnames(mat_cvd) <- .trim_labels(colnames(mat_cvd), max_hm_lbl_length)
    }
    if (any(nchar(row_labels_display) > max_hm_lbl_length)) {
      row_labels_display <- .trim_labels(row_labels_display, max_hm_lbl_length)
    }
  }
  
  max_dim <- gDRutils::get_settings_from_json("MAX_DIM_MATRIX_CLUSTER",
                                              system.file(package = "gDRplots", "settings.json"))
  gDR_cluster_condition <- any(dim(mat_cvd) < max_dim)
  
  if (cluster_rows) {
    cl_rows <- .get_pheatmap_cluster_param(mat_cvd, distfun, gDR_cluster_condition)
  } else {
    cl_rows <- FALSE
  }
  
  if (cluster_cols) {
    cl_cols <- .get_pheatmap_cluster_param(t(mat_cvd), distfun, gDR_cluster_condition)
  } else {
    cl_cols <- FALSE
  }
  
  min_val <- min(mat_cvd[!is.infinite(mat_cvd)], na.rm = TRUE)
  max_val <- max(mat_cvd[!is.infinite(mat_cvd)], na.rm = TRUE)
  
  if (is.infinite(min_val)) {
    min_val <- 0
  }
  if (is.infinite(max_val)) {
    max_val <- 0
  }
  if (min_val == max_val) {
    max_val <- max_val + 0.1
    min_val <- min_val - 0.1
  }
  
  breaks <- seq(min_val, max_val, length.out = no_breaks + 1)
  
  if (is.null(colors_vec)) {
    hm_colors <- .get_smooth_palette(no_breaks)
  } else {
    hm_colors <- grDevices::colorRampPalette(colors_vec)(no_breaks)
  }
  
  ls_output[["data"]][["matrix"]] <- data.table::as.data.table(mat_cvd_raw,
                                                               keep.rownames = "Treatment_Key")
  
  if (is.na(hm_title)) {
    hm_title <- paste0("Combination Metrics: ", metric)
  }
  
  ls_output[["heatmap"]] <- 
    pheatmap::pheatmap(
      mat = mat_cvd,
      scale = "none",
      display_numbers = !any(dim(mat_cvd) > c(30, 30)),
      color = hm_colors,
      breaks = breaks,
      border_color = "lightgrey",
      angle_col = 90,
      main = hm_title,
      fontsize = 8,
      fontsize_col = .get_pheatmap_fontsize(mat_cvd, "col"),
      fontsize_row = .get_pheatmap_fontsize(mat_cvd, "row"),
      na_col = "darkgray",
      cluster_rows = cl_rows,
      cluster_cols = cl_cols,
      annotation_row = anno_df,
      annotation_col = annotation_col,
      annotation_colors = annotation_colors,
      labels_row = row_labels_display,
      silent = TRUE
    )
  
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


#' Prep matrix with metric value based on the Metrics assay
#'
#' @param dt_response \code{data.table} representing data from the \code{Metrics} assay,
#'  outputted by \code{\link[gDRutils:convert_se_assay_to_dt]{gDRutils::convert_se_assay_to_dt}}
#'  and single-agent \code{SummarizedExperiment} or 
#'  \code{data.table} representing data from the \code{scores} assay,
#'  outputted by \code{\link[gDRutils:convert_se_assay_to_dt]{gDRutils::convert_se_assay_to_dt}}
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
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
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
    tryCatch({
      fm_string <- "get(cellline_name) ~ paste(get(drug_name))"
      data.table::dcast(
        data = tab_response,
        formula = stats::as.formula(fm_string),
        value.var = metric
      )
    }, warning = function(w) {
      .stop_on_aggregation("prep_pheatmap_matrix", fm_string)
    }, message = function(m) {
      .stop_on_aggregation("prep_pheatmap_matrix", fm_string)
    })
    
  } else if (experiment_type == gDRutils::get_supported_experiments("cd")) {
    tryCatch({
      fm_string <- 'get(cellline_name) ~ paste(get(drug_name), "x", paste0(get(drug_name_2), "__", get(conc_2)))'
      data.table::dcast(
        data = tab_response,
        formula = stats::as.formula(fm_string),
        value.var = metric
      )
    }, warning = function(w) {
      .stop_on_aggregation("prep_pheatmap_matrix", fm_string)
    }, message = function(m) {
      .stop_on_aggregation("prep_pheatmap_matrix", fm_string)
    })
  } else {
    tryCatch({
      fm_string <- 'get(cellline_name) ~ paste(get(drug_name), "x", get(drug_name_2))'
      data.table::dcast(
        data = tab_response,
        formula = stats::as.formula(fm_string),
        value.var = metric
      )
    }, warning = function(w) {
      .stop_on_aggregation("prep_pheatmap_matrix", fm_string)
    }, message = function(m) {
      .stop_on_aggregation("prep_pheatmap_matrix", fm_string)
    })
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


#' Change NA into given string
#'
#' @param x vector with items suspected of being NA
#' @param lbl_NA string - replacement for NA - as default "NA"
#'
#' @return character (for NA -> given string)
#' @keywords internal
#' 
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#' 
#' @export 
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
#' 
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#' 
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
#' 
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#' 
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
        col_na <- if (any(map_ann[[ann]] %in% c("black", "#000000"))) {
          "darkred"
        } else {
          "black"
        }
        map_ann[[ann]] <- c(map_ann[[ann]], "NA" = col_na)
      } else if (NROW(missing_lvl) > 0) {
        map_ann[[ann]] <- NULL # allow default coloring
      }
      map_ann[[ann]] <- map_ann[[ann]][required_lvl]
    }
  }
  # final
  map_ann
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
#'   defaults to \code{\link[stats:dist]{stats::dist}} using euclidean euclidean.
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
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
.get_pheatmap_cluster_param <- function(mat_to_cluster,
                                        distfun = stats::dist,
                                        additional_condition = TRUE) {
  
  checkmate::assert_matrix(mat_to_cluster, mode = "numeric", row.names = "unique")
  checkmate::assert_function(distfun)
  checkmate::assert_flag(additional_condition)
  
  if (additional_condition && NROW(mat_to_cluster) >= 2) {
    tryCatch(stats::hclust(distfun(mat_to_cluster)),
             error = function(e) { 
               FALSE 
             }) # if any problem with distfun -> no clustering
  } else {
    FALSE
  }
}


#' Prep annotation data.table acc to metric matrix for pheatmap::pheatmat
#' 
#' @param dt_anno \code{data.table} that specifies the annotations shown on left side of the heatmap 
#'   or shown above the heatmap - depending on the \code{anno_var}.
#'   Each row defines the features for a specific row. The rows in the data and in the annotation
#'   are matched using corresponding names from the required \code{anno_var} column.
#' @param mat_with_metric numeric matrix with metric values; must have named rows and columns
#' @param anno_var string with variable describing annotation dimension:
#'   one of: \code{CellLineName} for rows or \code{DrugName} for column.
#'
#' @return \code{data.table} with annotation updated to \code{mat_with_metric}
#' 
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#' 
#' @keywords internal
.fill_pheatmap_annotation <- function(
    dt_anno,
    mat_with_metric,
    anno_var = gDRutils::get_env_identifiers("cellline_name")) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  
  checkmate::assert_data_table(dt_anno)
  checkmate::assert_matrix(mat_with_metric, mode = "numeric", row.names = "unique", col.names = "unique")
  checkmate::assert_choice(anno_var, choices = c(cellline_name, drug_name))
  checkmate::assert_subset(anno_var, choices = names(dt_anno))
  
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

#' stop wrapper for `data.table::dcast` to handle unexpected aggregation
#' 
#' idea from: https://github.com/Rdatatable/data.table/issues/5386
#' 
#' @param fname string with the name of the function that failed to the \code{data.table::dcast} aggregation
#' @param formula string with the formula used in \code{data.table::dcast}
#' 
#' @keywords internal
#' 
#' @author Arkadiusz Gladki \email{arkadiusz.gladki@@contractors.roche.com}
#'
#' @return \code{NULL}
.stop_on_aggregation <- function(fname, formula) {
  checkmate::assert_string(fname)
  checkmate::assert_string(formula)
  
  stop(
    sprintf(
      "Unexpected data aggregation in function: '%s' with formula: '%s'",
      fname,
      formula
    )
  )
  invisible(NULL)
}

#' Compute color for number font in pheatmap::pheatmap based on given color palette and breaks
#'
#' @param mat_with_metric numeric matrix with metric values; must have named rows and columns
#' @param colors_vec character vector of colors (valid name or hex) used in heatmap;
#'   must to be one item shorter than \code{no_breaks}
#' @param breaks numeric vector of breaks on scale used for mapping values to colors
#' @param dark_color_font string with valid color name of font for field with dark background
#' @param light_color_font string with valid color name of font for field without dark background
#'
#' @return named \code{matrix} with number color
#' @examples
#' \dontrun{
#' mat <- matrix(-14:30, ncol = 5,
#'               dimnames = list(letters[1:9], LETTERS[1:5]))
#' no_breaks <- 15
#' breaks <- seq(from = min(mat), to = max(mat), length.out = no_breaks + 1)
#' ls_colors <- c("limegreen", "darkblue", "orange")
#' hm_colors <- grDevices::colorRampPalette(ls_colors)(no_breaks)
#' 
#' number_color <- .get_pheatmap_number_color(mat, hm_colors, breaks)
#' 
#' pheatmap::pheatmap(mat,
#'                    breaks = breaks,
#'                    color = hm_colors,
#'                    display_numbers = TRUE,
#'                    number_color = number_color,
#'                    cluster_rows = FALSE,
#'                    cluster_cols = FALSE)
#' }
#' 
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#' 
#' @keywords internal
.get_pheatmap_number_color <- function(mat_with_metric,
                                       colors_vec,
                                       breaks,
                                       dark_color_font = "white",
                                       light_color_font = "black") {
  
  checkmate::assert_matrix(mat_with_metric, mode = "numeric", row.names = "unique", col.names = "unique")
  checkmate::assert_numeric(breaks, any.missing = FALSE, min.len = 3)
  checkmate::assert_character(colors_vec, len = NROW(breaks) - 1) # required by pheatmap::pheatmap
  checkmate::assert_string(dark_color_font)
  checkmate::assert_string(light_color_font)
  
  # preserve some buggy name for font
  if (!is_valid_color(dark_color_font)) dark_color_font <- "white"
  if (!is_valid_color(light_color_font)) light_color_font <- "black"
  
  # fast end end when colors_vec does not contain valid color names or contain not dark colors only
  if (!all(vapply(colors_vec, is_valid_color, logical(1))) ||
      !any(vapply(colors_vec, is_color_dark, logical(1)))) return(light_color_font)
  
  # breaks must be arranged in ascending order
  breaks <- sort(breaks)
  
  # check dark colors
  ls_dark_area <- which(vapply(colors_vec, is_color_dark, logical(1)))
  # dark color stars...
  first <- min(ls_dark_area)
  last <- max(ls_dark_area) # dark color ends.
  # if the dark colors are in the middle of palette 
  middle <- if (any(diff(ls_dark_area) > 1)) {
    which(diff(ls_dark_area) > 1)
  } else {
    NULL
  }
  # final dark ranges (index of colors in hm_color_palette)
  # final dark ranges (index of colors in hm_color_palette)
  # colors_vec:  | color_1 | color_2 | color_3 | ... | color_n |
  # breaks:      1         2         3         4     n        n+1
  dark_idx <- c(first, ls_dark_area[middle] + 1, ls_dark_area[middle + 1], last + 1)
  # dark breaks (numeric value for dark range)
  breaks[1] <- -Inf
  breaks[NROW(breaks)] <- Inf
  dark_ranges <- breaks[dark_idx]
  
  # check whether matrix values are in dark ranges
  ls_range_condition <- lapply(seq_len(NROW(dark_ranges) / 2), function(i) {
    mat_min <- matrix(dark_ranges[2 * i - 1], nrow = NROW(mat_with_metric), ncol = NCOL(mat_with_metric))
    mat_max <- matrix(dark_ranges[2 * i], nrow = NROW(mat_with_metric), ncol = NCOL(mat_with_metric))
    mat_with_metric > mat_min & mat_with_metric <= mat_max
  })
  
  # prepare the final matrix with color names
  mat_number_color <- Reduce(pmax, ls_range_condition)
  mat_number_color[is.na(mat_number_color)] <- 0 # light_color_font for NA field (assumption: grey)
  mat_number_color[] <-
    vapply(mat_number_color, function(x) {
      if (x) {
        dark_color_font
      } else {
        light_color_font
      }
    }, character(1))
  
  return(mat_number_color)
}


#' Get fontsize for rownames or colnames in pheatmap::pheatmap 
#'
#' @param matrix numeric matrix with metric values.
#' @param dimension character value, either "row" or "col", indicating whether to compute fontsize for rows or columns.
#' @param threshold_count integer value of the number of rows/columns for which the font size remains standard.
#' 
#' @return numeric value of font size.
#' 
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#' 
#' @keywords internal
.get_pheatmap_fontsize <- function(matrix,
                                   dimension = c("row", "col"),
                                   threshold_count = 40L) {
  dimension <- match.arg(dimension)
  
  checkmate::assert_matrix(matrix)
  checkmate::assert_int(threshold_count, lower = 1)
  
  count <- if (dimension == "row") {
    NROW(matrix)
  } else {
    NCOL(matrix)
  }
  
  if (count > threshold_count) {
    0.6 * 8
  } else {
    8
  }
}


#' Trim labels to the required number of characters
#'
#' @param lbls_vec character vector with labels to be trimmed
#' @param max_lbl_length numeric value for the maximum number of characters in the label;
#'  if set to Inf, no trimming will be performed
#'
#' @returns character vectors with trimmed labels
#' @examples
#' \dontrun{
#' ls_lbls <- c(
#'   "short_lbl", "short_lbl", "veryveryverylong", 
#'   "long_duplicates|lbl_1AB", "long_duplicates|lbl_1AB", "long_duplicates|lbl_123"
#' )
#' 
#' .trim_labels(lbls_vec = ls_lbls)
#' .trim_labels(lbls_vec = ls_lbls, max_lbl_length = 15)
#' }
#' 
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#' 
#' @keywords internal
.trim_labels <- function(lbls_vec,
                         max_lbl_length = Inf) {
  
  checkmate::assert_character(lbls_vec, any.missing = FALSE)
  checkmate::assert_number(max_lbl_length, lower = 5)
  
  too_long_lbl <- which(nchar(lbls_vec) > max_lbl_length)
  
  if (NROW(too_long_lbl) == 0) {
    lbls_vec
  } else {
    trimed_vec_lbl <- lbls_vec
    names(trimed_vec_lbl) <- lbls_vec
    trimed_vec_lbl[too_long_lbl] <- 
      paste0(substr(trimed_vec_lbl[too_long_lbl], 1, max_lbl_length - 3), "...")
    # duplicates in the trimmed list
    if (any(duplicated(trimed_vec_lbl))) {
      trim_dup <- unique(trimed_vec_lbl[duplicated(trimed_vec_lbl)])
      
      for (trim_d in trim_dup) {
        i_dup_0 <- which(trimed_vec_lbl == trim_d)
        i_dup <- i_dup_0[sort(unique(names(i_dup_0)))]
        if (NROW(i_dup) == 1) next
        len_dup <- vapply(lbls_vec[i_dup], function(i) nchar(i), numeric(1))
        
        # find min number of character to distinguish strings
        min_distinguishing <-
          max(vapply(seq_len(NROW(i_dup) - 1), function(i) {
            str_1 <- 
              utils::head(c(strsplit(lbls_vec[i_dup[i]], split = "")[[1]], rep(NA, max(len_dup))), max(len_dup)) 
            str_2 <- 
              utils::head(c(strsplit(lbls_vec[i_dup[i + 1]], split = "")[[1]], rep(NA, max(len_dup))), max(len_dup)) 
            
            if (any(str_1 != str_2) && !is.na(any(str_1 != str_2))) {
              min(which(str_1 != str_2), na.rm = TRUE)
            } else {
              min(len_dup)
            }
          }, numeric(1)))
        
        # update duplicates in trimmed list
        trimed_vec_lbl[i_dup_0] <-
          vapply(i_dup_0, function(i) {
            if (nchar(lbls_vec[i]) <= min_distinguishing) {
              lbls_vec[i]
            } else {
              paste0(substr(lbls_vec[i], 1, min_distinguishing), "...")
            }
          }, character(1))
      }
    }
    # final trimmed
    trimed_vec_lbl
  }
}
