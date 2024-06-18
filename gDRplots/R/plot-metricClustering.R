#' Build labels for plotly tooltips for Metric Clustering heatmap
#'
#' @param data a data table in which the labels are to be constructed
#' 
#' @return A character vector with labels
#' 
#' @keywords utils_label
build_label_clustering <- function(data) {
  checkmate::assert_data_table(data)
  
  if (nrow(data) != nrow(unique(data))) {
    stop("non-unique rows identified in input")
  }
  
  pidfs <- gDRutils::get_prettified_identifiers(simplify = TRUE)
  
  # Handle co-treatment data.
  data <- coerce_cotreatment_data(data)
  exp_cotrt_ids <- c(pidfs[["drug_name2"]], pidfs[["concentration2"]])
  COTRT_DRUG_IDS <- intersect(exp_cotrt_ids, names(data))
  if (length(COTRT_DRUG_IDS) == length(exp_cotrt_ids)) {
    label <- sprintf("%s\n(%s at %.4g &mu;M)",
                     data[[pidfs[["drug_name"]]]], data[[pidfs[["drug_name2"]]]], data[[pidfs[["concentration2"]]]])
  } else {
    label <- data[[pidfs[["drug_name"]]]]
  }
  reformat_untreated_cases(label)
}


#' Prepare matrix for clustering and heat map construction.
#'
#' The function calls \code{dcast} to reshape the data.table.
#' \code{variable} is averaged at this stage if there are replicates.
#'
#' @param data a data.table containing data and metadata of metrics, conditions, and treatments
#' @param variable character string; name of variable to cluster on
#'
#' @return A named numeric matrix with drugs in rows and cell lines in columns.
#'
#' @keywords plugin_plot
#'
#' @seealso \code{plotly_metric_clustering}, \code{MetricClustering}
#'
#' @export
#'
prepareDataMH <- function(data, variable) {
  checkmate::assert_data_table(data)
  checkmate::assert_string(variable)
  
  pidfs <- gDRutils::get_prettified_identifiers(simplify = TRUE)
  
  DRUG_ID <- pidfs[["drug_name"]]
  CELL_ID <- pidfs[["cellline_name"]]
  
  CELL_METADATA <- c(pidfs[["cellline_tissue"]])
  DRUG_METADATA <- c(pidfs[["drug_moa"]])
  
  COTRT_DRUG_IDS <- intersect(c(pidfs[["drug_name2"]], pidfs[["concentration2"]]), names(data))
  COTRT_DRUG_METADATA <- intersect(c(pidfs[["drug_moa2"]]), names(data))
  
  keep <- c(CELL_ID, DRUG_ID, COTRT_DRUG_IDS, variable)
  
  # Transform data.
  form <- create_formula(CELL_ID, DRUG_ID, COTRT_DRUG_IDS)
  wide <-
    data.table::dcast(data = data[, .SD, .SDcols = keep], 
                      formula = form, 
                      fun.aggregate = mean, 
                      value.var = variable)
  
  # Drug annotations.
  row_ids <- c(DRUG_ID, COTRT_DRUG_IDS)
  # TODO: This can be removed once we operate on SE objects.
  drug_annotation <- unique(
    data[, .SD, .SDcols = c(DRUG_ID, DRUG_METADATA, COTRT_DRUG_IDS, COTRT_DRUG_METADATA)])
  # TODO: This can be called outside of a specific module in the management of the data.
  drug_annotation <- convert_factor_to_character(drug_annotation)
  identifiers <- wide[, .SD, .SDcols = row_ids, drop = FALSE]
  row_annotation <- map_annotations(identifiers, drug_annotation, row_ids)
  
  mat <- as.matrix(wide[, .SD, .SDcols = !row_ids]) # Enable coercion to matrix 
  rownames(mat) <- build_label_clustering(identifiers)
  
  if (variable %in% get_metrics_to_transform()) {
    mat <- log10(mat)
  }
  
  # Cell annotations.
  # TODO: This can be removed once we operate on SE objects.
  col_selected <- names(data)[names(data) %in% c(CELL_ID, CELL_METADATA)]
  cl_annotation <- unique(data[, .SD, .SDcols = col_selected])
  # TODO: This can be called outside of a specific module in the management of the data.
  cl_annotation <- .tidy_cell_metadata(cl_annotation)
  cl_df <- data.table::data.table("V1" = colnames(mat), check.names = FALSE)
  colnames(cl_df) <- CELL_ID
  col_annotation <- map_annotations(cl_df, cl_annotation, CELL_ID)
  
  annotations <- list(row_annotation = row_annotation, col_annotation = col_annotation)
  
  list(data_matrix = mat, annotations = annotations)
}

#' Draw interactive Metric Clustering heatmap
#'
#' Plot metric heat map with clustering
#'
#' Generate a heat map of a selected metric, clustering cell lines and drugs.
#' Distances are computed with \code{\link[gDRplots]{computeDistances}}.
#'
#' @param data a numeric matrix with dimnames
#' @param variable character string, naming metric to analyze
#' @param annotations list of of length 2,
#'                    containing row and column annotations as \code{data.table}s
#' @param transpose logical; whether or not to transpose the matrix
#' @param colors named list with colors for heatmap itself, rows and columns clusterings;
#'.   names: \code{heatmap}, \code{row} and \code{col}
#' 
#' @param dendrogram character string indicating whether to compute 'none',
#'    'row', 'column' or 'both' dendrograms. Defaults to 'both'.
#'
#' @return An interactive plot created by \code{iheatmapr}.
#'
#' @keywords plugin_plot
#'
#' @seealso \code{prepareDataMH} \code{MetricClustering}, \code{\link[gDRplots]{computeDistances}}
#'
#' @export
#'
plotly_metric_clustering <- function(data, 
                                     variable, 
                                     annotations, 
                                     transpose = FALSE, 
                                     colors, 
                                     dendrogram = "both") {
  
  checkmate::assert_matrix(data)
  checkmate::assert_numeric(data, finite = TRUE)
  checkmate::assert_string(variable)
  checkmate::assert_list(annotations)
  checkmate::assert_true(length(annotations) == 2)
  checkmate::assert_data_table(annotations$row_annotation)
  checkmate::assert_data_table(annotations$col_annotation)
  checkmate::assert_list(colors)
  checkmate::assert_names(names(colors), must.include = c("heatmap", "row", "col"))
  if (is.logical(dendrogram)) {
    dendrogram <- if (dendrogram)
      "both"
    else
      "none"
  }
  match.arg(dendrogram, c("both", "row", "column", "none"))
  
  pidfs <- gDRutils::get_prettified_identifiers(simplify = TRUE)
  row_name <- pidfs[["drug_name"]]
  col_name <- pidfs[["cellline_name"]]
  
  if (transpose) {
    data <- t(data)
    annotations[c("row_annotation", "col_annotation")] <- annotations[c("col_annotation", "row_annotation")]
    colors[c("row", "col")] <- colors[c("col", "row")]
    tmp <- row_name
    row_name <- col_name
    col_name <- tmp
  }
  
  main_label <- if (variable %in% get_metrics_to_transform()) {
    sub("50", "<sub>50</sub>", sprintf("log10(%s) [&mu;M]", variable))
  } else {
    variable
  }
  
  # create plot title
  plot_title <- sprintf("Clustering analysis of %s", main_label)
  
  # figure out layout
  # include names of annotations in computing margins
  # the "* 1.3 " is for "Secondary Drug Concentration" which still didn't fit
  longest_row_name <- max(nchar(c(rownames(data), names(annotations$row_annotation))))
  longest_col_name <- max(nchar(c(colnames(data), names(annotations$col_annotation)))) * 1.3
  font_size <- 225 / max(nrow(data), ncol(data))
  font_size <- if (font_size > 15) {
    15
  } else if (font_size < 9) {
    9
  } else {
    font_size
  }
  margin_left <- longest_row_name * font_size / 2
  margin_bottom <- longest_col_name * font_size / 2
  # adjust alignment of row/column labels
  row_labels <- adjust_label(rownames(data))
  col_labels <- adjust_label(colnames(data))
  
  # get heatmap range for specific variables, for others return min and max value from data
  heatmap_range_list <- get_visualization_range()
  heatmap_range <- if (!is.null(heatmap_range_list[[variable]])) {
    heatmap_range_list[[variable]]
  } else {
    c(min(data, na.rm = TRUE), max(data, na.rm = TRUE))
  }
  
  data_long <- data.table::as.data.table(as.table(data))
  data.table::setnames(data_long, c("V1", "V2", "N"), c("Var1", "Var2", "Freq"), skip_absent = TRUE)
  
  col_duplicates <- stats::aggregate(Freq ~ Var1, data_long, FUN = calc_duplicate_freq)
  row_duplicates <- stats::aggregate(Freq ~ Var2, data_long, FUN = calc_duplicate_freq)
  
  heatmap_hover_text <- matrix(sprintf(
    "%s: %s\n%s: %s\n%s: %.2f",
    row_name,
    data_long$Var1,
    col_name,
    data_long$Var2,
    main_label,
    data_long$Freq
  ),
  nrow = nrow(data),
  ncol = ncol(data))
  
  show_dendrogram <- rep(all(dim(data) >= 2), 2)
  data_index <- ifelse(all(show_dendrogram), 7, 3)
  
  # build plot
  heatmap_plot <- heatmaply::heatmaply(
    x = data,
    plot_method = "plotly",
    main = plot_title,
    limit = heatmap_range,
    Rowv = nrow(data) >= 2 && !any(row_duplicates$Freq > 0.95),
    Colv = ncol(data) >= 2 && !any(col_duplicates$Freq > 0.95),
    dendrogram = dendrogram,
    show_dendrogram = show_dendrogram,
    colors = create_color_palette(colors$heatmap, heatmap_range),
    key.title = main_label,
    row_side_colors = annotations$row_annotation,
    row_side_palette = colors$row,
    fontsize_row = font_size,
    col_side_colors = annotations$col_annotation,
    col_side_palette = colors$col,
    column_text_angle = 90,
    fontsize_col = font_size,
    distfun = gDRplots::computeDistances,
    dend_hoverinfo = FALSE,
    margins = c(margin_bottom, margin_left, NA, 0),
    seriate = "mean", # matrix sorting
    custom_hovertext = heatmap_hover_text,
    colorbar_xpos = 1.02,
    colorbar_len = 0.2,
    midpoint = 0
  )
  
  # Turn off hover for annotations
  heatmap_plot$x[[1]] <- lapply(heatmap_plot$x[[1]], function(x) {
    if (any(dim(x$z) != dim(data))) {
      x$hoverinfo <- "skip"
      x
    } else {
      x
    }
  })
  
  heatmap_plot$x$source <- "heatmap_plot"
  heatmap_plot$x$data_index <- data_index
  gDR_plotly_config(force_heatmaply_limits(heatmap_plot, heatmap_range))
}


#' map identifiers to their corresponding annotation metadata
#'
#' @param identifiers character vector of identifiers to map to annotations, i.e. cell line or drug identifiers
#' @param annotation data.table containing column \code{from} and other metadata on the identifiers
#' @param from character vector column names in \code{annotation} corresponding to columns in \code{identifiers}.
#' Order of columns in \code{from} should correspond to order in \code{identifiers}.
#'
#' @return character vector the same length as \code{identifiers} containing the corresponding annotations
#' for the \code{identifiers}
#' @importFrom S4Vectors match DataFrame
#' @keywords internal
map_annotations <- function(identifiers, annotation, from) {
  checkmate::assert_data_table(identifiers)
  checkmate::assert_data_table(annotation)
  
  if (length(identifiers) != length(unique(identifiers))) {
    stop("'identifiers' contains non-unique entries")
  }
  
  if (nrow(annotation) != nrow(unique(annotation))) {
    stop("'annotation' contains non-unique entries")
  }
  
  if (!all(names(identifiers) %in% names(annotation))) {
    stop(sprintf("missing required columns in 'annotation': '%s'", paste0(names(identifiers), collapse = ", ")))
  }
  
  # No annotations.
  if (length(annotation) <= length(identifiers)) {
    return(identifiers)
  }
  
  keep <- !colnames(annotation) %in% from
  out <- annotation[
    S4Vectors::match(S4Vectors::DataFrame(identifiers),
                     S4Vectors::DataFrame(annotation[, .SD, .SDcols = from, drop = FALSE])), 
    .SD, .SDcols = keep, drop = FALSE]
  out
}


# TODO: This can be done when the assay metrics are retrieved.
#' @keywords internal
.tidy_cell_metadata <- function(cl_annotation) {
  pt <- gDRutils::get_prettified_identifiers("cellline_tissue", simplify = TRUE)
  
  if (pt %in% colnames(cl_annotation)) {
    cl_annotation[[pt]] <- tolower(cl_annotation[[pt]])
  }
  cl_annotation
}

#' Create formula
#'
#' @param cell_id String of cell identifier column.
#' @param drug_id String of drug identifier column.
#' @param cotrt_drug_ids Character vector of any cotreatment drug experimental design identifier columns.
#' @return formula accounting for any cotreatment columns
#' @keywords internal
create_formula <- function(cell_id, drug_id, cotrt_drug_ids) {
  form <- sprintf("`%s` ~ `%s`", drug_id, cell_id)
  if (length(cotrt_drug_ids) > 0) {
    cotrt_drug_ids <-
      as.character(vapply(cotrt_drug_ids, function(x) {
        sprintf("`%s`", x)
      }, character(1)))
    form <- paste(paste(cotrt_drug_ids, collapse = " + "), form, sep = " + ")
  }
  stats::as.formula(form)
}

#' @keywords internal
calc_duplicate_freq <- function(x) {
  sum(duplicated(x)) / length(x)
}

#' Convert data.table factor columns to character
#'
#' @param tbl data.table
#' 
#' @examples
#' dt <- data.table::data.table(
#'   col1 = c("1","2","3"), 
#'   col2 = c("a", "b", "a"), 
#'   stringsAsFactors = TRUE)
#' dt
#' convert_factor_to_character(dt)
#' 
#' @return data.table with converted columns
#' @keywords internal
#' 
#' @export
convert_factor_to_character <- function(tbl) {
  checkmate::assert_data_table(tbl)
  
  factor_cols <- names(tbl)[vapply(tbl, is.factor, logical(1))]
  tbl[, (factor_cols) := lapply(.SD, as.character), .SDcols = factor_cols]
  
  tbl
}
