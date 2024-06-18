#' Draw interactive Drug Combo heatmap
#'
#' Build experiment overview heat map.
#'
#' Produces a clustered heat map of cell line vs. drug pair.
#' The app registers clicks on panes and opens a pop-up with
#' detailed exploration of the drug interaction.
#'
#' @param data a list of data.table(s) with combo data
#' @param metric_combo character string specifying which combo score metric is used;
#'                     one of \code{names(gDRutils::get_combo_score_field_names())}
#' @param metric_growth character string specifying which growth metric is used;
#'                      one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#'
#' @return plotly object - heatmap
#'
#' @examples
#' \dontrun{
#' mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
#' combo_data_l <- gDRutils::convert_combo_data_to_dt(mae[[gDRutils::get_supported_experiments("combo")]])
#' plotlyDCMother(data = combo_data_l, 
#'                metric_combo = "hsa_score", 
#'                metric_growth = "RV")
#'
#' plotlyDCMother(data = combo_data_l, 
#'                metric_combo = "CIScore_50", 
#'                metric_growth = "GR")
#' }
#' 
#' @keywords plugin_plot
#' @export
#'
plotlyDCMother <- function(data, 
                           metric_combo, 
                           metric_growth) {
  
  checkmate::assert_list(data)
  checkmate::assert_data_table(data[[1]])
  checkmate::assert_string(metric_combo)
  checkmate::assert_string(metric_growth)
  checkmate::assert_choice(metric_growth, c("GR", "RV"))
  
  cscore_anames <- gDRutils::get_combo_score_field_names()
  checkmate::assert_choice(metric_combo, names(cscore_anames))
  metric_combo <- match.arg(metric_combo, names(cscore_anames), FALSE)
  assay_name <- gDRutils::convert_combo_field_to_assay(metric_combo)
  
  data <- get_combo_score_matrix(data = data, c_assay = metric_combo, normalization_type = metric_growth)
  
  if (all(is.na(data))) {
    return(plotly::layout(plotly::plotly_empty(), title = "no data to display"))
  }
  
  # transform data, if necessary
  # TODO: assure data has not be transformed already
  ci_scores <- as.character(cscore_anames[c("CIScore_50", "CIScore_80")])
  if (is.element(metric_combo, ci_scores)) {
    data <- log2(data)
  }
  # prepare labels
  label_combo <- cscore_anames[[metric_combo]]
  label_growth <- switch(metric_growth,
                         "GR" = "GR value",
                         "RV" = "Relative Viability")
  main_label <- label_combo
  plot_title <- sprintf("%s with %s", main_label, label_growth)
  cs <- gDRutils::get_combo_col_settings(metric_growth, assay_name)
  
  # obtain layout properties
  longest_row_name <- max(nchar(rownames(data)))
  longest_col_name <- max(nchar(colnames(data)))
  font_size <- floor(225 / max(nrow(data), ncol(data)))
  font_size <- if (font_size > 15) {
    15
  } else if (font_size < 9) {
    9
  } else {
    font_size
  }
  margin_left <- ceiling(longest_row_name * font_size / 2)
  margin_bottom <- ceiling(longest_col_name * font_size / 2)
  
  #TODO:
  # 1. better color handling
  # https://github.com/talgalili/heatmaply/issues/156
  # https://github.com/talgalili/heatmaply/issues/1
  # 2. colorbar position
  # see xanchor and yanchor below
  
  data_long <- data.table::as.data.table(as.table(data))
  heatmap_hover_text <- matrix(
    sprintf(
      "%s: %s\n%s: %s\n%s: %.2f",
      "Drug Pair",
      data_long$V1,
      gDRutils::get_prettified_identifiers("cellline_name"),
      data_long$V2,
      main_label,
      data_long$N
    ),
    nrow = nrow(data),
    ncol = ncol(data)
  )
  
  show_dendrogram <- rep(all(dim(data) >= 2), 2)
  data_index <- ifelse(all(show_dendrogram), 3, 1)
  
  # build plot
  heatmap_combo <- heatmaply::heatmaply(
    x = data,
    plot_method = "plotly",
    main = plot_title,
    limit = cs$limits,
    Rowv = nrow(data) >= 2,
    Colv = ncol(data) >= 2,
    show_dendrogram = show_dendrogram,
    colors = create_color_palette(cs$colors, cs$limits),
    key.title = main_label,
    fontsize_row = font_size,
    column_text_angle = 90,
    fontsize_col = font_size,
    distfun = gDRplots::computeDistances,
    dend_hoverinfo = FALSE,
    margins = c(margin_bottom, margin_left, NA, 0),
    seriate = "mean", # matrix sorting
    custom_hovertext = heatmap_hover_text,
    colorbar_xpos = 1.02,
    colorbar_len = 0.2
  )
  heatmap_combo$x$source <- "heatmap_combo"
  heatmap_combo$x$data_index <- data_index
  gDR_plotly_config(force_heatmaply_limits(heatmap_combo, cs$limits))
  
}

#' Draw interactive Drug Combo details heatmap with isobolograms
#'
#' Functions for detailed visualization of combo drug treatments, i.e. 2-way concentration arrays.
#'
#' @param data a list of data.table(s) with combo data
#' @param drug1_name string a name of the first drug to be selected
#' @param drug2_name string a name of the second drug to be selected
#' @param cell_line string a name of cell line to be selected
#' @param c_assay character string specifying the combo base assay to be selected,
#'        any from \code{names(gDRutils::get_combo_excess_field_names())} will do
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param iso_levels charvec isobologram levels to be selected
#' @param with_iso_points logical flag whether draw the point for IC50/GR50 isobologram
#' @param height integer; heatmap height
#' @param width integer; heatmap width
#' 
#' @return plotly object - heatmap for 2-way concentration arrays with isobolograms
#'         
#' @examples
#' \dontrun{
#' mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
#' combo_data_l <- gDRutils::convert_combo_data_to_dt(mae[[gDRutils::get_supported_experiments("combo")]])
#'
#'plotlyDC1(
#'  data = combo_data_l,
#'  drug1_name = "drug_005",
#'  drug2_name = "drug_021",
#'  cell_line = "cellline_HB",
#'  normalization_type = "GR",
#'  c_assay = "hsa_excess"
#')
#'
#' plotlyDC1(
#'   data = combo_data_l,
#'   drug1_name = "drug_005",
#'   drug2_name = "drug_021",
#'   cell_line = "cellline_HB",
#'   normalization_type = "RV",
#'   c_assay = "smooth"
#' )
#'
#' plotlyDC1(
#'   data = combo_data_l,
#'   drug1_name = "drug_005",
#'   drug2_name = "drug_021",
#'   cell_line = "cellline_HB",
#'   normalization_type = "GR",
#'   c_assay = "bliss_excess",
#'   iso_levels = c("0.25", "0.5", "0.75")
#' )
#' }
#'
#' @keywords plugin_plot
#' @export
#'
plotlyDC1 <- function(data,
                      drug1_name,
                      drug2_name,
                      cell_line,
                      c_assay = names(gDRutils::get_combo_excess_field_names())[1],
                      normalization_type = "RV",
                      iso_levels =  c("0.5", "0.75"),
                      with_iso_points = FALSE,
                      height = 325L,
                      width = 425L) {
  # argument checks
  checkmate::assert_list(data)
  checkmate::assert_data_table(data[[1]])
  if (is.na(drug1_name)) { # Return NULL for cases with no drugs (reactive conflict for sa and combo matrix heatmaps)
    return(NULL)
  }
  checkmate::assert_string(drug1_name)
  checkmate::assert_string(drug2_name)
  checkmate::assert_string(cell_line)
  
  checkmate::assert_string(normalization_type)
  checkmate::assert_choice(normalization_type, c("GR", "RV"))
  checkmate::assert_string(c_assay)
  cbase_anames <- gDRutils::get_combo_excess_field_names()
  checkmate::assert_choice(c_assay, names(cbase_anames))
  checkmate::assert_numeric(height)
  checkmate::assert_number(width)
  
  checkmate::assert_character(iso_levels)
  # TODO: get values from isobologram assay (iso_level) column
  gDRutils::assert_choices(iso_levels, as.character(seq(0.01, 0.99, by = 0.01)))
  
  iso_colors <- gDRutils::get_iso_colors(normalization_type = normalization_type)
  
  # get prettified identfiers
  pidfs <- gDRutils::get_prettified_identifiers()
  pi_drug_name <- pidfs[["drug_name"]]
  pi_drug_name2 <- pidfs[["drug_name2"]]
  pi_drug <- pidfs[["drug"]]
  pi_drug2 <- pidfs[["drug2"]]
  pi_cell_name <- pidfs[["cellline_name"]]
  pi_conc_name <- pidfs[["concentration"]]
  pi_conc2_name <- pidfs[["concentration2"]]
  pi_duration_name <- pidfs[["duration"]]
  sdata <-
    get_combo_base_data(
      data,
      drug1_name,
      drug2_name,
      cell_line,
      c_assay,
      normalization_type = normalization_type
    )
  
  if (is.null(sdata)) {
    return(NULL)
  }
  
  drug1_axis <- sdata$drug1_axis
  drug2_axis <- sdata$drug2_axis
  
  
  # transform data for isobolograms
  data_iso <-
    get_isobologram_data(
      data,
      drug1_name,
      drug2_name,
      cell_line,
      c_assay = names(gDRutils::get_combo_assay_names(group = "combo_iso"))[1],
      normalization_type = normalization_type,
      drug1_axis = drug1_axis,
      drug2_axis = drug2_axis,
      iso_levels = iso_levels
    )
  
  condition <- sdata$condition
  checkmate::assert_data_table(drug1_axis)
  checkmate::assert_data_table(drug2_axis)
  checkmate::assert_character(normalization_type)
  checkmate::assert_character(condition)
  checkmate::assert_named(condition)
  
  # prepare data
  data <- sdata$matrix
  
  data_long <- data.table::melt(data.table::as.data.table(data, keep.rownames = TRUE),
                                measure.vars = colnames(data), variable.factor = FALSE)
  data.table::setnames(data_long, c("conc_1", "conc_2", "value"))
  data_long[["value"]] <- pmin(1.1, data_long[["value"]])
  data_long[, c("conc_1", "conc_2")] <-
    lapply(data_long[, c("conc_1", "conc_2")], function(x) {
      round(as.numeric(x), digits = 5)
    })
  data_long <- data_long[drug1_axis, on = "conc_1", nomatch = NULL]
  data_long <- data_long[drug2_axis, on = "conc_2", nomatch = NULL]
  # prepare plot title
  matrix_pretty <- gDRutils::get_combo_excess_field_names()[[c_assay]]
  title_plot <- sprintf("%s (%s for %s, T = %sh)", condition[[pi_cell_name]], matrix_pretty,
                        normalization_type, condition[[pi_duration_name]])
  # axis properties
  tickvals_x <- drug2_axis$pos_x
  tickvals_y <- drug1_axis$pos_y
  ticktext_x <- sprintf("%.2g", drug2_axis$conc_2)
  ticktext_y <- sprintf("%.2g", drug1_axis$conc_1)
  range_x <- range(drug2_axis$pos_x)
  range_y <- range(drug1_axis$pos_y)
  title_x <- sprintf("%s [&mu;M]", condition[[pi_drug_name2]])
  title_y <- sprintf("%s [&mu;M]", condition[[pi_drug_name]])
  hover_label <- list(bgcolor = "white")
  # line properties
  slope <- 1
  line_limits_x <- range_x
  line_limits_y <- range_y
  
  # color bar properties
  assay_name <- gDRutils::convert_combo_field_to_assay(c_assay)
  cs <- gDRutils::get_combo_col_settings(normalization_type, assay_name)
  title_colorbar <-
    gDRutils::get_assay_names(type = assay_name, prettify = TRUE)
  if (title_colorbar == "MX full") {
    title_colorbar <- normalization_type
  }
  
  # fixing column names
  isobologram_columns <- gDRutils::get_isobologram_columns(prettify = FALSE)
  names(data_long)[which(names(data_long) %in% isobologram_columns)] <-
    unlist(lapply(names(data_long), function(x) {
      if (x %in% isobologram_columns) {
        gDRutils::get_isobologram_columns(x)
      }
    }))
  pos_x <- gDRutils::get_isobologram_columns("pos_x")
  pos_y <- gDRutils::get_isobologram_columns("pos_y")
  pos_x_ref <- gDRutils::get_isobologram_columns("pos_x_ref")
  pos_y_ref <- gDRutils::get_isobologram_columns("pos_y_ref")
  iso_level <- gDRutils::get_isobologram_columns("iso_level")
  # build tooltips
  data_long[["label"]] <- build_label(data_long, "combo1-heatmap")
  # build plot
  plot_base <- plotly::plot_ly(type = "heatmap",
                               x = data_long[[pos_x]], y = data_long[[pos_y]],
                               z = data_long[["value"]],
                               colors = cs$colors, colorbar = list(title = title_colorbar,
                                                                   tickvals = cs$breaks,
                                                                   ticktext = sprintf("%.2f", cs$breaks)),
                               zmin = cs$limits[1], zmax = cs$limits[2],
                               text = data_long[["label"]], hoverinfo = "text",
                               hoverlabel = hover_label,
                               height = height, width = width)
  
  
  
  # draw isobologram data
  ## add points of the isobologram at GR/IC50
  if (with_iso_points && any("0.5" %in% data_iso[[iso_level]])) {
    data <- data_iso[data_iso[[iso_level]] == "0.5", ]
    data[["label"]] <- build_label(data, "combo1-points")
    plot_base <- plotly::add_markers(plot_base, inherit = FALSE,
                                     x = data[[pos_x]], y = data[[pos_y]],
                                     text = data[["label"]], color = I("black"),
                                     symbol = data[["fit_type"]], symbols = c(1, 5),
                                     showlegend = FALSE)
  }
  
  ## add isobolograms as lines
  if (NROW(data_iso) > 0) {
    metric <- switch(normalization_type, "GR" = "GR", "RV" = "IC")
    data_iso[["name"]] <- sprintf("%s%i", metric, 100 - 100 * as.numeric(data_iso[[iso_level]]))
    data_iso[["label_ref"]] <- build_label(data_iso, "combo1-lines_ref")
    data_iso[["label"]] <- build_label(data_iso, "combo1-lines")
    
    # add lines
    for (i in seq_along(unique(data_iso[[iso_level]]))) {
      # subset data for single iso value
      data_subset <- data_iso[data_iso[[iso_level]] == unique(data_iso[[iso_level]])[i], ]
      if (all(is.na(data_subset[[pos_x]]))) {
        next
      }
      # workaround #1 for labels
      # see https://github.com/plotly/plotly.R/issues/2069 for more details
      # use labels from the first layer (heatmap) for the ticks defined in the subsequent layers
      # (i.e. ticks for isobolograms and reference isobolograms)
      ml <- map_coords(data_long, data_subset)
      data_subset[ml$data_idx, "label"] <- data_long[ml$ref_idx, "label"]
      
      ml_ref <- map_coords(data_long, data_subset, d_cols = c(pos_x_ref, pos_y_ref))
      data_subset[ml_ref$data_idx, "label_ref"] <- data_long[ml_ref$ref_idx, "label"]
      
      # set color
      iso_color <- as.character(iso_colors[data_subset[[1, iso_level]]])
      # add reference isobol
      plot_base <- plotly::add_paths(plot_base, inherit = FALSE,
                                     x = data_subset[[pos_x_ref]], y = data_subset[[pos_y_ref]],
                                     split = data_subset[[iso_level]],
                                     line = list(color = iso_color, width = 1, dash = "10px"),
                                     text = data_subset[["label_ref"]],
                                     hoverinfo = "text",
                                     showlegend = FALSE,
                                     hoverlabel = hover_label)
      # add actual isobol
      plot_base <- plotly::add_paths(plot_base, inherit = FALSE,
                                     x = data_subset[[pos_x]], y = data_subset[[pos_y]],
                                     split = data_subset[[iso_level]], name = data_subset[["name"]],
                                     line = list(color = iso_color, width = 2),
                                     text = data_subset[["label"]],
                                     hoverinfo = "text",
                                     showlegend = TRUE,
                                     hoverlabel = hover_label)
    }
  }
  
  # apply layout
  plot_laidout <- plotly::layout(plot_base,
                                 title = title_plot,
                                 xaxis = list(title = title_x, range = range_x,
                                              tickvals = tickvals_x, ticktext = ticktext_x,
                                              zeroline = FALSE, showgrid = FALSE, showline = TRUE, mirror = TRUE),
                                 yaxis = list(title = title_y, range = range_y,
                                              tickvals = tickvals_y, ticktext = ticktext_y,
                                              zeroline = FALSE, showgrid = FALSE, showline = TRUE, mirror = TRUE),
                                 legend = list(title = "Isobol")
  )
  # modify config
  plot_final <- gDR_plotly_config(plot_laidout,
                                  editable = TRUE,
                                  showAxisRangeEntryBoxes = FALSE)
  return(plot_final)
  
}

#' get base combo data of given filters
#'
#' @examples
#' \dontrun{
#' mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
#' combo_data_l <- gDRutils::convert_combo_data_to_dt(mae[[gDRutils::get_supported_experiments("combo")]])
#'
#' get_combo_base_data(
#'  data = combo_data_l,
#'  drug1_name = "drug_004",
#'  drug2_name = "drug_021",
#'  cell_line = "cellline_GB",
#'  c_assay = "bliss_excess"
#' )
#' }
#' @inheritParams plotlyDC1
#'
#' @author Arkadiusz Gładki \email{arkadiusz.gladki@@contractors.roche.com}
#'
#' @keywords internal
#' @export
#'
get_combo_base_data <-
  function(data,
           drug1_name,
           drug2_name,
           cell_line,
           c_assay = names(gDRutils::get_combo_base_assay_names())[[1]],
           normalization_type = "RV") {
    
    s_ntype <- normalization_type
    
    pidfs <- gDRutils::get_prettified_identifiers()
    pi_drug_name <- pidfs[["drug_name"]]
    pi_drug_name2 <- pidfs[["drug_name2"]]
    pi_drug <- pidfs[["drug"]]
    pi_drug2 <- pidfs[["drug2"]]
    pi_cell_name <- pidfs[["cellline_name"]]
    pi_conc_name <- pidfs[["concentration"]]
    pi_conc2_name <- pidfs[["concentration2"]]
    pi_duration_name <- pidfs[["duration"]]
    
    ntype_name <- gDRutils::prettify_flat_metrics("normalization_type", human_readable = TRUE)
    
    smooth_mx <- "smooth"
    assay_name <- gDRutils::convert_combo_field_to_assay(c_assay)
    mx_assay_name <- gDRutils::convert_combo_field_to_assay(smooth_mx)
    p_c_assay <- gDRutils::prettify_flat_metrics(c_assay, human_readable = TRUE)
    
    # get axis list
    ## get name of the "SmoothMatrix" assay
    mx_dt <- data[[mx_assay_name]]
    
    mx_dt <-
      mx_dt[mx_dt[[pi_drug_name]] == drug1_name &
              mx_dt[[pi_drug_name2]] == drug2_name &
              mx_dt[[ntype_name]] %in% s_ntype &
              mx_dt[[pi_cell_name]] == cell_line, ]
    
    if (!checkmate::test_data_table(mx_dt, min.rows = 1L)) {
      return(NULL)
    }
    
    axes_l <- get_drug_axes(mx_dt, smooth_mx)
    
    if (is.null(axes_l)) {
      return(NULL)
    }
    
    # different data model for 'BlissExcess' vs 'SmoothMatrix'
    # standardized in gDRutils::convert_combo_data_to_dt<-
    sdt <- if (c_assay == "smooth") {
      mx_dt
    } else {
      dt <- data[[assay_name]]
      dt[dt[[pi_drug_name]] == drug1_name &
           dt[[pi_drug_name2]] == drug2_name &
           dt[[ntype_name]] %in% s_ntype &
           dt[[pi_cell_name]] == cell_line, ]
    }
    
    dt_ <- data.table::dcast(
      sdt,
      get(pi_conc_name) ~ get(pi_conc2_name),
      value.var = p_c_assay,
      fun.aggregate = function(x) mean(x, na.rm = TRUE)
    )
    
    mat <- as.matrix(dt_[, -"pi_conc_name"])
    rownames(mat) <- dt_$pi_conc_name
    
    ids <- intersect(c(pi_drug_name, pi_drug_name2, pi_drug, pi_drug2, pi_duration_name, pi_cell_name), names(sdt))
    cond_dt <- sdt[1, ids, with = FALSE]
    cond_nv <- structure(as.character(cond_dt), names = names(cond_dt))
    list(matrix = mat, drug1_axis = axes_l[[1]], drug2_axis = axes_l[[2]], condition = cond_nv)
  }

#' get drug axes
#'
#' define the position to nicely plot the combination matrix such that the marginal (single agent)
#'   are slightly offset from the combination matrix and that the tiles are properly spaced.
#'   This is also required to get the position to plot the isobolograms
#'
#' @param dt data.table
#' @param value_col string default value is `smooth`
#'
#' @author Arkadiusz Gładki \email{arkadiusz.gladki@@contractors.roche.com}
#'
#' @keywords internal
#' 
get_drug_axes <- function(dt, value_col = "smooth") {
  pidfs <- gDRutils::get_prettified_identifiers()
  
  pi_conc_name <- pidfs[["concentration"]]
  pi_conc2_name <- pidfs[["concentration2"]]
  pi_value_col <-
    gDRutils::prettify_flat_metrics(value_col, human_readable = TRUE)
  dt_response <- data.table::dcast(
    dt,
    get(pi_conc_name) ~ get(pi_conc2_name),
    value.var = pi_value_col,
    fun.aggregate = function(x) {
      mean(x, na.rm = TRUE)
    }
  )
  
  mx_response <- as.matrix(dt_response[, -"pi_conc_name"])
  rownames(mx_response) <- dt_response$pi_conc_name
  
  # fill in cases where conc == 0 and conc_2 == 0 which would not exist in the data, but would be 1
  mx_response[1, 1] <- 1
  #remove empty rows/columns
  mx_response <-
    mx_response[rowSums(!is.na(mx_response)) > 2, colSums(!is.na(mx_response)) > 2]
  
  if (is.null(rownames(mx_response))) {
    return(NULL)
  }
  
  # drug_1 is diluted along the rows and will be the y-axis of the matrix plots
  drug1_axis <-
    data.table::data.table(
      conc_1 = .round_conc(rownames(mx_response)),
      log10conc1 = 0,
      pos_y = 0,
      marks_y = 0
    )
  drug1_axis$log10conc1 <- log10(drug1_axis$conc)
  drug1_axis$pos_y <- drug1_axis$log10conc1
  drug1_axis$pos_y[1] <-
    2 * drug1_axis$pos_y[2] - drug1_axis$pos_y[3] - log10(1.5)
  drug1_axis$marks_y <- sprintf("%.2g", drug1_axis$conc_1)
  
  # drug_2 is diluted along the columns and will be the x-axis of the matrix plots
  drug2_axis <-
    data.table::data.table(
      conc_2 = .round_conc(colnames(mx_response)),
      log10conc2 = 0,
      pos_x = 0,
      marks_x = 0
    )
  drug2_axis$log10conc2 <- log10(drug2_axis$conc_2)
  drug2_axis$pos_x <- drug2_axis$log10conc2
  drug2_axis$pos_x[1] <-
    2 * drug2_axis$pos_x[2] - drug2_axis$pos_x[3] - log10(1.5)
  drug2_axis$marks_x <- sprintf("%.2g", drug2_axis$conc_2)
  list(drug1_axis, drug2_axis)
}

#' get score combo data for given filters
#' 
#' @examples
#' \dontrun{
#' mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' combo_data_l <- gDRutils::convert_combo_data_to_dt(se)
#' get_combo_score_matrix(combo_data_l, "bliss_score")
#' }
#'
#' @inheritParams plotlyDC1
#' @param c_assay character string specifying the combo score assay to be selected,
#'        any from \code{names(gDRutils::get_combo_score_field_names())} will do
#'
#' @return matrix
#'
#' @author Arkadiusz Gładki \email{arkadiusz.gladki@@contractors.roche.com}
#'
#' @keywords internal
#' @export
#'
get_combo_score_matrix <-
  function(data,
           c_assay = names(gDRutils::get_combo_score_field_names())[[1]],
           normalization_type = "RV") {
    
    p_c_assay <- gDRutils::prettify_flat_metrics(c_assay, human_readable = TRUE)
    checkmate::assert_choice(normalization_type, c("RV", "GR"))
    
    
    idfs <- gDRutils::get_env_identifiers()
    pidfs <- gDRutils::get_prettified_identifiers()
    
    # TODO: change to pidfs once columns in combo data are prettified
    drug_name <- pidfs[["drug_name"]]
    drug2_name <- pidfs[["drug_name2"]]
    cell_name <- pidfs[["cellline_name"]]
    
    
    assay_name <- gDRutils::convert_combo_field_to_assay(c_assay)
    checkmate::assert_list(data)
    dt <- data[[assay_name]]
    dt[, "drug_combo" := paste(get(drug_name), get(drug2_name), sep = " x ")]
    norm_col <- gDRutils::prettify_flat_metrics("normalization_type", human_readable = TRUE)
    data.table::setkeyv(dt, norm_col)
    
    out_dt <- data.table::dcast(dt[normalization_type],
                                drug_combo ~ get(cell_name),
                                value.var = p_c_assay,
                                fun.aggregate = mean)
    
    out_matrix <- as.matrix(out_dt[, -"drug_combo"])
    rownames(out_matrix) <- out_dt$drug_combo
    
    out_matrix
  }

#' get isobologram combo data for given filters
#'
#' @inheritParams plotlyDC1
#' @param c_assay character string specifying the combo iso assay to be selected,
#'        any from \code{gDRutils::get_combo_assay_names(group = "combo_iso")} will do
#' @param drug1_axis data.table with drug1 data from \code{get_combo_base_data}
#' @param drug2_axis data.table with drug2 data from \code{get_combo_base_data}
#'
#' @author Arkadiusz Gładki \email{arkadiusz.gladki@@contractors.roche.com}
#'
#' @keywords internal
#' @export
#'
get_isobologram_data <-
  function(data,
           drug1_name,
           drug2_name,
           cell_line,
           drug1_axis,
           drug2_axis,
           c_assay = names(gDRutils::get_combo_assay_names(group = "combo_iso"))[[1]],
           normalization_type = "RV",
           iso_levels = c("0.5", "0.75")) {
    
    
    checkmate::assert_choice(normalization_type, c("GR", "RV"))
    
    s_ntype <- normalization_type
    ntype_name <- gDRutils::prettify_flat_metrics("normalization_type", human_readable = TRUE)
    int_cols <- gDRutils::get_isobologram_columns()
    iso_level_name <- gDRutils::get_isobologram_columns("iso_level")
    
    ciso_anames <- gDRutils::get_combo_assay_names(group = "combo_iso")
    checkmate::assert_choice(c_assay, names(ciso_anames))
    
    idfs <- gDRutils::get_env_identifiers()
    pidfs <- gDRutils::get_prettified_identifiers()
    pi_drug_name <- pidfs[["drug_name"]]
    pi_drug_name2 <- pidfs[["drug_name2"]]
    pi_cell_name <- pidfs[["cellline_name"]]
    
    pos_x <- gDRutils::get_isobologram_columns("pos_x")
    pos_y <- gDRutils::get_isobologram_columns("pos_y")
    
    checkmate::assert_list(data)
    dt <- data[[ciso_anames[[c_assay]]]]
    
    # temporary fix: internal columns in isobologram
    cols_idx <- match(gDRutils::prettify_flat_metrics(int_cols, human_readable = TRUE), colnames(dt))
    if (all(is.na(cols_idx))) {
      warning("No isobolgram data available for this combination")
      return(NULL)
    }
    colnames(dt)[cols_idx] <- gsub(" ", "_", int_cols)
    dt <- stats::na.omit(dt, col = pos_x)
    
    
    # rotate 45 degree to calculate smooth curve:
    dt$x1 <- (dt[[pos_x]] - min(drug2_axis$pos_x) -
                (dt[[pos_y]] - min(drug1_axis$pos_y))) / sqrt(2) # conc_ratio
    dt$x2 <- (dt[[pos_x]] - min(drug2_axis$pos_x) +
                (dt[[pos_y]] - min(drug1_axis$pos_y))) / sqrt(2) # new response value
    x2_extra_offset <- 1 / 4 # offset helps with smoothing of the edges
    dt$x2_off <- dt$x2 + abs(dt$x1) * x2_extra_offset
    
    # TODO: iso_level added to prettify functions?
    data.table::setkeyv(dt, c(pi_drug_name,
                              pi_drug_name2, pi_cell_name, iso_level_name, ntype_name))
    
    dt[list(drug1_name, drug2_name, cell_line, iso_levels, s_ntype)]
  }

#' Map x,y coordinates from one given data.table to the x,y coords of reference data.table
#'
#' @param ref data.table with refrence data
#' @param data data.table with data of interest
#' @param r_cols two-element charvec with column names for x,y position in reference data.table
#' @param d_cols two-element charvec with column names for x,y position in data data.table
#' 
#' @return row indices of data.table with x,y coords mapped to x,y coords of data data.table
#' 
#' @keywords internal
#' 
map_coords <-
  function(ref,
           data,
           r_cols = c("Pos_x", "Pos_y"),
           d_cols = c("Pos_x", "Pos_y")) {
    
    checkmate::assert_data_table(ref, min.rows = 2)
    checkmate::assert_data_table(data, min.rows = 2)
    checkmate::assert_character(r_cols, min.len = 2, max.len = 2)
    checkmate::assert_character(d_cols, min.len = 2, max.len = 2)
    
    # map only coords from data within the x,y ranges of ref
    # omit the records from data not fullfilling these criteria
    d1_min_idx <-
      which(data[[d_cols[1]]] < min(ref[[r_cols[[1]]]]))
    d1_max_idx <-
      which(data[[d_cols[1]]] > max(ref[[r_cols[[1]]]]))
    d2_min_idx <-
      which(data[[d_cols[2]]] < min(ref[[r_cols[[2]]]]))
    d2_max_idx <-
      which(data[[d_cols[2]]] > max(ref[[r_cols[[2]]]]))
    data_idx <-
      setdiff(seq_len(nrow(data)), unique(c(
        d1_min_idx, d1_max_idx, d2_min_idx, d2_max_idx
      )))
    data <- data[data_idx, ]
    
    ref[["orig_order"]] <- seq_len(nrow(ref))
    data.table::setorderv(ref, r_cols)
    # get up limits
    ulx <- calc_up_limes(ref[[r_cols[1]]])
    ref <- merge(ref, ulx, by.x = r_cols[1], by.y = "ux")
    uly <- calc_up_limes(ref[[r_cols[2]]], cnames = c("uy", "uly"))
    ref <- merge(ref, uly, by.x = r_cols[2], by.y = "uy")
    # assure default sorting (by pos_x first and then by pos_y)
    data.table::setorderv(ref, r_cols)
    
    # calc indexes
    idx <- vapply(seq_len(nrow(data)), function(x) {
      xdif_v <- ref$ulx - data[[d_cols[1]]][x]
      xd_idx <- which(xdif_v >= 0)
      xd_idx2 <- which.min(xdif_v[xd_idx])
      xm_pos <- xd_idx[xd_idx2]
      
      ydif_v <- ref$uly - data[[d_cols[2]]][x]
      yd_idx <- which(ydif_v >= 0)
      # start from xm_pos
      yd_idx <- yd_idx[yd_idx >= xm_pos]
      yd_idx2 <- which.min(ydif_v[yd_idx])
      ym_pos <- yd_idx[yd_idx2]
      if (length(ym_pos)) {
        ym_pos
      } else {
        NA
      }
    }, integer(1))
    
    list(data_idx = data_idx, ref_idx = ref[idx, ]$orig_order)
  }

#' Calculate superior limits for series of sorted values
#'
#' @param x numeric vector
#' @param cnames charvec names of the columns in the output data.table
#' 
#' @keywords internal
#' 
#' @return data.table with superior limit for series of sorted values
#'
calc_up_limes <- function(x, cnames = c("ux", "ulx")) {
  
  checkmate::assert_numeric(x)
  checkmate::assert_character(cnames, min.len = 2, max.len = 2)
  
  ux <- sort(unique(x))
  tv <- ux[c(2:length(ux), length(ux))]
  ul <- ux + (tv - ux) / 2
  df <- data.table::data.table(col1 = ux, col2 = ul)
  colnames(df) <- cnames
  df
}


#' Draw interactive Drug Combo plot with drug ratios
#' 
#' @inheritParams plotlyDC1
#' 
#' @returns plotly object - line chart with lines for different ration of drugs on isolevels
#' 
#' @examples
#' \dontrun{
#' mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
#' combo_data_l <- gDRutils::convert_combo_data_to_dt(mae[[gDRutils::get_supported_experiments("combo")]])
#'
#' plotlyDCRatio(
#'   data = combo_data_l,
#'   drug1_name = "drug_005",
#'   drug2_name = "drug_021",
#'   cell_line = "cellline_HB",
#'   normalization_type = "GR",
#'   c_assay = "hsa_excess"
#' )
#'
#' plotlyDCRatio(
#'   data = combo_data_l,
#'   drug1_name = "drug_005",
#'   drug2_name = "drug_021",
#'   cell_line = "cellline_HB",
#'   normalization_type = "RV",
#'   c_assay = "smooth"
#' )
#'
#' plotlyDCRatio(
#'   data = combo_data_l,
#'   drug1_name = "drug_005",
#'   drug2_name = "drug_021",
#'   cell_line = "cellline_HB",
#'   normalization_type = "GR",
#'   c_assay = "bliss_excess",
#'   iso_levels = c("0.25", "0.5", "0.75")
#' )
#' }
#' 
#' @keywords plugin_plot
#' @export
#'
plotlyDCRatio <- function(data,
                          drug1_name,
                          drug2_name,
                          cell_line,
                          c_assay = names(gDRutils::get_combo_excess_field_names())[[1]],
                          normalization_type = "RV",
                          iso_levels =  c("0.5", "0.75"),
                          height = 325,
                          width = 425) {
  
  
  # argument checks
  checkmate::assert_list(data)
  checkmate::assert_data_table(data[[1]])
  if (is.na(drug1_name)) {
    return(NULL)
  }
  checkmate::assert_string(drug1_name)
  checkmate::assert_string(drug2_name)
  checkmate::assert_string(cell_line)
  
  checkmate::assert_string(normalization_type)
  checkmate::assert_choice(normalization_type, c("GR", "RV"))
  checkmate::assert_string(c_assay)
  cbase_anames <- gDRutils::get_combo_excess_field_names()
  checkmate::assert_choice(c_assay, names(cbase_anames))
  checkmate::assert_character(iso_levels)
  # TODO: get values from isobologram assay (iso_level) column
  gDRutils::assert_choices(iso_levels, as.character(seq(0.01, 0.99, by = 0.01)))
  
  sdata <- get_combo_base_data(
    data,
    drug1_name,
    drug2_name,
    cell_line,
    c_assay,
    normalization_type = normalization_type
  )
  
  if (is.null(sdata)) {
    return(NULL)
  }
  
  # get prettified identifiers
  pidfs <- gDRutils::get_prettified_identifiers()
  pi_drug_name <- pidfs[["drug_name"]]
  pi_drug_name2 <- pidfs[["drug_name2"]]
  pi_drug <- pidfs[["drug"]]
  pi_drug2 <- pidfs[["drug2"]]
  pi_cell_name <- pidfs[["cellline_name"]]
  pi_conc_name <- pidfs[["concentration"]]
  pi_conc2_name <- pidfs[["concentration2"]]
  pi_duration_name <- pidfs[["duration"]]
  
  
  drug1_axis <- sdata$drug1_axis
  drug2_axis <- sdata$drug2_axis
  
  condition <- sdata$condition
  
  checkmate::assert_data_table(drug1_axis)
  checkmate::assert_data_table(drug2_axis)
  checkmate::assert_character(normalization_type)
  checkmate::assert_character(condition)
  checkmate::assert_named(condition)
  checkmate::assert_names(names(condition), must.include = c(pi_drug_name, pi_drug_name2))
  
  
  matrix_pretty <- gDRutils::get_combo_excess_field_names()[c_assay]
  
  title_plot <- sprintf("%s (%s for %s, T = %sh)", condition[[pi_cell_name]], matrix_pretty,
                        normalization_type, condition[[pi_duration_name]])
  
  # transform data for isobolograms
  data_iso <- get_isobologram_data(
    data,
    drug1_name,
    drug2_name,
    cell_line,
    c_assay = names(gDRutils::get_combo_assay_names(group = "combo_iso"))[1],
    normalization_type = normalization_type,
    drug1_axis = drug1_axis,
    drug2_axis = drug2_axis,
    iso_levels = iso_levels
  )
  
  if (NROW(data_iso) == 0) { # Use NROW for cases with data_iso == NULL
    return(NULL)
  }
  
  ## extract assay components
  
  # TODO(foltynsk): can colors be more contrasting?
  iso_colors <- gDRutils::get_iso_colors(normalization_type = normalization_type)
  
  log10_ratio <- gDRutils::get_isobologram_columns("log10_ratio")
  log2_CI <- gDRutils::get_isobologram_columns("log2_CI")
  
  iso_level <- gDRutils::get_isobologram_columns("iso_level")
  
  # double check assay components
  checkmate::assert_data_table(data_iso, min.rows = 1L)
  checkmate::assert_names(names(data_iso), must.include = log10_ratio)
  checkmate::assert_character(iso_colors)
  checkmate::assert_named(iso_colors)
  checkmate::assert_string(normalization_type)
  normalization_type <- match.arg(normalization_type, choices = c("GR", "RV"), several.ok = FALSE)
  
  # axis properties
  range_x <- range(data_iso[[log10_ratio]], na.rm = TRUE)
  range_y <- c(-5.5, 4.5)
  tickvals_x <- -3:3
  tickvals_y <- -5:4
  ticktext_x <- c(paste0("1/", 10 ^ (3:1)), 10 ^ (0:3))
  ticktext_y <- c(paste0("1/", 2 ^ (5:1)), 2 ^ (0:4))
  title_x <- sprintf("%s / %s ratio", condition[pi_drug_name2], condition[pi_drug_name])
  title_y <- "CI"
  # line_properties
  line_horizontal <- list(type = "line", line = list(width = 1, color = "#000000"), layer = "above",
                          x0 = range_x[1], x1 = range_x[2], y0 = 0, y1 = 0)
  # color bar properties
  title_colorbar <- switch(normalization_type, "GR" = "GR", "RV" = "IC")
  tickvals_colorbar <- stats::quantile(as.numeric(names(iso_colors)), seq(0, 1, .25))
  ticktext_colorbar <- 100 - round(100 * (stats::quantile(as.numeric(names(iso_colors)), seq(0, 1, .25))))
  
  # build tooltips
  data_iso[["label"]] <- build_label(data_iso, "combo-ratios")
  
  data.table::setorderv(data_iso, c(log10_ratio, log2_CI))
  # build plot
  plot_base <- plotly::plot_ly(
    x = data_iso[[log10_ratio]],
    y = data_iso[[log2_CI]],
    split = data_iso[[iso_level]],
    color = data_iso[[iso_level]],
    colors = iso_colors,
    text = data_iso[["label"]],
    hoverinfo = "text",
    type = "scatter",
    mode = "lines",
    height = height,
    width = width
  )
  # apply layout
  plot_laidout <- plotly::layout(plot_base,
                                 showlegend = TRUE,
                                 title = title_plot,
                                 plot_bgcolor = "#FFF0",
                                 xaxis = list(title = title_x, zeroline = FALSE, showline = TRUE, mirror = TRUE,
                                              range = range_x, tickvals = tickvals_x, ticktext = ticktext_x),
                                 yaxis = list(title = title_y, zeroline = FALSE, showline = TRUE, mirror = TRUE,
                                              range = range_y, tickvals = tickvals_y, ticktext = ticktext_y),
                                 shapes = list(line_horizontal))
  # modify config
  plot_final <- gDR_plotly_config(plot_laidout,
                                  editable = TRUE,
                                  showAxisRangeEntryBoxes = FALSE)
  
  return(plot_final)
}

#' Round concentrations used in the heatmap
#' 
#' @keywords internal
.round_conc <- function(x) {
  x <- as.numeric(x)
  round(x, ceiling(-log10(max(x))) + 5)
}
