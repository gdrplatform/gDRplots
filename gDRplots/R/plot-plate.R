#' Plot plate data (Heatmap + Dose Ranks + Smart Legend + Explicit QC)
#' 
#' @param dt_plate \code{data.table}. Input data containing plate layout, measurements, and standard gDR identifiers.
#' @param ctrl_fail_threshold \code{numeric} (Default: 0.6). Flags control wells with a readout below
#' \code{Mean(Controls) * ctrl_fail_threshold}.
#' @param n_sd \code{numeric} (Default: 1). Sets the upper limit for flagging high signals using
#' the formula \code{Mean(Controls) + (n_sd * SD(Controls))}.
#' @param use_sd_threshold \code{logical} (Default: TRUE). If \code{TRUE}, uses the dynamic
#' SD-based limit defined by \code{n_sd}. If \code{FALSE}, uses a fixed limit of \code{Mean(Controls) * 1.1}.
#' @param items_per_line \code{integer}. Number of doses to show per line in the legend. Default 6.
#' @return A named list of ggplot objects for each barcode
#' @keywords QC_plot
#' 
#' @author Bartosz Czech \email{bartosz.czech@@contractors.roche.com}
#' 
#' @export
#' @examples
#' conc_series <- c(0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10)
#' 
#' test_data <- data.table::data.table(
#'    WellColumn = sprintf("%02d", rep(1:12, each = 8)),
#'    WellRow = rep(LETTERS[1:8], times = 12),
#'    Barcode = rep(c("Plate_1", "Plate_2"), each = 48),
#'    clid = "CellLineA"
#' )
#' 
#' invisible(test_data[, Concentration := rep(rep(conc_series, length.out = 48), 2)])
#' invisible(test_data[, Gnumber := ifelse(Concentration == 0, "vehicle", "Drug_A")])
#' invisible(test_data[Gnumber == "vehicle", Concentration := 0])
#' 
#' invisible(test_data[, ReadoutValue := ifelse(Gnumber == "vehicle", 
#'                                     rnorm(.N, 1000, 50), 
#'                                     rnorm(.N, 500, 100))])
#' library(ggtext)
#' plots <- plot_plate_stack_info(test_data)
#' plots[[1]]
#' 
#' combo_data <- data.table::copy(test_data)
#' invisible(combo_data[, Barcode := "Plate_Combo"])
#' invisible(combo_data[, Concentration_2 := rep(rev(conc_series), length.out = .N)])
#' invisible(combo_data[, Gnumber_2 := ifelse(Concentration_2 == 0, "vehicle", "Drug_B")])
#' invisible(combo_data[Gnumber_2 == "vehicle", Concentration_2 := 0])
#' 
#' combo_plots <- plot_plate_stack_info(combo_data)
#' combo_plots[[1]]
#'
plot_plate_stack_info <- function(dt_plate, 
                                  ctrl_fail_threshold = 0.6, 
                                  n_sd = 1,
                                  use_sd_threshold = TRUE,
                                  items_per_line = 6) {
  
  checkmate::assert_data_table(dt_plate)
  checkmate::assert_number(ctrl_fail_threshold, lower = 0, upper = 1)
  checkmate::assert_number(n_sd, lower = 0)
  checkmate::assert_flag(use_sd_threshold)
  checkmate::assert_int(items_per_line, lower = 1)
  
  barcode <- gDRutils::get_env_identifiers("barcode")[1]
  
  req_cols <- c("ReadoutValue", "WellRow", "WellColumn", barcode)
  checkmate::assert_names(names(dt_plate), must.include = req_cols)
  
  unique_plates <- unique(dt_plate[[barcode]])
  
  plate_list <- gDRutils::loop(unique_plates, function(pid) {
    dt_subset <- dt_plate[get(barcode) == pid]
    
    plot_single_plate_stack_info(
      dt_subset = dt_subset,
      plate_id = pid,
      ctrl_fail_threshold = ctrl_fail_threshold,
      n_sd = n_sd,
      use_sd_threshold = use_sd_threshold,
      items_per_line = items_per_line
    )
  })
  
  names(plate_list) <- unique_plates
  return(plate_list)
}

#' Plot a single plate's stack info
#' 
#' Generates the QC and dose-rank heatmap for a single plate subset.
#' Handles all necessary data transformations (factor conversion, ranking, color mapping) internally.
#' 
#' @param dt_subset \code{data.table}. Subset of the plate data corresponding to a single barcode/plate.
#' @param plate_id \code{character} (Optional). The identifier (barcode) of the plate.
#' If \code{NULL}, attempts to extract it from the data.
#' @param drug_color_mapping \code{character} (Optional). Named vector of color codes.
#' If \code{NULL}, generated automatically.
#' @inheritParams plot_plate_stack_info
#' 
#' @return A \code{ggplot} object.
#' @keywords QC_plot
#' 
#' @author Bartosz Czech \email{bartosz.czech@@contractors.roche.com}
#' 
#' @export
#' @examples
#' conc_series <- c(0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10)
#' 
#' test_data <- data.table::data.table(
#'    WellColumn = rep(1:12, each = 8),
#'    WellRow = rep(LETTERS[1:8], times = 12),
#'    clid = "CellLineA",
#'    Barcode = "Plate_1"
#' )
#' 
#' invisible(test_data[, Concentration := rep(conc_series, length.out = .N)])
#' invisible(test_data[, Gnumber := ifelse(Concentration == 0, "vehicle", "Drug_A")])
#' invisible(test_data[Gnumber == "vehicle", Concentration := 0])
#' 
#' invisible(test_data[, ReadoutValue := ifelse(Gnumber == "vehicle", 
#'                                     rnorm(.N, 1000, 50), 
#'                                     rnorm(.N, 1000 * (1 / (1 + Concentration)), 50))])
#'                                      
#' library(ggtext)
#' plot_single_plate_stack_info(test_data)
#' 
#' invisible(test_data[, Concentration_2 := rep(rev(conc_series), length.out = .N)])
#' invisible(test_data[, Gnumber_2 := ifelse(Concentration_2 == 0, "vehicle", "Drug_B")])
#' 
#' invisible(test_data[Gnumber_2 == "vehicle", Concentration_2 := 0])
#' 
#' plot_single_plate_stack_info(test_data)
#'
plot_single_plate_stack_info <- function(dt_subset,
                                         plate_id = NULL,
                                         drug_color_mapping = NULL,
                                         ctrl_fail_threshold = 0.6,
                                         n_sd = 1,
                                         use_sd_threshold = TRUE,
                                         items_per_line = 6) {
  
  checkmate::assert_data_table(dt_subset)
  
  dt_subset <- data.table::copy(dt_subset)
  
  concentration <- gDRutils::get_env_identifiers("concentration")
  concentration2 <- gDRutils::get_env_identifiers("concentration2")
  drug <- gDRutils::get_env_identifiers("drug")
  drug2 <- gDRutils::get_env_identifiers("drug2")
  cellline <- gDRutils::get_env_identifiers("cellline")
  barcode <- gDRutils::get_env_identifiers("barcode")[1]
  
  if (is.null(plate_id)) {
    if (barcode %in% names(dt_subset)) {
      plate_id <- as.character(unique(dt_subset[[barcode]])[1])
    } else {
      plate_id <- "Unknown"
    }
  }
  
  if (is.character(dt_subset$WellColumn) || is.factor(dt_subset$WellColumn)) {
    dt_subset[, WellColumn := as.numeric(as.character(WellColumn))]
  }
  
  row_levels <- sort(unique(as.character(dt_subset$WellRow)), decreasing = TRUE)
  dt_subset[, WellRow := factor(as.character(WellRow), levels = row_levels)]
  
  if (!"rank_1" %in% names(dt_subset)) {
    dt_subset[, rank_1 := .calc_dose_rank(get(concentration))]
  }
  
  has_combo <- concentration2 %in% names(dt_subset) && drug2 %in% names(dt_subset)
  if (has_combo) {
    if (!"rank_2" %in% names(dt_subset)) {
      dt_subset[, rank_2 := .calc_dose_rank(get(concentration2))]
    }
  }
  
  if (is.null(drug_color_mapping)) {
    drug_color_mapping <- generate_color_mappings(dt_subset)
  }
  
  solid_shapes <- c(16, 17, 15, 18, 19, 20, 8) 
  
  is_ctrl <- dt_subset[[concentration]] == 0 
  if (concentration2 %in% names(dt_subset)) {
    is_ctrl <- is_ctrl & (dt_subset[[concentration2]] == 0 | is.na(dt_subset[[concentration2]]))
  }
  
  dt_controls <- dt_subset[is_ctrl]
  dt_treated <- dt_subset[!is_ctrl]
  
  mean_ctrl <- mean(dt_controls$ReadoutValue, na.rm = TRUE)
  sd_ctrl <- stats::sd(dt_controls$ReadoutValue, na.rm = TRUE)
  qc_valid <- is.finite(mean_ctrl)
  
  if (use_sd_threshold) {
    limit_desc <- paste0("Mean(Ctrl) + ", n_sd, "*SD")
  } else {
    limit_desc <- "Mean(Ctrl) + 10%"
  }
  low_desc <- paste0("Mean(Ctrl) * ", ctrl_fail_threshold)
  
  if (qc_valid) {
    if (is.na(sd_ctrl)) {
      sd_ctrl <- 0
    }
    
    if (use_sd_threshold) {
      limit_value <- mean_ctrl + (n_sd * sd_ctrl)
    } else {
      limit_value <- mean_ctrl * 1.1
    }
    
    low_limit_value <- mean_ctrl * ctrl_fail_threshold
    
    bad_ctrl_rows <- dt_controls[ReadoutValue < low_limit_value]
    suspicious_treated <- dt_treated[ReadoutValue > limit_value]
    
    count_low <- NROW(bad_ctrl_rows)
    count_high <- NROW(suspicious_treated)
  } else {
    count_low <- 0
    count_high <- 0
    limit_value <- NA
    low_limit_value <- NA
    bad_ctrl_rows <- dt_controls[0]
    suspicious_treated <- dt_treated[0]
  }
  
  doses_1 <- sort(unique(dt_subset[[concentration]][dt_subset[[concentration]] > 0]))
  dose_key_str <- .format_dose_list(doses_1, items_per_line)
  dose_legend <- paste0("<b>Dose Key (inside dots):</b><br>", dose_key_str)
  
  if (has_combo) {
    doses_2 <- sort(unique(dt_subset[[concentration2]][dt_subset[[concentration2]] > 0]))
    if (!identical(doses_1, doses_2) && length(doses_2) > 0) {
      dose_key_str_2 <- .format_dose_list(doses_2, items_per_line)
      dose_legend <- paste0(dose_legend, "<br><b>Drug 2 Key:</b><br>", dose_key_str_2)
    }
  }
  
  if (!qc_valid) {
    subtitle_text <- "<span style='color:orange;'>QC Status: <b>No Controls</b></span>"
    qc_legend <- "N/A"
  } else {
    if (count_low == 0 && count_high == 0) {
      subtitle_text <- "<span style='color:darkgreen;'>QC Status: <b>OK</b></span>"
    } else {
      msgs <- c()
      if (count_low > 0) {
        msgs <- c(msgs, paste0("<b>", count_low, "</b> Low"))
      }
      if (count_high > 0) {
        msgs <- c(msgs, paste0("<b>", count_high, "</b> High"))
      }
      subtitle_text <- paste0("<span style='color:red;'>QC ALERTS: ", paste(msgs, collapse = " | "), "</span>")
    }
    
    qc_legend <- paste0(
      "<span style='color:red;'>&#9632;</span> <b>Low</b> ",
      "<span style='color:#555; font-size:8pt;'>(< ", round(low_limit_value, 1), " [", low_desc, "])</span> | ",
      "<span style='color:darkorange;'>&#9632;</span> <b>High</b> ",
      "<span style='color:#555; font-size:8pt;'>(> ", round(limit_value, 1), " [", limit_desc, "])</span>"
    )
  }
  
  caption_text <- paste0(dose_legend, "<br><br><b>QC Flags:</b> ", qc_legend)
  
  if (has_combo) {
    offset_1 <- -0.18
  } else {
    offset_1 <- 0
  }
  offset_2 <- 0.18
  
  combo_layers <- NULL
  if (has_combo) {
    combo_layers <- list(
      ggplot2::geom_point(ggplot2::aes(
        x = WellColumn + offset_2, y = WellRow, color = !!rlang::sym(drug2)
      ), size = 4),
      ggplot2::geom_text(ggplot2::aes(
        x = WellColumn + offset_2, y = WellRow, label = rank_2
      ), color = "white", size = 2.5, fontface = "bold")
    )
  }
  
  p <- ggplot2::ggplot(dt_subset) +
    ggplot2::geom_tile(ggplot2::aes(x = WellColumn, y = WellRow, fill = ReadoutValue)) +
    ggplot2::geom_point(ggplot2::aes(
      x = WellColumn + offset_1, y = WellRow,
      color = !!rlang::sym(drug), shape = !!rlang::sym(cellline)
    ), size = 4) +
    ggplot2::geom_text(ggplot2::aes(
      x = WellColumn + offset_1, y = WellRow, label = rank_1
    ), color = "white", size = 2.5, fontface = "bold") +
    combo_layers + 
    ggplot2::geom_tile(ggplot2::aes(x = WellColumn, y = WellRow), fill = NA, color = "black", linewidth = 0.5) +
    ggplot2::geom_tile(data = bad_ctrl_rows, ggplot2::aes(x = WellColumn, y = WellRow),
                       color = "red", fill = NA, linewidth = 1, width = 1, height = 1) +
    ggplot2::geom_tile(data = suspicious_treated, ggplot2::aes(x = WellColumn, y = WellRow),
                       color = "darkorange", fill = NA, linewidth = 1, width = 1, height = 1) +
    ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#6FA8DC", name = "Readout") +
    ggplot2::scale_x_continuous(breaks = sort(unique(dt_subset$WellColumn)), 
                                labels = sort(unique(dt_subset$WellColumn)), position = "top") +
    ggplot2::scale_color_manual(values = drug_color_mapping) +
    ggplot2::scale_shape_manual(values = solid_shapes) +
    ggplot2::labs(
      title = paste0("Plate ", plate_id, " QC Analysis"), 
      subtitle = subtitle_text,
      caption = caption_text,
      x = "Column", y = "Row"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   plot.subtitle = ggtext::element_markdown(hjust = 0.5, size = 11),
                   plot.caption = ggtext::element_markdown(hjust = 0, size = 8, color = "#444444", lineheight = 1.2),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank())
  
  return(p)
}


#' Plot data from a specific column
#'
#' @param dt_plate The data table containing plate information
#' @param column_name The name of the column to plot
#'
#' @return A named list of ggplot objects for each barcode and column
#' 
#' @keywords QC_plot
#' 
#' @author Bartosz Czech \email{bartosz.czech@@contractors.roche.com}
#' 
#' @examples
#' test_data <- data.table::data.table(
#'   WellColumn = rep(1:12, each = 8),
#'   WellRow = rep(LETTERS[1:8], times = 12),
#'   clid = "A",
#'   Gnumber = c(rep("untreated", 48), sample(1:5, size = 48, replace = TRUE)),
#'   Gnumber_2 = c(rep("untreated", 48), sample(1:5, size = 48, replace = TRUE)),
#'   Concentration = runif(96, min = 0, max = 100),
#'   ReadoutValue = runif(96, min = 0, max = 100),
#'   Barcode = rep(c("A", "B"), 48)
#'   )
#' plot_plate(test_data, "Gnumber")[["A"]]
#' @export
plot_plate <- function(dt_plate, column_name) {
  
  dt_plate_copy <- data.table::copy(dt_plate)
  concentration <- gDRutils::get_env_identifiers("concentration")
  concentration2 <- gDRutils::get_env_identifiers("concentration2")
  concentration3 <- gDRutils::get_env_identifiers("concentration3")
  drug <- gDRutils::get_env_identifiers("drug")
  drug2 <- gDRutils::get_env_identifiers("drug2")
  cellline <- gDRutils::get_env_identifiers("cellline")
  well_position <- gDRutils::get_env_identifiers("well_position")
  barcode <- gDRutils::get_env_identifiers("barcode")
  
  plate_palette <- gDRutils::get_settings_from_json(
    "PLATE_PALETTE",
    system.file(package = "gDRplots", "settings.json"))
  
  checkmate::assert_data_table(dt_plate_copy)
  checkmate::assert_names(names(dt_plate_copy),
                          must.include = c("ReadoutValue",
                                           concentration,
                                           drug,
                                           cellline,
                                           well_position))
  checkmate::assert_choice(column_name, choices = names(dt_plate_copy))
  
  dt_plate_copy[, WellColumn := as.numeric(WellColumn)]
  dt_plate_copy[, WellRow := factor(WellRow, levels = sort(unique(WellRow), decreasing = TRUE))]
  
  barcode_idf <- intersect(barcode, names(dt_plate_copy))
  if (NROW(barcode_idf) == 0) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }
  
  continuous <- is.numeric(dt_plate_copy[[column_name]]) &&
    !column_name %in% c(concentration, concentration2, concentration3)
  
  dt_plate_copy[, (column_name) := {
    col <- .SD[[column_name]]
    if (is.numeric(col) && !column_name %in% c(concentration, concentration2, concentration3)) {
      round(col, 5)
    } else if (is.numeric(col)) {
      factor(round(col, 5))
    } else {
      factor(col)
    }
  }]
  
  plate_list <- gDRutils::loop(unique(dt_plate_copy[[barcode_idf]]), function(x) {
    dt_plate_copy_subset <- dt_plate_copy[get(barcode_idf) == x]
    
    dt_plate_copy_subset[, WellRow := droplevels(WellRow)]
    if (!continuous) {
      dt_plate_copy_subset[, (column_name) := droplevels(get(column_name))]
    }
    
    unique_val <- sort(unique(dt_plate_copy_subset[[column_name]]))
    
    if (!continuous) {
      gradient_colors <- grDevices::colorRampPalette(plate_palette)(length(unique_val))
      names(gradient_colors) <- unique_val
    }
    
    p <- ggplot2::ggplot(dt_plate_copy_subset) +
      ggplot2::geom_rect(ggplot2::aes(xmin = WellColumn - 0.5, xmax = WellColumn + 0.5,
                                      ymin = as.numeric(WellRow) - 0.5,
                                      ymax = as.numeric(WellRow) + 0.5,
                                      fill = !!rlang::sym(column_name)),
                         color = "black", linewidth = 0.2) +
      ggrepel::geom_text_repel(data = dt_plate_copy_subset, 
                               ggplot2::aes(x = WellColumn, y = WellRow, label = !!rlang::sym(column_name)), 
                               color = "black", size = 3, bg.color = "white", bg.r = 0.08, force = 0) +
      ggplot2::scale_x_continuous(breaks = sort(unique(dt_plate_copy_subset$WellColumn)),
                                  labels = sort(unique(dt_plate_copy_subset$WellColumn)),
                                  position = "top") +
      ggplot2::labs(
        title = paste0(length(unique(dt_plate_copy_subset$WellColumn)) * length(unique(dt_plate_copy_subset$WellRow)),
                       "-Well Visualization of Plate ", x, " for ", column_name),
        x = "Column",
        y = "Row"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
        plot.title = ggplot2::element_text(hjust = 0.5),
        plot.subtitle = ggplot2::element_text(hjust = 0.5),
        panel.grid.minor = ggplot2::element_blank()
      )
    
    
    if (continuous) {
      p <- p + ggplot2::scale_fill_gradient(low =  plate_palette[[1]],
                                            high = plate_palette[[2]])
    } else {
      p <- p + ggplot2::scale_fill_manual(values = gradient_colors,
                                          name = column_name,
                                          limits = names(gradient_colors))
    }
    
    p + ggplot2::geom_tile(ggplot2::aes(x = WellColumn, y = WellRow),
                           fill = NA, color = "black", linewidth = 0.5)
    
  })
  
  names(plate_list) <- unique(dt_plate_copy[[barcode_idf]])
  plate_list
}


#' Helper function to generate color mappings for Gnumber and Gnumber_2
#' 
#' @param dt_plate_subset The subset of data
#' @param untrt_tag Untreated tag identifier
#' 
#' @keywords internal
#' 
#' @author Bartosz Czech \email{bartosz.czech@@contractors.roche.com}
#' 
#' @return A list containing color mappings for Gnumber and Gnumber_2
generate_color_mappings <- function(dt_plate_subset,
                                    untrt_tag = gDRutils::get_env_identifiers("untreated_tag")) {
  checkmate::assert_data_table(dt_plate_subset)
  checkmate::assert_character(untrt_tag)
  special_colors <- rep("darkgray", length(untrt_tag))
  names(special_colors) <- untrt_tag
  
  drug <- gDRutils::get_env_identifiers("drug")
  drug2 <- gDRutils::get_env_identifiers("drug2")
  drug3 <- gDRutils::get_env_identifiers("drug3")
  
  other_gnumbers <- setdiff(unique(unlist(dt_plate_subset[, intersect(names(dt_plate_subset),
                                                                      c(drug, drug2, drug3)), with = FALSE])),
                            untrt_tag)
  if (NROW(other_gnumbers) > 0) {
    palette_colors <- scales::hue_pal()(length(other_gnumbers))
    names(palette_colors) <- other_gnumbers
    c(special_colors, palette_colors)
  } else {
    special_colors 
  }
}

### helpers

#' Calculate Dose Ranks
#'
#' Helper function to rank unique positive concentration values.
#' Zero values or NAs are assigned a rank of "-".
#'
#' @param vals \code{numeric} vector of concentration values.
#' @return A \code{character} vector of ranks.
#' @keywords internal
#' @author Bartosz Czech \email{bartosz.czech@@contractors.roche.com}
.calc_dose_rank <- function(vals) {
  vals <- round(vals, 6)
  unique_doses <- sort(unique(vals[vals > 0]))
  ranks <- character(length(vals))
  ranks[vals == 0 | is.na(vals)] <- "-" 
  
  if (length(unique_doses) > 0) {
    positive_indices <- vals > 0 & !is.na(vals)
    ranks[positive_indices] <- match(vals[positive_indices], unique_doses)
  }
  return(ranks)
}

#' Format Dose List for Legend
#'
#' Helper function to create a formatted string of dose ranks and values 
#' for the plot legend, breaking lines based on `n_per_line`.
#'
#' @param dose_vec \code{numeric} vector of unique positive doses.
#' @param n_per_line \code{integer} number of items to display per line.
#' @return A single \code{character} string with HTML line breaks.
#' @keywords internal
#' @author Bartosz Czech \email{bartosz.czech@@contractors.roche.com}
.format_dose_list <- function(dose_vec, n_per_line) {
  if (length(dose_vec) == 0) return("")
  
  lbls <- paste(seq_along(dose_vec), signif(dose_vec, 3), sep = ": ")
  grp <- ceiling(seq_along(lbls) / n_per_line)
  lines <- tapply(lbls, grp, paste, collapse = " | ")
  
  paste(lines, collapse = "<br>")
}
