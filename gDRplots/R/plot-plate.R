#' Plot plate data
#' 
#' @param dt_plate The data table containing plate information
#' @param ctrl_fail_threshold Numeric (0-1). Flags controls below `Mean(Ctrl) * this`. Default 0.6.
#' @param n_sd Numeric. Number of SDs above mean to flag high signals. Default 1 (High Sensitivity).
#' @param use_sd_threshold Logical. If TRUE, uses dynamic limit (`Mean + n_sd*SD`). If FALSE, uses fixed limit (`Mean * 1.1`).
#' @return A named list of ggplot objects for each barcode
#' @keywords QC_plot
#' @examples
#' test_data <- data.table::data.table(
#'   WellColumn = rep(1:12, each = 8),
#'   clid = c("A", "A", "B", "B"),
#'   WellRow = rep(LETTERS[1:8], times = 12),
#'   Gnumber = c(rep("untreated", 48), sample(1:5, size = 48, replace = TRUE)),
#'   Gnumber_2 = c(rep("untreated", 48), sample(1:5, size = 48, replace = TRUE)),
#'   Concentration = sample(seq(0, 2.5, 1.25), 12, replace = TRUE),
#'   Concentration_2 = sample(seq(0, 2.5, 1.25), 12, replace = TRUE),
#'   ReadoutValue = runif(96, min = 0, max = 100),
#'   Barcode = rep(c("A", "B"), 48)
#'   )
#' plot_plate_stack_info(test_data)[["A"]]
#' @export
plot_plate_stack_info <- function(dt_plate, 
                                     ctrl_fail_threshold = 0.6, 
                                     n_sd = 1,
                                     use_sd_threshold = TRUE) {
  dt_plate_copy <- data.table::copy(dt_plate)
  
  concentration <- gDRutils::get_env_identifiers("concentration")
  concentration2 <- gDRutils::get_env_identifiers("concentration2")
  drug <- gDRutils::get_env_identifiers("drug")
  drug2 <- gDRutils::get_env_identifiers("drug2")
  cellline <- gDRutils::get_env_identifiers("cellline")
  barcode <- gDRutils::get_env_identifiers("barcode")
  
  plate_palette <- gDRutils::get_settings_from_json("PLATE_PALETTE", system.file(package = "gDRplots", "settings.json"))
  
  checkmate::assert_data_table(dt_plate_copy)
  dt_plate_copy[, WellColumn := as.numeric(WellColumn)]
  dt_plate_copy[, WellRow := factor(WellRow, levels = sort(unique(WellRow), decreasing = TRUE))]
  
  barcode_idf <- intersect(barcode, names(dt_plate_copy))
  if (NROW(barcode_idf) == 0) return(ggplot2::ggplot() + ggplot2::theme_void())
  
  overall_color_mapping <- generate_color_mappings(dt_plate_copy) 
  
  plate_list <- gDRutils::loop(unique(dt_plate_copy[[barcode_idf]]), function(x) {
    dt_plate_copy_subset <- dt_plate_copy[get(barcode_idf) == x]
    dt_plate_copy_subset[, WellRow := droplevels(WellRow)]
    
    is_ctrl <- dt_plate_copy_subset[[concentration]] == 0 
    if (concentration2 %in% names(dt_plate_copy_subset)) {
      is_ctrl <- is_ctrl & (dt_plate_copy_subset[[concentration2]] == 0 | is.na(dt_plate_copy_subset[[concentration2]]))
    }
    
    dt_controls <- dt_plate_copy_subset[is_ctrl]
    dt_treated  <- dt_plate_copy_subset[!is_ctrl]
    
    mean_ctrl <- mean(dt_controls$ReadoutValue, na.rm = TRUE)
    sd_ctrl   <- sd(dt_controls$ReadoutValue, na.rm = TRUE)
    if (is.na(sd_ctrl)) sd_ctrl <- 0
    
    if (use_sd_threshold) {
      limit_value <- mean_ctrl + (n_sd * sd_ctrl)
      high_formula_desc <- paste0("Mean_Ctrl + ", n_sd, "*SD")
    } else {
      limit_value <- mean_ctrl * 1.1 
      high_formula_desc <- "Mean_Ctrl + 10%"
    }
    
    low_limit_value <- mean_ctrl * ctrl_fail_threshold
    low_formula_desc <- paste0("Mean_Ctrl * ", ctrl_fail_threshold)
    
    bad_ctrl_rows <- dt_controls[ReadoutValue < low_limit_value]
    suspicious_treated <- dt_treated[ReadoutValue > limit_value]
    
    count_low <- nrow(bad_ctrl_rows)
    count_high <- nrow(suspicious_treated)
    
    if (count_low == 0 && count_high == 0) {
      subtitle_text <- "<span style='color:darkgreen;'>QC Status: <b>OK</b></span>"
    } else {
      msgs <- c()
      if (count_low > 0) msgs <- c(msgs, paste0("<b>", count_low, "</b> Low Ctrls"))
      if (count_high > 0) msgs <- c(msgs, paste0("<b>", count_high, "</b> High Artifacts"))
      subtitle_text <- paste0("<span style='color:red;'>QC ALERTS: ", paste(msgs, collapse = " | "), "</span>")
    }
    
    caption_text <- paste0(
      "<b>QC DIAGNOSTICS:</b><br>",
      "<span style='color:red; font-size:14pt;'>&#9632;</span> ",
      "<b>Low Control</b> (Edge Effect): value < ", round(low_limit_value, 1), 
      " <span style='color:#777777; font-size:9pt;'>(", low_formula_desc, ")</span>",
      "<br>",
      "<span style='color:darkorange; font-size:14pt;'>&#9632;</span> ",
      "<b>High Signal</b> (Artifact): value > ", round(limit_value, 1), 
      " <span style='color:#777777; font-size:9pt;'>(", high_formula_desc, ")</span>"
    )
    
    has_combo <- concentration2 %in% colnames(dt_plate_copy_subset) && drug2 %in% colnames(dt_plate_copy_subset)
    doses <- sort(unique(unlist(dt_plate_copy_subset[, intersect(names(dt_plate_copy_subset), c(concentration, concentration2)), with = FALSE])))
    gradient_colors <- grDevices::colorRampPalette(plate_palette)(length(doses))
    names(gradient_colors) <- doses
    
    dt_plate_copy_subset[, (concentration) := factor(get(concentration), levels = doses)]
    if (has_combo) dt_plate_copy_subset[, (concentration2) := factor(get(concentration2), levels = doses)]
    
    xmax_offset <- if (has_combo) 0 else 0.5
    x_point_offset <- if (has_combo) -0.25 else 0
    
    p <- ggplot2::ggplot(dt_plate_copy_subset) +
      ggplot2::geom_rect(ggplot2::aes(
        xmin = WellColumn - 0.5, xmax = WellColumn + xmax_offset,
        ymin = as.numeric(WellRow) - 0.5, ymax = as.numeric(WellRow) + 0.5,
        fill = !!rlang::sym(concentration)
      ), color = NA) +
      ggplot2::geom_point(ggplot2::aes(
        x = WellColumn + x_point_offset, y = WellRow,
        color = !!rlang::sym(drug), shape = !!rlang::sym(cellline)
      ), size = 3)
    
    if (has_combo) {
      p <- p + ggplot2::geom_rect(ggplot2::aes(
        xmin = WellColumn, xmax = WellColumn + 0.5,
        ymin = as.numeric(WellRow) - 0.5, ymax = as.numeric(WellRow) + 0.5,
        fill = !!rlang::sym(concentration2)
      ), color = NA) +
        ggplot2::geom_point(ggplot2::aes(x = WellColumn + 0.25, y = WellRow, color = !!rlang::sym(drug2)), size = 3)
    }
    
    p <- p +
      ggrepel::geom_text_repel(
        data = dt_plate_copy_subset,
        ggplot2::aes(x = WellColumn, y = WellRow, label = round(ReadoutValue, 1)),
        color = "black", size = 3, bg.color = "white", bg.r = 0.05,
        force = 0, nudge_y = -0.3, segment.color = NA
      ) +
      ggplot2::geom_tile(ggplot2::aes(x = WellColumn, y = WellRow), fill = NA, color = "black", linewidth = 0.5) +
      ggplot2::geom_tile(data = bad_ctrl_rows, ggplot2::aes(x = WellColumn, y = WellRow),
                         color = "red", fill = NA, linewidth = 0.75, width = 1, height = 1) +
      ggplot2::geom_tile(data = suspicious_treated, ggplot2::aes(x = WellColumn, y = WellRow),
                         color = "darkorange", fill = NA, linewidth = 0.75, width = 1, height = 1) +
      
      ggplot2::scale_fill_manual(values = gradient_colors, name = concentration, limits = names(gradient_colors)) +
      ggplot2::scale_x_continuous(breaks = sort(unique(dt_plate_copy_subset$WellColumn)),
                                  labels = sort(unique(dt_plate_copy_subset$WellColumn)), position = "top") +
      ggplot2::scale_color_manual(values = overall_color_mapping) +
      
      ggplot2::labs(
        title = paste0("Plate ", x, " QC Analysis"), 
        subtitle = subtitle_text,
        caption = caption_text,
        x = "Column", y = "Row"
      ) +
      
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                     plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                     plot.subtitle = ggtext::element_markdown(hjust = 0.5, size = 11),
                     plot.caption = ggtext::element_markdown(hjust = 0, size = 9, color = "#333333", lineheight = 1.3),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank()) +
      
      ggplot2::scale_shape_manual(values = scales::shape_pal()(length(unique(dt_plate_copy[[cellline]]))))
    
    p
  })
  names(plate_list) <- unique(dt_plate_copy[[barcode_idf]])
  plate_list
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
