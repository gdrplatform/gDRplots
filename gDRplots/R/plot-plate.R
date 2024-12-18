#' Plot plate data
#' 
#' @param dt_plate The data table containing plate information
#' 
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
plot_plate_stack_info <- function(dt_plate) {
  concentration <- gDRutils::get_env_identifiers("concentration")
  concentration2 <- gDRutils::get_env_identifiers("concentration2")
  drug <- gDRutils::get_env_identifiers("drug")
  drug2 <- gDRutils::get_env_identifiers("drug2")
  cellline <- gDRutils::get_env_identifiers("cellline")
  well_position <- gDRutils::get_env_identifiers("well_position")
  barcode <- gDRutils::get_env_identifiers("barcode")
  
  checkmate::assert_data_table(dt_plate)
  checkmate::assert_names(names(dt_plate),
                          must.include = c("ReadoutValue",
                                           concentration,
                                           drug,
                                           cellline,
                                           well_position))
  
  dt_plate[, WellColumn := as.numeric(WellColumn)]
  dt_plate[, WellRow := factor(WellRow, levels = sort(unique(WellRow), decreasing = TRUE))]
  
  barcode_idf <- intersect(barcode, names(dt_plate))
  if (NROW(barcode_idf) == 0) {
    return(ggplot2::ggplot())
  }
  
  overall_color_mapping <- generate_color_mappings(dt_plate)
  
  plate_list <- lapply(unique(dt_plate[[barcode_idf]]), function(x) {
    dt_plate_subset <- dt_plate[get(barcode_idf) == x]
    
    dt_plate_subset[, WellRow := droplevels(WellRow)]
    
    has_combo <- concentration2 %in% colnames(dt_plate_subset) && drug2 %in% colnames(dt_plate_subset)
    
    doses <- sort(unique(unlist(dt_plate_subset[, intersect(names(dt_plate_subset),
                                                            c(concentration, concentration2)), with = FALSE])))
    gradient_colors <- grDevices::colorRampPalette(c("#c6dbef", "white", "#08306b"))(length(doses))
    names(gradient_colors) <- doses
    
    dt_plate_subset[, (concentration) := factor(get(concentration), levels = doses)]
    if (has_combo) {
      dt_plate_subset[, (concentration2) := factor(get(concentration2), levels = doses)]
    }
    
    
    p <- ggplot2::ggplot(dt_plate_subset) +
      ggplot2::geom_rect(ggplot2::aes(
        xmin = WellColumn - 0.5,
        xmax = WellColumn + ifelse(has_combo, 0, 0.5),
        ymin = as.numeric(WellRow) - 0.5,
        ymax = as.numeric(WellRow) + 0.5,
        fill = !!rlang::sym(concentration)
      ),
      color = "black", linewidth = 0.2
      ) + {
        if (has_combo) ggplot2::geom_rect(ggplot2::aes(
        xmin = WellColumn,
        xmax = WellColumn + 0.5,
        ymin = as.numeric(WellRow) - 0.5,
        ymax = as.numeric(WellRow) + 0.5,
        fill = !!rlang::sym(concentration2)
      ),
      color = "black", linewidth = 0.2
      )
        }
    +
      ggplot2::geom_point(ggplot2::aes(
        x = WellColumn + ifelse(has_combo, -0.25, 0),
        y = WellRow,
        color = !!rlang::sym(drug),
        shape = !!rlang::sym(cellline)
      ),
      size = 3
      ) + {
        if (has_combo) ggplot2::geom_point(ggplot2::aes(
        x = WellColumn + 0.25,
        y = WellRow,
        color = !!rlang::sym(drug2)
      ),
      size = 3
      )
        }
    +
      ggplot2::geom_text(ggplot2::aes(
        x = WellColumn, y = WellRow, label = round(ReadoutValue, 1)
      ),
      color = "black", size = 3, vjust = 1.2
      ) +
      ggplot2::scale_fill_manual(values = gradient_colors, name = concentration, limits = names(gradient_colors)) +
      ggplot2::scale_x_continuous(breaks = sort(unique(dt_plate_subset$WellColumn)),
                                  labels = sort(unique(dt_plate_subset$WellColumn)),
                                  position = "top") +
      ggplot2::labs(
        title = paste0(length(unique(dt_plate_subset$WellColumn)) * length(unique(dt_plate_subset$WellRow)),
                       "-Well Visualization of Plate ", x),
        x = "Column",
        y = "Row"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
        plot.title = ggplot2::element_text(hjust = 0.5),
        plot.subtitle = ggplot2::element_text(hjust = 0.5),
        panel.grid.minor = ggplot2::element_blank()
      ) +
      ggplot2::scale_color_manual(values = overall_color_mapping) +
      ggplot2::scale_shape_manual(values = scales::shape_pal()(length(unique(dt_plate[[cellline]])))) +
      ggplot2::geom_tile(ggplot2::aes(x = WellColumn, y = WellRow), fill = NA, color = "black", linewidth = 0.5)
    
  })
  
  names(plate_list) <- unique(dt_plate[[barcode_idf]])
  plate_list
}


#' Plot data from a specific column
#'
#' @param dt_plate The data table containing plate information
#' @param column_name The name of the column to plot
#'
#' @return A named list of ggplot objects for each barcode and column
#' @keywords QC_plot
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
  concentration <- gDRutils::get_env_identifiers("concentration")
  concentration2 <- gDRutils::get_env_identifiers("concentration2")
  concentration3 <- gDRutils::get_env_identifiers("concentration3")
  drug <- gDRutils::get_env_identifiers("drug")
  drug2 <- gDRutils::get_env_identifiers("drug2")
  cellline <- gDRutils::get_env_identifiers("cellline")
  well_position <- gDRutils::get_env_identifiers("well_position")
  barcode <- gDRutils::get_env_identifiers("barcode")
  
  checkmate::assert_data_table(dt_plate)
  checkmate::assert_names(names(dt_plate),
                          must.include = c("ReadoutValue",
                                           concentration,
                                           drug,
                                           cellline,
                                           well_position))
  checkmate::assert_choice(column_name, choices = names(dt_plate))
  
  dt_plate[, WellColumn := as.numeric(WellColumn)]
  dt_plate[, WellRow := factor(WellRow, levels = sort(unique(WellRow), decreasing = TRUE))]
  
  barcode_idf <- intersect(barcode, names(dt_plate))
  if (NROW(barcode_idf) == 0) {
    return(ggplot2::ggplot())
  }
  
  unique_val <- sort(unique(dt_plate[[column_name]]))
  continuous <- is.numeric(unique_val) && !column_name %in% c(concentration, concentration2, concentration3)
  
  if (!continuous) {
    gradient_colors <- grDevices::colorRampPalette(c("#c6dbef", "#08306b"))(length(unique_val))
    names(gradient_colors) <- unique_val
    dt_plate[, (column_name) := factor(get(column_name), levels = unique_val)]
  } else {
    if (is.numeric(dt_plate[[column_name]])) {
      dt_plate[, (column_name) := round(get(column_name), 5)]
    }
  }
  
  
  plate_list <- lapply(unique(dt_plate[[barcode_idf]]), function(x) {
    dt_plate_subset <- dt_plate[get(barcode_idf) == x]
    
    dt_plate_subset[, WellRow := droplevels(WellRow)]
    
    p <- ggplot2::ggplot(dt_plate_subset) +
      ggplot2::geom_rect(ggplot2::aes(xmin = WellColumn - 0.5, xmax = WellColumn + 0.5,
                                      ymin = as.numeric(WellRow) - 0.5,
                                      ymax = as.numeric(WellRow) + 0.5,
                                      fill = !!rlang::sym(column_name)),
                         color = "black", linewidth = 0.2) +
      ggrepel::geom_text_repel(data = dt_plate_subset, 
                               ggplot2::aes(x = WellColumn, y = WellRow, label = !!rlang::sym(column_name)), 
                               color = "black", size = 3, bg.color = "white", bg.r = 0.08, force = 0) +
      ggplot2::scale_x_continuous(breaks = sort(unique(dt_plate_subset$WellColumn)),
                                  labels = sort(unique(dt_plate_subset$WellColumn)),
                                  position = "top") +
      ggplot2::labs(
        title = paste0(length(unique(dt_plate_subset$WellColumn)) * length(unique(dt_plate_subset$WellRow)),
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
      p <- p + ggplot2::scale_fill_gradient(low =  "#c6dbef", high = "#08306b")
    } else {
      p <- p + ggplot2::scale_fill_manual(values = gradient_colors,
                                          name = column_name,
                                          limits = names(gradient_colors))
    }
    
    p + ggplot2::geom_tile(ggplot2::aes(x = WellColumn, y = WellRow),
                           fill = NA, color = "black", linewidth = 0.5)
    
  })
  
  names(plate_list) <- unique(dt_plate[[barcode_idf]])
  plate_list
}


#' Helper function to generate color mappings for Gnumber and Gnumber_2
#' 
#' @param dt_plate_subset The subset of data
#' @param untrt_tag Untreated tag identifier
#' 
#' @keywords internal
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
