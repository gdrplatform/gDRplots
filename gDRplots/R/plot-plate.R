#' Plot plate data
#' 
#' @param dt_plate The data table containing plate information
#' 
#' @return A list of ggplot objects for each barcode
#' @keywords QC_plot
#' @examples
#' test_data <- data.table::data.table(
#'   WellColumn = rep(1:12, each = 8),
#'   clid = c("A", "A", "B", "B"),
#'   WellRow = rep(LETTERS[1:8], times = 12),
#'   Gnumber = c(rep("untreated", 48), sample(1:5, size = 48, replace = TRUE)),
#'   Gnumber_2 = c(rep("untreated", 48), sample(1:5, size = 48, replace = TRUE)),
#'   Concentration = runif(96, min = 0, max = 100),
#'   ReadoutValue = runif(96, min = 0, max = 100),
#'   Barcode = rep(c("A", "B"), 48)
#'   )
#' plot_plate_stack_info(test_data)[[1]]
#' @export
plot_plate_stack_info <- function(dt_plate) {
  checkmate::assert_data_table(dt_plate)
  checkmate::assert_names(names(dt_plate),
                          must.include = c("ReadoutValue",
                                           unlist(gDRutils::get_env_identifiers(c("drug",
                                                                                  "concentration",
                                                                                  "cellline",
                                                                                  "well_position"),
                                                                                simplify = FALSE))))
  
  barcode_idf <- intersect(gDRutils::get_env_identifiers("barcode"), names(dt_plate))
  if (NROW(barcode_idf) == 0) {
    ggplot2::ggplot()
  } else {
    barcodes <- stringr::str_sort(unique(dt_plate[[barcode_idf]]), numeric = TRUE)
    
    plate_list <- lapply(barcodes, function(x) {
      dt_plate_subset <- filter_data_by_barcode(dt_plate, x, barcode_idf)
      
      
      conc_cols <- unlist(gDRutils::get_env_identifiers(c("concentration", "concentration2"),
                                                        simplify = FALSE))
      drug_cols <- unlist(gDRutils::get_env_identifiers(c("drug", "drug2"),
                                                        simplify = FALSE))
      cl <- gDRutils::get_env_identifiers("cellline")
      
      has_combo <- conc_cols[["concentration2"]] %in% colnames(dt_plate_subset) &&
        drug_cols[["drug2"]] %in% colnames(dt_plate_subset)
      color_mapping <- generate_color_mappings(dt_plate_subset)
      
      doses <- unique(unlist(dt_plate_subset[, intersect(names(dt_plate_subset),
                                                         conc_cols), with = FALSE]))
      doses <- sort(doses)
      
      gradient_colors <- grDevices::colorRampPalette(c("#c6dbef", "white", "#08306b"))(length(doses))
      colors <- gradient_colors
      names(colors) <- doses
      
      dt_plate_subset[[conc_cols[["concentration"]]]] <- factor(dt_plate_subset[[conc_cols[["concentration"]]]],
                                                                levels = doses)
      if (has_combo) {
        dt_plate_subset[[conc_cols[["concentration2"]]]] <- factor(dt_plate_subset[[conc_cols[["concentration2"]]]],
                                                                   levels = doses)
      }
      
      if (has_combo) {
        p <- ggplot2::ggplot() +
          ggplot2::geom_rect(data = dt_plate_subset, 
                             ggplot2::aes(xmin = WellColumn - 0.5,
                                          xmax = WellColumn,
                                          ymin = as.numeric(WellRow) - 0.5,
                                          ymax = as.numeric(WellRow) + 0.5,
                                          fill = !!rlang::sym(conc_cols[["concentration"]])), 
                             color = "black", size = 0.2) + 
          ggplot2::geom_rect(data = dt_plate_subset,
                             ggplot2::aes(xmin = WellColumn,
                                          xmax = WellColumn + 0.5,
                                          ymin = as.numeric(WellRow) - 0.5,
                                          ymax = as.numeric(WellRow) + 0.5,
                                          fill = !!rlang::sym(conc_cols[["concentration2"]])), 
                             color = "black", size = 0.2) +  
          ggplot2::geom_point(data = dt_plate_subset, 
                              ggplot2::aes(x = WellColumn - 0.25,
                                           y = WellRow,
                                           color = !!rlang::sym(drug_cols[["drug"]]),
                                           shape = !!rlang::sym(cl)),
                              size = 3) +
          ggplot2::geom_point(data = dt_plate_subset, 
                              ggplot2::aes(x = WellColumn + 0.25,
                                           y = WellRow,
                                           color = !!rlang::sym(drug_cols[["drug2"]]),
                                           shape = !!rlang::sym(cl)),
                              size = 3)
        
      } else {
        p <- ggplot2::ggplot() +
          ggplot2::geom_rect(data = dt_plate_subset, 
                             ggplot2::aes(xmin = WellColumn - 0.5,
                                          xmax = WellColumn + 0.5,
                                          ymin = as.numeric(WellRow) - 0.5,
                                          ymax = as.numeric(WellRow) + 0.5,
                                          fill = !!rlang::sym(conc_cols[["concentration"]])), 
                             color = "black", size = 0.2) + 
          ggplot2::geom_point(data = dt_plate_subset, 
                              ggplot2::aes(x = WellColumn,
                                           y = WellRow,
                                           color = !!rlang::sym(drug_cols[["drug"]]),
                                           shape = !!rlang::sym(cl)),
                              size = 3)
      }
      
      p <- p +
        ggrepel::geom_text_repel(data = dt_plate_subset, 
                                 ggplot2::aes(x = WellColumn, y = WellRow, label = round(ReadoutValue, 1)), 
                                 color = "black", size = 3, bg.color = "white", bg.r = 0.05, force = 0, nudge_y = -0.3, segment.color = NA) +
        ggplot2::scale_fill_manual(values = colors, name = conc_cols[["concentration"]], limits = names(colors)) +
        ggplot2::scale_x_continuous(breaks = sort(unique(dt_plate_subset$WellColumn)),
                                    labels = sort(unique(dt_plate_subset$WellColumn)),
                                    position = "top") +
        ggplot2::labs(
          title = paste0(length(unique(dt_plate_subset$WellColumn)) * length(unique(dt_plate_subset$WellRow)), "-Well Visualization of Plate ", x),
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
        ggplot2::scale_color_manual(values = color_mapping) +
        ggplot2::scale_shape_manual(values = scales::shape_pal()(length(unique(dt_plate_subset[[cl]]))))  # Automatically generate distinct shapes
      
      p + ggplot2::geom_tile(data = dt_plate_subset,
                             ggplot2::aes(x = WellColumn,
                                          y = WellRow),
                             fill = NA,
                             color = "black",
                             size = 0.5)
    })
    names(plate_list) <- barcodes
    plate_list
  }
}

#' Plot data from a specific column
#' 
#' @param dt_plate The data table containing plate information
#' @param column_name The name of the column to plot
#' 
#' @return A list of ggplot objects for each barcode and column
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
#' plot_plate(test_data, "Gnumber")[[1]]
#' @export
plot_plate <- function(dt_plate, column_name) {
  checkmate::assert_data_table(dt_plate)
  checkmate::assert_names(names(dt_plate),
                          must.include = c("WellColumn",
                                           "WellRow",
                                           "ReadoutValue",
                                           unlist(gDRutils::get_env_identifiers(c("drug",
                                                                                  "concentration",
                                                                                  "cellline",
                                                                                  "well_position"),
                                                                                simplify = FALSE))))
  checkmate::assert_choice(column_name, choices = names(dt_plate))
  
  conc_cols <- unlist(gDRutils::get_env_identifiers(c("concentration", "concentration2", "concentration3"),
                                                    simplify = FALSE))
  
  barcode_idf <- intersect(gDRutils::get_env_identifiers("barcode"), names(dt_plate))
  if (NROW(barcode_idf) == 0) {
    ggplot2::ggplot()
  } else {
    barcodes <- stringr::str_sort(unique(dt_plate[[barcode_idf]]), numeric = TRUE)
    
    plate_list <- lapply(barcodes, function(x) {
      dt_plate_subset <- filter_data_by_barcode(dt_plate, x, barcode_idf)
      
      if (is.numeric(dt_plate_subset[[column_name]])) {
        dt_plate_subset[[column_name]] <- round(dt_plate_subset[[column_name]], 5)
      }
      
      unique_val <- unique(dt_plate_subset[[column_name]])
      
      continuous <- is.numeric(unique_val) && 
        !column_name %in% conc_cols
      
      unique_val <- sort(unique_val)
      
      if (!continuous) {
        gradient_colors <- grDevices::colorRampPalette(c("#c6dbef", "#08306b"))(length(unique_val))
        colors <- gradient_colors
        names(colors) <- unique_val
        
        dt_plate_subset[[column_name]] <- factor(dt_plate_subset[[column_name]], levels = unique_val)
      }
      
      p <- ggplot2::ggplot() +
        ggplot2::geom_rect(data = dt_plate_subset, 
                           ggplot2::aes(xmin = WellColumn - 0.5, xmax = WellColumn + 0.5, ymin = as.numeric(WellRow) - 0.5, ymax = as.numeric(WellRow) + 0.5, fill = .data[[column_name]]), 
                           color = "black", size = 0.2)
      
      p <- p +
        ggrepel::geom_text_repel(data = dt_plate_subset, 
                                 ggplot2::aes(x = WellColumn, y = WellRow, label = dt_plate_subset[[column_name]]), 
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
        p <- p + ggplot2::scale_fill_manual(values = colors, name = column_name, limits = names(colors))
      }
      
      p + ggplot2::geom_tile(data = dt_plate_subset,
                             ggplot2::aes(x = WellColumn,
                                          y = WellRow),
                             fill = NA,
                             color = "black",
                             size = 0.5)
    })
    names(plate_list) <- barcodes
    plate_list
  }
}

#' Helper function to filter data by barcode and prepare it for plotting
#' 
#' @param dt_plate The original data table
#' @param barcode The barcode to filter the data by
#' @param barcode_idf The identifier for the barcode column
#' 
#' @return A subset of the data filtered by the specified barcode
filter_data_by_barcode <- function(dt_plate, barcode, barcode_idf) {
  checkmate::assert_data_table(dt_plate)
  checkmate::assert_character(barcode)
  checkmate::assert_string(barcode_idf)
  checkmate::assert_names(names(dt_plate), must.include = barcode_idf)
  
  dt_plate_subset <- data.table::copy(dt_plate)
  dt_plate_subset <- dt_plate_subset[get(barcode_idf) == barcode]
  
  dt_plate_subset$WellColumn <- as.numeric(dt_plate_subset$WellColumn)
  dt_plate_subset$WellRow <- factor(dt_plate_subset$WellRow,
                                    sort(unique(dt_plate_subset$WellRow), decreasing = TRUE))
  return(dt_plate_subset)
}

#' Helper function to generate color mappings for Gnumber and Gnumber_2
#' 
#' @param dt_plate_subset The subset of data
#' @param untrt_tag Untreated tag identifier
#' 
#' @return A list containing color mappings for Gnumber and Gnumber_2
generate_color_mappings <- function(dt_plate_subset,
                                    untrt_tag = gDRutils::get_env_identifiers("untreated_tag")) {
  checkmate::assert_data_table(dt_plate_subset)
  checkmate::assert_character(untrt_tag)
  special_colors <- rep("darkgray", length(untrt_tag))
  names(special_colors) <- untrt_tag
  
  
  gnumber_cols <- unlist(gDRutils::get_env_identifiers(c("drug", "drug2", "drug3"), simplify = FALSE))
  other_gnumbers <- setdiff(unique(unlist(dt_plate_subset[, intersect(names(dt_plate_subset),
                                                                      gnumber_cols), with = FALSE])),
                            untrt_tag)
  palette_colors <- scales::hue_pal()(length(other_gnumbers))
  names(palette_colors) <- other_gnumbers
  color_mapping <- c(special_colors, palette_colors)
  
  return(color_mapping)
}
