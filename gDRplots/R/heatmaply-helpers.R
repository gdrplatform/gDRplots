#' Force setting limits for colorbar in heatmaply
#'
#' @param heatmaply object of class `heatmapr`
#' @param limits a two dimensional numeric vector specifying the data range for the scale
#' 
#' @examples
#' heatmaply <- heatmaply::heatmaply(mtcars)
#' heatmaply$x$data_index <- 4
#' force_heatmaply_limits(heatmaply, c(100, 400))
#' 
#'
#' @return object of `heatmapr` class with specified data range
#' @keywords heatmaply
#' @export
force_heatmaply_limits <- function(heatmaply, limits) {
  checkmate::assert_class(heatmaply, "plotly")
  checkmate::assert_numeric(limits, len = 2)
  
  range_index <- which(vapply(heatmaply$x$data, function(x) "zmin" %in% names(x),
                              logical(1)))
  
  if (length(range_index) == 0) {
    range_index <- heatmaply$x$data_index
  }
  
  heatmaply$x$data[[range_index]]$zmin <- min(limits)
  heatmaply$x$data[[range_index]]$zmax <- max(limits)
  heatmaply
}

#' Create color palette that takes into account 0 as a mid point with neutral color
#' 
#' Please note that for a completely positive or completely negative range, the neutral color 
#' will be applied to the value closest to zero,
#'
#' @param colors_vec colors to interpolate; must be a valid argument to `col2rgb()`. 
#' Color in the middle of the vector will be used as a neutral color that represent zero value
#' @param limits a two dimensional numeric vector specifying the data range for the scale 
#' to generate vector of colors
#' @param breaks a number of breaks used to divide limits into part for which colors will be generated
#'
#' @return vector of RGB colors
#' 
#' @examples
#' test_colors <- c("blue", "white", "red")
#' create_color_palette(colors_vec = test_colors, limits = c(1, 4), breaks = 1)
#' # "#FFFFFF" "#FF7F7F" "#FF0000"
#' 
#' create_color_palette(colors_vec = test_colors, limits = c(0, 3), breaks = 1)
#' # "#FFFFFF" "#FF7F7F" "#FF0000"
#' 
#' create_color_palette(colors_vec = test_colors, limits = c(-1.5, 2), breaks = 0.5)
#' # "#0000FF" "#5555FF" "#AAAAFF" "#FFFFFF" "#FFAAAA" "#FF5555" "#FF0000"
#' 
#' @keywords heatmaply
#' @export
create_color_palette <- function(colors_vec, 
                                 limits, 
                                 breaks = 0.05) {
  checkmate::assert_character(colors_vec, min.len = 2)
  checkmate::assert_numeric(limits, len = 2)
  checkmate::assert_number(breaks)
  
  mid_color_idx <- ceiling(length(colors_vec) / 2)
  breaks <- abs(breaks)
  
  if (all(limits == 0)) { 
    no_col_below <- 0
    no_col_under <- 1 # return neutral color
  } else if (all(limits <= 0)) { 
    no_col_below <- ceiling(stats::dist(limits)[1] / breaks)
    no_col_under <- 0
  } else if (all(limits >= 0)) { 
    no_col_below <- 0
    no_col_under <- ceiling(stats::dist(limits)[1] / breaks)
  } else { 
    no_col_below <- ceiling(stats::dist(c(min(limits), 0))[1] / breaks) + 1
    no_col_under <- ceiling(stats::dist(c(0, max(limits)))[1] / breaks)
  }
  
  # preparing palette
  colors_below <-
    grDevices::colorRampPalette(colors = colors_vec[seq_len(mid_color_idx)])(no_col_below)
  colors_under <-
    grDevices::colorRampPalette(colors = colors_vec[mid_color_idx:length(colors_vec)])(no_col_under)
  
  unique(c(colors_below, colors_under))
}

