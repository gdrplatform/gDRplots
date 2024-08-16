#' Determine whether or not a color name is valid
#' 
#' A name of color is valid when either is a color name listed by \code{\link[grDevices]{colors}}
#' or a hexadecimal string of the form \code{#rrggbb} 
#' 
#' @param col_name string representing a valid color
#' 
#' @keywords utils_color
#' @return logical flag
#' 
#' @examples
#' is_valid_color("darkblue")
#' is_valid_color("#FF8C00")
#' is_valid_color("#FF8C00DC")
#' is_valid_color("RED")
#' 
#' @export
is_valid_color <- function(col_name) {
  checkmate::assert_string(col_name)

  if (grepl("#", col_name)) {
    grepl("^#?([A-Fa-f0-9]{6}|[A-Fa-f0-9]{8})$", col_name)
  } else {
    col_name %in% grDevices::colors()
  }
}

#' get_iso_colors
#' 
#' 
#' @param  normalization_type charvec normalization_types expected in the data
#' @keywords internal
#'
#' @return named charvec with iso colors
#' 
#' @examples 
#' get_iso_colors()
#' 
#' @export
get_iso_colors <-
  function(normalization_type = c("RV", "GR")) {
    normalization_type <- match.arg(normalization_type)
    iso_cutoff <- seq(0, 1, 0.05)
    breaks <- iso_cutoff
    if (normalization_type == "GR") {
      colors <- vapply(iso_cutoff, function(x) {
        color_vector <- c(70, round((0.85 - x * 0.7) * 170), round((1.1 - x * 0.7) * 200))
        assert_RGB_format(color_vector)
        sprintf("#%s", paste(as.hexmode(color_vector), collapse = ""))
      },
      character(1)
      )
    } else {              
      colors <- vapply(iso_cutoff, function(x) {
        color_vector <- c(70, round((1 - x * .85) * 170), round((1.1 - x * .85) * 232))
        assert_RGB_format(color_vector)
        sprintf("#%s", paste(as.hexmode(color_vector), collapse = ""))
      },
      character(1)
      )
    } 
    names(colors) <- iso_cutoff
    colors
  }

#' Assert whether number may code color in rgb
#'
#' @param x numeric vector describing rgb color 
#' @keywords internal
#'
#' @return NULL
assert_RGB_format <- function(x) {
  if (any(x > 255)) {
    stop("Some value is greater than 255. Not valid RGB format.")
  }
}

