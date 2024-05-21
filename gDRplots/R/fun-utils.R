#' Determine whether or not a colour is dark
#' 
#' @param col_name string representing a valid colour
#' 
#' @examples
#' isColDark("blue")
#' isColDark("red")
#' isColDark("#000000")
#' 
#' @keywords utils_color
#' @return logical flag
#' 
#' @export
isColDark <- function(col_name) {
  checkmate::assert_string(col_name)
  stopifnot("Must be valid color name" = isValidColor(col_name))
  
  getColLuminance(col_name) <= 0.22
}


#' Calculate the luminance of a colour
#'
#' @param col_name string representing a valid colour
#' 
#' @examples
#' getColLuminance("blue")
#' getColLuminance("red")
#' getColLuminance("#000000")
#' getColLuminance("#906090")
#' getColLuminance("#906090F2")
#'
#' @keywords utils_color
#' @return single element numeric vector
#' 
#' @export
getColLuminance <- function(col_name) {
  checkmate::assert_string(col_name)
  stopifnot("Must be valid color name" = isValidColor(col_name))
  
  colrgb <- grDevices::col2rgb(col_name)
  lum <- lapply(colrgb, function(x) {
    x <- x / 255
    if (x <= 0.03928) {
      x <- x / 12.92
    } else {
      x <- ((x + 0.055) / 1.055)^2.4
    }
  })
  lum <- lum[[1]] * 0.2326 + lum[[2]] * 0.6952 + lum[[3]] * 0.0722
  lum
}


#' Determine whether or not a color name is valid
#' 
#' A name of color is valid when either is a color name listed by \code{\link[grDevices]{colors}}
#' or a hexadecimal string of the form \code{#rrggbb} 
#' 
#' @param col_name string representing a valid colour
#' 
#' @keywords utils_color
#' @return logical flag
#' 
#' @examples
#' isValidColor("darkblue")
#' isValidColor("#FF8C00")
#' isValidColor("#FF8C00DC")
#' isValidColor("RED")
#' 
#' @export
isValidColor <- function(col_name) {
  checkmate::assert_string(col_name)
  
  if (grepl("#", col_name)) {
    grepl("^#?([A-Fa-f0-9]{6}|[A-Fa-f0-9]{8})$", col_name)
  } else {
    col_name %in% grDevices::colors()
  }
}

#' Change color name into hexadecimal string
#' 
#' @param col_name string representing a valid color listed by \code{\link[grDevices]{colors}}
#' 
#' @keywords utils_color
#' @return hexadecimal string representing color
#' 
#' @examples
#' colorToHex("darkblue")
#' colorToHex("indianred2")
#' colorToHex("seagreen")
#' 
#' @export
colorToHex <- function(col_name) {
  checkmate::assert_string(col_name)
  stopifnot("Must be valid color name" = col_name %in% grDevices::colors())
  
  rgb <- grDevices::col2rgb(col_name)
  grDevices::rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
}
