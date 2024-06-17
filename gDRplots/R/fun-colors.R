#' Produce a color palette
#'
#' Build a color palette of any length from RColorBrewer templates.
#'
#' @param n numeric representing number of colors
#' @param name string representing name of a qualitative RColorBrewer palette
#' @param shuffle logical representing whether or not to shuffle the colors in the provided palette
#' 
#' @examples
#' paletteBrew(8, "Accent")
#' paletteBrew(20, "Accent")
#' paletteBrew(10, "Set1")
#'
#' @return A character vector of hex codes of length \code{n}.
#'
#' The chosen palette is passed through \code{colorRampPalette}.
#' Unless \code{n} is higher than the maximum length of the given palette,
#' the return value depends on the palette itself and the rules
#' set up in \code{RColorBrewer::brewer.pal}.
#' Otherwise it is extended and includes intermediate colors.
#' Visual properties of the palette are likely to be suffer.
#'
#' @keywords utils_color
#'
#' @seealso \code{RColorBrewer}, \code{colorRampPalette}
#'
#' @export
paletteBrew <- function(n, 
                        name,
                        shuffle = FALSE) {
  
  checkmate::assert_number(n, lower = 1)
  checkmate::assert_string(name)
  checkmate::assert_choice(name, choices = c("Accent", "Dark2", "Paired", "Pastel1",
                                             "Pastel2", "Set1", "Set2", "Set3"))
  checkmate::assert_logical(shuffle)
  
  # brewer.pal throws warnings if minimum/maximum palette lengths are exceeded
  # let's avoid these warnings
  b_pals <- if (n < 3) {
    RColorBrewer::brewer.pal(3, name)[seq(n)]
  } else {
    bpi <-
      data.table::data.table(RColorBrewer::brewer.pal.info, keep.rownames = TRUE)
    n_max <- bpi[bpi$rn == name, "maxcolors"][[1]]
    rep(RColorBrewer::brewer.pal(n_max, name), length.out = n)
  } 
  
  ans <- grDevices::colorRampPalette(b_pals)(n)
  
  if (shuffle) {
    ans <- sample(ans)
  }
  return(ans)
}

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
