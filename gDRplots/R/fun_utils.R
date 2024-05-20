#' Produce a color palette
#'
#' Build a color palette of any length from RColorBrewer templates.
#'
#' @param n number of colors
#' @param name name of a qualitative RColorBrewer palette
#' @param shuffle whether or not to shuffle the colors in the provided palette
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
#' @keywords internal
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


#' Display color palette
#'
#' Produce an image of a set of colors for easy viewing.
#'
#' @param colors character vector of color names or hex codes
#' 
#' @examples
#' paletteDisplay("red")
#' paletteDisplay("#DD0000")
#'
#' @keywords internal
#' @return `NULL`
#' 
#' @export
paletteDisplay <- function(colors) {
  checkmate::assert_character(colors, any.missing = TRUE) # NA will be shown as white
  stopifnot("Must be valid color name" = 
              all(vapply(stats::na.omit(colors), isValidColor, logical(1))))
  
  n <- length(colors)
  graphics::image(seq_len(n), 1, as.matrix(seq_len(n)), col = colors,
                  xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
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
#' @keywords utils
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
#' @keywords utils
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
#' @keywords utils
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
#' @keywords utils
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


#' Get Legend Title
#' 
#' prepare legend title formatted as HTML code
#'
#' @param var character name of the variable to be shown in legend
#' @param has_codrug_data logical whether the data is combo type
#' @param default_var character name of the variable to be shown in legend when \code{var} is not selected
#' 
#' @examples
#' get_legend_title("Drug")
#' get_legend_title("Drug 2", has_codrug_data = TRUE, default_var = "Conc_2")
#'
#' @keywords utils
#' @return single element list with character in HTML format
#' 
#' @export
get_legend_title <-
  function(var,
           has_codrug_data = FALSE,
           default_var = gDRutils::get_prettified_identifiers("concentration2", simplify = TRUE)) {
    checkmate::assert_string(var, null.ok = TRUE)
    checkmate::assert_flag(has_codrug_data)
    checkmate::assert_string(default_var)
    
    if (is.null(var) || tolower(var) == "none") {
      if (has_codrug_data) {
        list(text = sprintf("<b>%s</b>", default_var))
      } else {
        NULL
      }
    } else {
      if (has_codrug_data && var != default_var) {
        list(text = sprintf("<b>%s</b> && <b>%s</b>",
                            default_var,
                            var))
      } else {
        list(text = sprintf("<b>%s</b>", var))
      }
    }
  }


#' Helper function to evaluate if legend should be added to the figure
#'
#' Current logic is as follows: show legend only if:
#' (1) user colors by "Drug MOA" or "Primary Tissue" or
#' (2) there is a codrug data with at least two unique values
#' 
#' @param var_col name of column used for coloring dots
#' @param data \code{data.table} prepared by \code{prepareDataMC}
#' 
#' @examples
#' dt <- data.table::data.table(
#'   x = seq_len(10),
#'   y = seq_len(10) + 0.5,
#'   `Tissue` = c("lung", "brain")
#' )
#' do_show_legend("none", dt)
#' do_show_legend("Tissue", dt)
#' 
#' @keywords utils
#' @return logical value
#' 
#' @export
do_show_legend <- function(var_col, 
                           data) {
  
  checkmate::assert_string(var_col, null.ok = FALSE)
  checkmate::assert_data_table(data)
  
  pidfs <- gDRutils::get_prettified_identifiers(simplify = TRUE)
  
  var_col %in% c(pidfs[["drug_moa"]], pidfs[["cellline_tissue"]]) ||
    length(unique(data[[pidfs[["concentration2"]]]])) >= 2
}


#' Build labels for plotly tooltips
#'
#' @param data a data table in which the labels are to be constructed
#' @param view which view to prepare label for; character string
#' 
#' @examples
#' SE <- gDRutils::get_synthetic_data("small")[[1]]
#' dt <- convert_se_assay_to_custom_dt(SE, assay_name = "Metrics")
#' buildLabel(dt, "grid")
#'
#' @keywords utils
#' @return A character vector the same length as the row number of \code{data}.
#'
#' @export
buildLabel <- function(data, 
                       view) {
  
  pidfs <- gDRutils::get_prettified_identifiers(simplify = TRUE)
  cell_name <- pidfs[["cellline_name"]]
  drug_name <- pidfs[["drug_name"]]
  concentration_name <- pidfs[["concentration"]]
  drug_moa_name <- pidfs[["drug_moa"]]
  cell_tissue_name <- pidfs[["cellline_tissue"]]
  concentration2_name <- pidfs[["concentration2"]]
  drug2_name <- pidfs[["drug_name2"]]
  
  iso_level <- gDRcomponenst::get_isobologram_columns("iso_level")
  pos_x <- gDRcomponenst::get_isobologram_columns("pos_x")
  pos_x_ref <- gDRcomponenst::get_isobologram_columns("pos_x_ref")
  pos_y <- gDRcomponenst::get_isobologram_columns("pos_y")
  pos_y_ref <- gDRcomponenst::get_isobologram_columns("pos_y_ref")
  log10_ratio <- gDRcomponenst::get_isobologram_columns("log10_ratio")
  log2_CI <- gDRcomponenst::get_isobologram_columns("log2_CI")
  
  checkmate::assert_data_table(data)
  checkmate::assert_string(view)
  checkmate::assert_choice(view,
                           choices = c("clustering", "distribution", "ranking",
                                       "contrast", "curve", "grid",
                                       "combo1-heatmap", "combo1-points",
                                       "combo1-lines_ref", "combo1-lines",
                                       "combo-ratios", "combo3"))
  
  data <- coerce_cotreatment_data(data)
  vars_cotreatment <- intersect(c(drug2_name, concentration2_name), names(data))
  
  if (view  == "grid") {
    ans <- sprintf("%s: %s\n%s: %s", cell_name, data[[cell_name]], drug_name, data[[drug_name]])
    
  } else if (view == "curve") {
    # capture variables from enclosing call
    var_y <- dynGet("var_y", inherits = TRUE)
    var_col <- dynGet("var_col", inherits = TRUE)
    var_not_col <- dynGet("var_not_col", inherits = TRUE)
    if (length(vars_cotreatment) > 0 && !is.null(data[[concentration2_name]])) {
      ans <- sprintf("%s: %s\n%s: %s\n%s: %.4g &mu;M\n(%s at %.4g &mu;M)\n%s: %.2f",
                     var_col, data[["var_col"]],
                     var_not_col, data[["var_not_col"]],
                     concentration_name, data[[concentration_name]],
                     data[[drug2_name]], data[[concentration2_name]],
                     var_y, data[["var_y"]])
    } else {
      ans <- sprintf("%s: %s\n%s: %s\n%s: %.4g &mu;M\n%s: %.2f",
                     var_col, data[["var_col"]],
                     var_not_col, data[["var_not_col"]],
                     concentration_name, data[[concentration_name]],
                     var_y, data[["var_y"]])
    }
    
  } else if (view == "distribution") {
    # capture variables from enclosing call
    var_x <- dynGet("var_x", inherits = TRUE)
    var_y <- dynGet("var_y", inherits = TRUE)
    var_col <- dynGet("var_col", inherits = TRUE)
    title_x <- dynGet("title_x", inherits = TRUE)
    if (length(vars_cotreatment) > 0 && !is.null(data[[concentration2_name]])) {
      ans <- sprintf("%s: %s\n%s: %s\n%s(%s at %.4g &mu;M)\n%s: %.2f",
                     cell_name, data[[cell_name]],
                     drug_name, data[[drug_name]],
                     switch(var_col, "none" = "", sprintf("%s: %s\n", var_col, data[[var_col]])),
                     data[[drug2_name]], data[[concentration2_name]],
                     title_x, data[[var_y]])
    } else {
      ans <- sprintf("%s: %s\n%s: %s\n%s%s: %.2f",
                     cell_name, data[[cell_name]],
                     drug_name, data[[drug_name]],
                     switch(var_col, "none" = "", sprintf("%s: %s\n", var_col, data[[var_col]])),
                     title_x, data[[var_y]])
    }
    
  } else if (view == "contrast") {
    # capture variables from enclosing call
    var_x <- dynGet("var_x", inherits = TRUE)
    var_y <- dynGet("var_y", inherits = TRUE)
    var_txt <- dynGet("var_txt", inherits = TRUE)
    if (length(vars_cotreatment) > 0 && !is.null(data[[concentration2_name]])) {
      ans <- sprintf("%s: %s\n(%s at %.4g &mu;M)\n%s: %.4g\n%s: %.4g",
                     var_txt, data[[var_txt]],
                     data[[drug2_name]], data[[concentration2_name]],
                     var_x, data[[var_x]],
                     var_y, data[[var_y]])
    } else {
      ans <- sprintf("%s: %s\n%s: %.4g\n%s: %.4g",
                     var_txt, data[[var_txt]],
                     var_x, data[[var_x]],
                     var_y, data[[var_y]])
    }
    
  } else if (view == "ranking") {
    # capture variables from enclosing call
    var_x <- dynGet("var_x", inherits = TRUE)
    var_y <- dynGet("var_y", inherits = TRUE)
    var_col <- dynGet("var_col", inherits = TRUE)
    var_grp <- dynGet("var_grp", inherits = TRUE)
    title_x <- dynGet("title_x", inherits = TRUE)
    
    ## in case color and group are unspecified, throw in MOA and tissue information
    if (var_col == var_grp) {
      if (var_col == "none") {
        var_col <- if (var_x == drug_name)  {
          drug_moa_name
        } else if (var_x == cell_name)  {
          cell_tissue_name
        } else {
          stop("bad value provided for 'var_col'")
        }
      }
      var_grp <- if (var_col == drug_moa_name)  {
        cell_tissue_name
      } else if (var_col == cell_tissue_name)  {
        drug_moa_name
      } else {
        stop("bad value provided for 'var_col'")
      }
    }
    ## end filling duplicates
    if (var_col == "none") {
      var_col <- if (var_grp == drug_moa_name)  {
        cell_tissue_name
      } else if (var_grp == cell_tissue_name)  {
        drug_moa_name
      } else {
        stop("bad value provided for 'var_grp'")
      }
    }
    if (var_grp == "none") {
      var_grp <- if (var_col == drug_moa_name)  {
        cell_tissue_name
      } else if (var_col == cell_tissue_name)  {
        drug_moa_name
      } else {
        stop("bad value provided for 'var_grp'")
      }
    }
    ## end filling nones
    
    # build labels proper
    if (length(vars_cotreatment) > 0 && !is.null(data[[concentration2_name]])) {
      ans <- sprintf("%s: %s\n%s: %s\n%s: %s\n(%s at %.4g &mu;M)\n%s: %.2f",
                     var_x, data[[var_x]],
                     var_col, data[[var_col]],
                     var_grp, data[[var_grp]],
                     data[[drug2_name]], data[[concentration2_name]],
                     title_x, data[[var_y]])
    } else {
      ans <- sprintf("%s: %s\n%s: %s\n%s: %s\n%s: %.2f",
                     var_x, data[[var_x]],
                     var_col, data[[var_col]],
                     var_grp, data[[var_grp]],
                     title_x, data[[var_y]])
    }
    
  } else if (view == "combo1-heatmap") {
    # capture variables from enclosing call
    condition <- dynGet("condition", inherits = TRUE)
    matrix_pretty <- dynGet("matrix_pretty", inherits = TRUE)
    ans <- sprintf("Cell Line: %s\n%s: %.2g\n%s: %.2g\n%s: %.2g",
                   condition[cell_name],
                   condition[drug2_name], 10 ^ data[[pos_x]],
                   condition[drug_name], 10 ^ data[[pos_y]],
                   matrix_pretty, data[["value"]])
    
  } else if (view == "combo1-lines_ref") {
    # capture variables from enclosing call
    condition <- dynGet("condition", inherits = TRUE)
    matrix_pretty <- dynGet("matrix_pretty", inherits = TRUE)
    ans <- sprintf("Cell Line: %s\n%s: %.2g\n%s: %.2g\nreference %s",
                   condition[cell_name],
                   condition[drug2_name], 10 ^ data[[pos_x_ref]],
                   condition[drug_name], 10 ^ data[[pos_y_ref]],
                   data[["name"]])
    
  } else if (view == "combo1-lines") {
    # capture variables from enclosing call
    condition <- dynGet("condition", inherits = TRUE)
    matrix_pretty <- dynGet("matrix_pretty", inherits = TRUE)
    ans <- sprintf("Cell Line: %s\n%s: %.2g\n%s: %.2g\nIsobol: %s",
                   condition[cell_name],
                   condition[drug2_name], 10 ^ data[[pos_x]],
                   condition[drug_name], 10 ^ data[[pos_y]],
                   data[["name"]])
    
  } else if (view == "combo1-points") {
    # capture variables from enclosing call
    condition <- dynGet("condition", inherits = TRUE)
    matrix_pretty <- dynGet("matrix_pretty", inherits = TRUE)
    ans <- sprintf("Cell Line: %s\n%s: %.2g\n%s: %.2g\n%s: %.2g",
                   condition[cell_name],
                   condition[drug2_name], 10 ^ data[[pos_x]],
                   condition[drug_name], 10 ^ data[[pos_y]],
                   matrix_pretty, data[["x2_off"]])
    
  } else if (view == "combo-ratios") {
    # capture variables from enclosing call
    condition <- dynGet("condition", inherits = TRUE)
    ans <- sprintf("%s: %s\nlog10_ratio_conc: %.2g\nlog2_CI: %.2g\niso: %s",
                   cell_name,
                   condition[cell_name],
                   data[[log10_ratio]],
                   data[[log2_CI]],
                   data[[iso_level]])
    
  } else if (view == "combo3") {
    ans <- sprintf("level: %.2g\nlog2_CI: %.2g",
                   data[["level"]],
                   data[[log2_CI]])
    
  }
  
  reformat_untreated_cases(ans)
}

