#' Change text alignment in row/column labels in Clustering View.
#'
#' @param x character vector
#' @param pattern character string specifying the places in which to intervene; see \code{Details}
#' 
#' @examples
#' adjustLabel("Cell Line Name: X Drug Name: A (untreated at 0 &mu;M) Title: 1.00")
#' # Cell Line Name: X Drug Name: A &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(untreated at 0 &mu;M) Title: 1.00
#'
#' @return A modified character vector.
#'
#' This is an internal function that inserts whitespace into strings.
#' In cases where a label contains a drug/cell line name with co-treatment information
#' added in parentheses, some non-breaking spaces are added to separate the two.
#' Additionally, a number of non-breaking spaces is added to each string so that
#' they are all of the same length.
#'
#' The whitespace is inserted in places identified by \code{pattern}.
#' This must be a regular expression consisting of two parts (delimited by parentheses)
#' and the non-breaking spaces are inserted between the two.
#'
#' Alignment will be perfect only if a monospaced font is used.
#'
#' @export
#' @keywords utils_label
#'
#' @seealso \code{MetricClustering}, \code{plotlyMH}
#'
adjustLabel <- function(x, pattern = "(.*? )(\\(.*? at .*?\\))") {
  
  checkmate::assert_character(x)
  checkmate::assert_character(pattern, pattern = "^\\(.*?\\)\\(.*?\\)$")
  
  # only act if this pattern is present
  if (!any(grepl(pattern, x))) {
    return(x)
  }
  
  # how long are the strings
  lengths <- nchar(x)
  # how many nbsps to insert
  missings <- max(lengths) - lengths + 5
  # build sequences of nbsps to insert
  fillings <- vapply(missings, function(x) paste(rep("&nbsp;", x), collapse = ""), character(1))
  # insert whitespace
  ans <- mapply(
    function(string, filling) sub(pattern, sprintf("\\1%s\\2", filling), string),
    string = x, filling = fillings)
  
  return(ans)
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
#' @keywords utils_label
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
