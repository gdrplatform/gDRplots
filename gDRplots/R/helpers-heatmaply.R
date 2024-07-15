#' Compute distance between rows of a matrix.
#'
#' This function was lifted from the package \code{factoextra} and slightly adjusted.
#' The default method was changed to "spearman" and an option was added to replace
#' missing values in the resulting distance matrix with an arbitrary value.
#' This option defaults to 0 so that the function can be called by
#' \code{iheatmapr::add_row_clustering} and \code{iheatmapr::add_col_clustering}.
#'
#' @section Control:
#' The function offers different distance measures but since it is called internally,
#' it will always use default arguments. As a result, should you want to do
#' anything other than compute Spearman correlation distance and replace NAs with 0,
#' you must edit this function definition accordingly. At the moment (December 2020)
#' this is sufficient. If the control were to be given to the user,
#' remove defaults for this definition and create a wrapper within \code{plotly_metric_clustering}.
#' 
#' The function was originally named \code{computeDistances}
#'
#' @param x numeric matrix
#' @param method character string specifying the distance measure to be used;
#'               must be on of:
#'               "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski",
#'               "pearson", "spearman" or "kendall"
#' @param use character string specifying the method for computing covariances
#'                in the presence of missing values. Only used when
#'                "pearson", "spearman" or "kendall" chosen as a distance measure
#' @param stand logical flag specifying whether the data should be standardized;
#'              if TRUE, columns are converted to z scores using \code{scale}
#' @param dummy value to substitute for missing values; defaults to 0
#' @param ... arguments passed to internally called functions
#' 
#' @examples
#' x <- matrix(1:9, nrow = 3, ncol = 3)
#' rownames(x) <- letters[seq(nrow(x))]
#' compute_distances(x)
#'
#' @return Object of class \code{dist}. NA and NaN are substituted by the value of \code{dummy}.
#'
#' @keywords heatmaply
#'
#' @seealso [factoextra::get_dist()]
#'
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#'
#' @export
#'
compute_distances <- function(x, 
                              method = "spearman", 
                              use = "pairwise.complete.obs", 
                              stand = FALSE, 
                              dummy = 0, ...) {
  
  # row names are required to properly handle x matrices with non-variance rows
  checkmate::assert_matrix(x, row.names = "unique")
  checkmate::assert_numeric(x)
  checkmate::assert_string(method)
  distance_methods <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski",
                        "pearson", "spearman", "kendall")
  checkmate::assert_choice(method, distance_methods)
  checkmate::assert_string(use)
  use_v <- c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs")
  checkmate::assert_choice(use, use_v)
  checkmate::assert_flag(stand)
  checkmate::assert_number(dummy, na.ok = TRUE)
  
  if (stand) x <- scale(x)
  if (method %in% c("pearson", "spearman", "kendall")) {
    
    # calculate correlation only for rows with non-zero variance
    # `rowsMaxs == rowMins` works for Inf and -Inf (opposite to `rowsSds > 0`)
    no_variance_rows_v <- which(matrixStats::rowMaxs(x) == matrixStats::rowMins(x))
    valid_rows_v <- setdiff(seq_len(nrow(x)), no_variance_rows_v)
    
    res_cor <- if (length(no_variance_rows_v) > 0) {
      
      # initialize correlation matrix to zeros
      mat_v <- rep(0, nrow(x) * nrow(x))
      mat_zero <-
        matrix(
          data = mat_v,
          nrow = nrow(x),
          ncol = nrow(x),
          dimnames = list(rownames(x), rownames(x))
        )
      
      # if there are at least two rows with variance
      if (length(valid_rows_v) > 1) {
        
        y <- x[valid_rows_v, , drop = FALSE]
        # protect from the flood of potential warnings
        res_cor <- purrr::quietly(stats::cor)(t(y), method = method, use = use)$result
        
        # matrices have same columns and rows
        common_names <-
          intersect(rownames(mat_zero), rownames(res_cor))
        mat_zero[common_names, common_names] <- res_cor
      }
      mat_zero
    } else {
      # protect from the flood of potential warnings
      res_cor <-
        purrr::quietly(stats::cor)(t(x), method = method, use = use)$result
      res_cor
    }
    res_dist <- stats::as.dist(1 - res_cor, ...)
  } else {
    res_dist <- stats::dist(x, method = method, ...)
  }
  if (!is.na(dummy)) {
    res_dist[is.na(res_dist)] <- dummy
  }
  return(res_dist)
}


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

  # prevent issue with the same vale of both limits
  if (length(unique(limits)) == 1) {
    limits[2] <- limits[1] + breaks
  }

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

#' Get visualization range for each metrics
#' 
#' @examples
#' get_visualization_range()
#' 
#' @keywords heatmaply
#' @return list of numeric vectors with visualization range
#' @export
get_visualization_range <- function() {
  gDRutils::get_settings_from_json("VIS_RANGE",
                                   system.file(package = "gDRplots", "settings.json"))
}
