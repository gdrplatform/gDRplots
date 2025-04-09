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
#' @keywords utils
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

#' create a log-sequence
#'
#' Create a sequence of numbers growing in log-domain.
#'
#' The result is a numeric vector of length \code{length}.
#' Differences between items are constant in logarithmic domain
#' and therefore geometrically increase in linear domain.
#'
#' @param start,end numeric, lower and upper margins of the sequence
#' @param length integer, resulting sequence length
#' 
#' @keywords utils
#' @return A numeric vector, see \code{Details}.
#'
#' @examples
#' create_log_seq(1, 2, 5)
#'
#' @export 
create_log_seq <- function(start, end, length) {

  checkmate::assert_number(start, lower = 0, finite = TRUE)
  checkmate::assert_number(end, lower = 0, finite = TRUE)
  checkmate::assert_number(length, lower = 1, finite = TRUE)
  
  limits <- c(start, end)
  limits_log <- log10(limits)
  sequence_log <- seq(from = limits_log[1], limits_log[2], length.out = length)
  sequence <- 10 ^ sequence_log
  # final
  sequence
}                     

