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
    valid_rows_v <- setdiff(seq_len(NROW(x)), no_variance_rows_v)
    
    res_cor <- if (length(no_variance_rows_v) > 0) {
      
      # initialize correlation matrix to zeros
      mat_v <- rep(0, NROW(x) * NROW(x))
      mat_zero <-
        matrix(
          data = mat_v,
          nrow = NROW(x),
          ncol = NROW(x),
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

#' Round numbers to unique string
#' 
#' Rounds a numeric vector to the minimum precision needed to ensure uniqueness
#'
#' @param num_vec a numeric vector to be rounded and change into character
#' @param initial_digits numeric value for number of decimal places to start rounding with
#'
#' @keywords utils
#' @returns a character vector of unique numeric strings.
#' 
#' @examples
#' vec <- c(0.00000000, 0.00000256, 0.00001280, 0.00006400, 0.00032000, 
#'          0.00160000, 0.00800000, 0.04000000, 0.20000000, 1.00000000) 
#' 
#' round_to_unique_string(vec)
#' 
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#' 
#' @export
round_to_unique_string <- function(num_vec,
                                   initial_digits = 4) {
  
  checkmate::assert_numeric(num_vec, any.missing = FALSE)
  checkmate::assert_integerish(initial_digits, lower = 1)
  
  initial_rounded_vec <- round(num_vec, digits = initial_digits)
  
  final_char_vec <- format(initial_rounded_vec,
                           scientific = FALSE,
                           trim = TRUE)
  
  if (any(duplicated(final_char_vec))) {
    unique_duplicates <- 
      unique(initial_rounded_vec[duplicated(initial_rounded_vec) | duplicated(initial_rounded_vec, fromLast = TRUE)])
    
    for (d_val in unique_duplicates) {
      i_dup <- which(initial_rounded_vec == d_val)
      
      vec_subset <- num_vec[i_dup]
      max_digits <- max(nchar(format(vec_subset, scientific = FALSE)))
      
      # start the precision search from the next digit
      d <- initial_digits + 1
      while (d <= max_digits) {
        vec_subset_rounded <- round(vec_subset, digits = d)
        
        if (NROW(unique(vec_subset_rounded)) == NROW(vec_subset)) {
          # update the elements in the final_char_vec
          final_char_vec[i_dup] <- format(vec_subset_rounded, 
                                          scientific = FALSE,
                                          trim = TRUE)
          break
        }
        
        d <- d + 1
      }
    }
  }
  
  # final
  return(final_char_vec)
}
