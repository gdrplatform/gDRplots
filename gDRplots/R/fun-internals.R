#' internal function
#'
#' @param data data.table potentially containing cotreatment data
#' 
#' @examples
#' # lets prepare some dummy date which will mimic combination data
#' dt <- data.table::data.table(x = c("1","0.1"), y = c("a", "b"))
#' names(dt) <- 
#'    gDRutils::get_prettified_identifiers(c("concentration2", "drug_name2"), simplify = FALSE)
#' dt
#' coerce_cotreatment_data(dt)
#' dt
#'
#' @return side effect of coercing any relevant cotreatment variables to the appropriate 'type'
#' @export
#' @keywords utils
#'
#' @seealso \code{MetricClustering}, \code{buildLabelClustering}
#'
coerce_cotreatment_data <- function(data) {
  checkmate::assert_data_table(data)
  
  pidfs <- gDRutils::get_prettified_identifiers(simplify = TRUE)
  
  vars_cotreatment <- intersect(c(pidfs[["drug_name2"]], pidfs[["concentration2"]]), names(data))
  if (length(vars_cotreatment) > 0L) {
    # column Concentration_2 should be numeric and not character/factor
    if (is.character(data[[pidfs[["concentration2"]]]]) ||
        is.factor(data[[pidfs[["concentration2"]]]])) {
      data[[pidfs[["concentration2"]]]] <- as.numeric(as.character(data[[pidfs[["concentration2"]]]]))
    }
  }
  data
}


#' cleanup untreated tags
#'
#' @param label character/numeric vector with untreated tags to correct
#' 
#' @examples
#' reformat_untreated_cases("Cell Line Name: d\nDrug Name: a\n(untreated at 0 &mu;M)\nTitle: 1.00")
#' # Cell Line Name: d\nDrug Name: a\nTitle: 1.00
#' reformat_untreated_cases(c(
#'   "Cell Line Name: d\nDrug Name: a\n(untreated at 0 &mu;M)\nTitle: 1.00", 
#'   "Cell Line Name: d\nDrug Name: a\n(untreated at 4 &mu;M)\nTitle: 1.00"
#'  ))
#'  # [1] Cell Line Name: d\nDrug Name: a\nTitle: 1.00
#'  # [2] Cell Line Name: d\nDrug Name: a\n(untreated at 4 &mu;M)\nTitle: 1.00
#'
#' @return character/numeric vector corrected untreated tags
#' @export
#' @keywords utils
#'
#' @seealso \code{MetricClustering}, \code{buildLabelClustering}
#' 
reformat_untreated_cases <- function(label) {
  
  checkmate::assert_multi_class(label, c("character", "numeric"))
  sub("\\\n\\(.*? at 0\\.?0* &mu;M\\)", "", label)
}


#' Replace one value with one or more other values in a column of a data table.
#'
#' @param x a \code{data.table}
#' @param column character string; name of column to do the replacement in
#' @param replacee character string; value to be replaced
#' @param replacement character vector of values to replace \code{replacee} with; see \code{Details}
#' 
#' @examples
#' vals <- data.table::data.table(value = runif(10, 1, 2),
#' group = rep(c("A", "B"), each = 10))
#' nts <- data.table::data.table(
#'   value = runif(2, 0, 0.1),
#'   group = rep("nt", 2)
#' )
#' dt <- rbind(vals, nts)
#' dtr <- replaceValues(x = dt, column = "group", replacee = "nt", replacement = c("A", "B"))
#'                               
#'
#' @return A modified data table Keep in mind that if \code{replacement} has length of more than one,
#'         the number of rows will be greater than that of \code{x}.
#'
#' The part of \code{x} where \code{column} is equal to \code{replacee} is removed.
#' A copy of that part of the data is then appended where \code{replacee} has been substituted
#' with a single item of \code{replacement} and this repeated for each item in \code{replacement}.
#' The resulting data frame will be longer than \code{x} if \code{replacement} is longer than 1.
#'
#' @export
#' @keywords utils
#'
#' @seealso \code{MetricRanking}, \code{plotlyMR}
#'
replaceValues <- function(x, column, replacee, replacement) {
  checkmate::assert_data_table(x)
  checkmate::assert_string(column)
  checkmate::assert_choice(column, choices = names(x))
  checkmate::assert_string(replacee)
  checkmate::assert_choice(replacee, unique(as.character(x[[column]])))
  checkmate::assert_character(replacement)
  
  # heavy lifting is done by this function
  replacer <- function(x, column, replacee, replacement) {
    x[[column]] <- ifelse(x[[column]] == replacee, replacement, x[[column]])
    return(x)
  }
  
  # isolate part of x where replacement occurs
  data_replace <- x[x[[column]] == replacee, ]
  # isolate part of x that remains unchanged
  data_remain <- x[x[[column]] != replacee, ]
  # perform replacement
  data_replaced_list <- lapply(replacement, replacer, x = data_replace, column = column, replacee = replacee)
  data_replaced <- do.call(rbind, data_replaced_list)
  
  ans <- rbind(data_remain, data_replaced)
  # drop replaced factor level
  if (is.factor(ans[[column]])) {
    ans[[column]] <- factor(ans[[column]],
                            levels = levels(ans[[column]])[-which(levels(ans[[column]]) == replacee)])
  }
  return(ans)
}


#' Compute distance between rows of a matrix.
#'
#' This function was lifted from the package \code{factoextra} and slightly adjusted.
#' The default method was changed to "spearman" and an option was added to replace
#' missing values in the resulting distance matrix with an arbitrary value.
#' This option defaults to 0 so that the function can be called by
#' \code{iheatmapr::add_row_clustering} and code{iheatmapr::add_col_clustering}.
#'
#' @section Control:
#' The function offers different distance measures but since it is called internally,
#' it will always use default arguments. As a result, should you want to do
#' anything other than compute Spearman correlation distance and replace NAs with 0,
#' you must edit this function definition accordingly. At the moment (December 2020)
#' this is sufficient. If the control were to be given to the user,
#' remove defaults for this definition and create a wrapper within \code{plotlyMH}.
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
#' computeDistances(x)
#'
#' @return Object of class \code{dist}. NA and NaN are substituted by the value of \code{dummy}.
#'
#' @keywords internal
#'
#' @seealso [factoextra::get_dist()]
#'
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#'
#' @export
#'
computeDistances <- function(x, 
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
#' @keywords utils
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
