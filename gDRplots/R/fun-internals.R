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
#' The resulting data.table will be longer than \code{x} if \code{replacement} is longer than 1.
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
