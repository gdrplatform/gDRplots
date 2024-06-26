#' Prepare markdown chunk based on the plots list
#'
#' Function output should be generated with \code{knitr::knit(text = unlist(<result>))}
#'
#' @param plt_list named list with generated plots to be shown in tabs
#' @param chunk_name string name of markdown chunk; preferable without spaces
#' @param header_level numeric level of markdown header
#'
#' @examples
#' plotlist <- lapply(unique(iris$Species), function(iris_name) {
#'   plot(iris[iris$Species == iris_name, c("Sepal.Length", "Sepal.Width")])
#' })
#' names(plotlist) <- unique(iris$Species)
#' 
#' prep_plot_chunk(plotlist, "iris")  
#'
#' @return list of character vector - input for \code{knitr::knit}
#' @keywords internal
#' 
#' @seealso \code{\link[knitr]{knit}}
#' 
#' @export
prep_plot_chunk <- function(plt_list, 
                            chunk_name,
                            header_level = 3) {
  checkmate::assert_list(plt_list)
  checkmate::assert_named(plt_list)
  checkmate::assert_string(chunk_name)
  checkmate::assert_int(header_level, lower = 1)
  
  lvl <- paste0(rep("#", header_level), collapse =  "")
  plt_list_name <- deparse(substitute(plt_list))
  template <- c(
    sprintf("%s `r names(%s)[{{nm}}]`\n", lvl, plt_list_name),
    sprintf("```{r %s {{nm}}, echo = FALSE}\n", chunk_name),
    sprintf("%s[[{{nm}}]] \n", plt_list_name),
    "```\n",
    "\n"
  )
  lapply(seq_along(plt_list), function(nm) {
    knitr::knit_expand(text = template)
  })
}

#' Escape colon
#'
#' @param x String
#'
#' @examples
#' escape_special_characters("ABC:123", "AD_12")
#'
#' @return Original string with \code{:}s escaped
#' @keywords internal
#'
#' @export
escape_special_characters <- function(x) {
  checkmate::assert_string(x)
  gsub(pattern = "\\:", replacement = "\\\\:", x = x)
}
