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
#' escape_special_characters("ABC:123")
#' escape_special_characters("AD_12")
#'
#' @return Original string with \code{:}s escaped
#' @keywords internal
#'
#' @export
escape_special_characters <- function(x) {
  checkmate::assert_string(x)
  gsub(pattern = "\\:", replacement = "\\\\:", x = x)
}

#' Estimate the optimal plot size (either ggplot or pheatmap) for saving plots
#'
#' @param plt a ggplot or pheatmap object
#' @param base_width an integer with default base_width
#' @param base_height an integer with default base_height
#' @param scale_factor an integer with default scale_factor
#' @param max_size an integer with the maximum size of the plot (either width or height)
#'
#' @return named vector with optimal width and height used in ggsave function
#' @keywords internal
#' @export
estimate_plot_size <- function(plt,
                               base_width = 10,
                               base_height = 6,
                               scale_factor = 0.5,
                               max_size = 49.9) {
  
  checkmate::assert_multi_class(plt, c("ggplot", "pheatmap"))
  checkmate::assert_numeric(base_width, lower = 0, finite = TRUE)
  checkmate::assert_numeric(base_height, lower = 0, finite = TRUE)
  checkmate::assert_numeric(scale_factor, lower = 0, finite = TRUE)
  checkmate::assert_numeric(max_size, lower = 0, finite = TRUE)
  
  
  if (inherits(plt, "ggplot")) {
    # For ggplot2 objects
    plot_data <- ggplot2::ggplot_build(plt)$data
    num_elements <- length(unique(plot_data[[1]]$group))
    estimated_width <- base_width + num_elements * scale_factor
    estimated_height <- base_height + num_elements * scale_factor
  } else if (inherits(plt, "pheatmap")) {
    # For pheatmap objects
    matrix_position <- which(plt$gtable$layout$name == "matrix")
    matrix_dim <- dim(plt$gtable$grobs[[matrix_position]]$children[[1]]$gp$fill)
    num_rows <- matrix_dim[[1]]
    num_cols <- matrix_dim[[2]]
    estimated_width <- base_width + num_cols * scale_factor
    estimated_height <- base_height + num_rows * scale_factor
  } else {
    stop("Unsupported plot type. Only ggplot2 and pheatmap objects are supported.")
  }
  
  # Cap the width and height at max_size
  estimated_width <- min(estimated_width, max_size)
  estimated_height <- min(estimated_height, max_size)
  
  return(c(width = estimated_width, height = estimated_height))
}
