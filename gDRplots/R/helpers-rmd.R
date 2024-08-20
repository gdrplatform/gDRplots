#' Prepare markdown chunk based on the plots list
#'
#' Function output should be generated with \code{knitr::knit(text = unlist(<result>))}
#'
#' @param plt_list named list with generated plots to be shown in tabs
#'     (for list without name - only ordinal numbers will be generated)
#' @param chunk_name string name of markdown chunk; preferable without spaces
#' @param header_level numeric level of markdown header
#'
#' @examples
#' plotlist <- lapply(unique(iris$Species), function(iris_name) {
#'   ggplot2::ggplot(iris[iris$Species == iris_name, c("Sepal.Length", "Sepal.Width")]) +
#'   ggplot2::geom_point(ggplot2::aes(x = Sepal.Length, y = Sepal.Width))
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

#' Escape colon and hash
#'
#' @param x String
#'
#' @examples
#' escape_special_characters("ABC:123")
#' escape_special_characters("AD_12")
#' escape_special_characters("AD#12")
#'
#' @return Original string with \code{:}s and \code{#}s escaped
#' @keywords internal
#'
#' @export
escape_special_characters <- function(x) {
  checkmate::assert_string(x)
  if (grepl("\\:", x)) x <- gsub(pattern = "\\:", replacement = "\\\\:", x = x)
  if (grepl("#", x)) x <- gsub(pattern = "#", replacement = "[hash]", x = x)
  x
}

#' Replace spaces with another character
#'
#' @param x String where matches are sought
#' @param replacement String replacement for spaces
#'
#' @examples
#' neutralize_spaces("GDC-123|Abc x G01234")
#' neutralize_spaces("MNO-321P 789R YY#1 ")
#' neutralize_spaces("drug_001 x drug_002", ".")
#' 
#' @return String with spaces replaced by the specified character
#' @keywords internal
#'
#' @export
neutralize_spaces <- function(x,
                              replacement = "_") {
  checkmate::assert_string(x)
  checkmate::assert_string(replacement)
  gsub(" ", replacement, trimws(x))
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
  
  # # Cap the width and height at max_size
  # estimated_width <- min(estimated_width, max_size)
  # estimated_height <- min(estimated_height, max_size)
  
  return(c(width = estimated_width, height = estimated_height))
}

#' Save gDR plots to a specified path
#'
#' @param plt A plot object; either a ggplot2 or pheatmap object.
#' @param path A string specifying the path where the plot should be saved.
#' @param format A string specifying the format for saving the plot; either "svg", "png", or "pdf". Default is "svg".
#'
#' @return NULL
#' @keywords internal
#' @export
#'
#' @examples
#' tmp_dir <- file.path(tempdir(), "plot_dir")
#' dir.create(tmp_dir, showWarnings = FALSE)
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) + ggplot2::geom_point()
#' save_plot(plt = p, path = paste(tmp_dir, "mtcars_scatter", sep = "/"), format = "png")
save_plot <- function(plt, path, format = "svg") {
  checkmate::assert_multi_class(plt, c("ggplot", "pheatmap"))
  checkmate::assert_string(path)
  checkmate::assert_choice(format, choices = c("svg", "png", "pdf"))
  
  # Check if the directory exists and has write access
  dir_path <- dirname(path)
  if (!dir.exists(dir_path)) {
    stop("The specified directory does not exist.")
  }
  
  if (file.access(dir_path, 2) != 0) {
    stop("The specified directory does not have write access.")
  }
  
  # Estimate plot size
  plot_size <- estimate_plot_size(plt)
  
  filename <- paste(path, format, sep = ".")
  
  
  # Save the plot in the specified format
  ggplot2::ggsave(filename = filename,
                  plot = plt,
                  units = "in",
                  width = plot_size[["width"]],
                  height = plot_size[["height"]],
                  dpi = 300,
                  limitsize = FALSE,
                  device = format)
  
  invisible(NULL)
}
