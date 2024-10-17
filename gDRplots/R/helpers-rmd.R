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
#' \dontrun{
#' plotlist <- lapply(unique(iris$Species), function(iris_name) {
#'   ggplot2::ggplot(iris[iris$Species == iris_name, c("Sepal.Length", "Sepal.Width")]) +
#'   ggplot2::geom_point(ggplot2::aes(x = Sepal.Length, y = Sepal.Width))
#' })
#' names(plotlist) <- unique(iris$Species)
#' 
#' prep_plot_chunk(plotlist, "iris")
#'
#' }
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

#' Prepare markdown chunk based on the nested plots list
#'
#' Function output should be generated with \code{knitr::knit(text = unlist(<result>))}
#'
#' @param plt_list named list with generated plots to be shown in tabs; list of plots in nested 
#'   hierarchy, where last 4th level is plot and 3rd level is \code{normalization_type} describe by
#'   one of: "GR" ("GR Value") or "RV" ("Relative Viability")#'     
#' @param chunk_name string name of markdown chunk; preferable without spaces
#' @param header_level numeric level of markdown header - only for the first level
#'
#' @examples
#' \dontrun{
#' mae <- gDRutils::get_synthetic_data("small")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
#' 
#' # help function
#' plot_col <- function(tab_plt, norm_type, col = "red") {
#'   tab_plt <- data.table::melt(
#'     data = tab_plt[normalization_type == norm_type][, c("rId", "xc50", "x_mean", "x_max")],
#'     id = "rId")
#'   plt <- ggplot2::ggplot(tab_plt, ggplot2::aes(x = variable, y = value)) +
#'     ggplot2::geom_col(fill = col)
#'   return(plt)
#' }
#' 
#' # creating nested list with plots
#' plotlist <- list()
#' ls_color <- c("darkred", "orange", "darkcyan")
#' for (drug in unique(dt_metrics$DrugName)) {
#'   for (cl in unique(dt_metrics$CellLineName)) {
#'     tab_plot <- dt_metrics[DrugName == drug & CellLineName == cl]
#'     
#'     plt_GR <- lapply(ls_color, function(col) plot_col(tab_plot, "RV", col))
#'     names(plt_GR) <- sprintf("%s_%s", "GR", ls_color)
#'     plt_RV <- lapply(ls_color, function(col) plot_col(tab_plot, "RV", col))
#'     names(plt_RV) <- sprintf("%s_%s", "RV", ls_color)
#'     
#'     plotlist[[drug]][[cl]][["RV"]] <- plt_RV
#'     plotlist[[drug]][[cl]][["GR"]] <- plt_GR
#'   }
#' }
#' 
#' prep_nested_plot_chunk(plotlist, "metric_value")
#' 
#' prep_nested_plot_chunk(plotlist, "metric_value")
#' 
#' }
#' @return list of character vector - input for \code{knitr::knit}
#' @keywords internal
#' 
#' @seealso \code{\link[knitr]{knit}}
#' 
#' @export
prep_nested_plot_chunk <- function(plt_list,
                                   chunk_name,
                                   header_level = 2) {
  checkmate::assert_list(plt_list)
  checkmate::assert_named(plt_list)
  checkmate::assert_string(chunk_name)
  checkmate::assert_int(header_level, lower = 1)
  
  lvl_1 <- paste0(rep("#", header_level), collapse =  "")
  lvl_2 <- paste0(rep("#", header_level + 1), collapse =  "")
  lvl_3 <- paste0(rep("#", header_level + 2), collapse =  "")
  lvl_4 <- paste0(rep("#", header_level + 3), collapse =  "")
  plt_list_name <- deparse(substitute(plt_list))
  
  lapply(names(plt_list), function(nm_1) {
    c(
      sprintf("%s %s {.tabset}\n\n", lvl_1, nm_1),
      unlist(
        lapply(names(plt_list[[nm_1]]), function(nm_2) {
          c(
            sprintf("%s %s {.tabset .tabset-fade .tabset-pills}\n\n", lvl_2, nm_2),
            unlist(
              lapply(names(plt_list[[nm_1]][[nm_2]]), function(nm_norm) {
                norm_title <- switch(nm_norm,
                                     "RV" = "Relative Viability",
                                     "GR" = "GR Value")
                c(
                  sprintf("%s %s {.tabset .tabset-dropdown}\n\n", lvl_3, norm_title),
                  unlist(
                    lapply(names(plt_list[[nm_1]][[nm_2]][[nm_norm]]), function(nm_vis) {
                      
                      chunk_name <- sprintf("%s__%s_%s_%s",
                                            chunk_name, nm_1, nm_2, nm_norm)
                      
                      plt_list_name <- sprintf('%s[["%s"]][["%s"]][["%s"]]', 
                                               plt_list_name, nm_1, nm_2, nm_norm)
                      
                      template <- c(
                        sprintf("%s {{nm_vis}} \n\n", lvl_4),
                        sprintf("```{r %s {{nm_vis}}, echo = FALSE}\n", chunk_name),
                        sprintf('%s[["{{nm_vis}}"]] \n', plt_list_name),
                        "```\n",
                        "\n"
                      )
                      knitr::knit_expand(text = template)
                    })
                  )
                )
              })
            )
          )
        })
      )
    )
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
#'
#' @return named vector with optimal width and height used in \code{\link[ggplot2]{ggsave}} function
#' @keywords internal
#' 
#' @examples
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) + ggplot2::geom_point()
#' estimate_plot_size(p)
#' 
#' @export
estimate_plot_size <- function(plt,
                               base_width = 10,
                               base_height = 6,
                               scale_factor = 0.5) {
  
  checkmate::assert_multi_class(plt, c("ggplot", "pheatmap"))
  checkmate::assert_numeric(base_width, lower = 0, finite = TRUE)
  checkmate::assert_numeric(base_height, lower = 0, finite = TRUE)
  checkmate::assert_numeric(scale_factor, lower = 0, finite = TRUE)

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
  return(c(width = estimated_width, height = estimated_height))
}

#' Save gDR plots to a specified path
#'
#' @param plt A plot object; either a ggplot2 or pheatmap object.
#' @param path A string specifying the path where the plot should be saved.
#' @param format A string specifying the format for saving the plot; either "svg", "png", or "pdf". Default is "svg".
#'
#' @return \code{NULL}
#' @keywords internal
#' 
#' @seealso \code{\link[ggplot2]{ggsave}}
#'
#' @examples
#' tmp_dir <- file.path(tempdir(), "plot_dir")
#' dir.create(tmp_dir, showWarnings = FALSE)
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) + ggplot2::geom_point()
#' save_plot(plt = p, path = paste(tmp_dir, "mtcars_scatter", sep = "/"), format = "png")
#' 
#' @export 
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
