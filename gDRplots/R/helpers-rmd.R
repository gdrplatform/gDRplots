#' Prepare markdown chunk based on a list of plots
#'
#' Generates markdown code for displaying plots in a document using `knitr::knit()`.
#' The function handles both simple lists of plots and nested lists, allowing for the creation of tabbed
#' sections for grouped plots.
#'
#' @param plt_list A named list of plots.  Names will be used as headings for plots/tab groups.
#' If unnamed, ordinal numbers will be used.  Can be nested lists for tabbed output.
#' If a nested list is provided, the inner lists should also be named.
#' @param chunk_name A character string specifying the base name for the generated code chunks.  Avoid spaces.
#' @param header_level An integer specifying the markdown header level to use (e.g., 1 for `#`, 2 for `##`, etc.).
#' @param tabset_options A character vector of options for the tabset. This is only used when `plt_list`
#' is a nested list.
#' Possible values are "unnumbered", "tabset", and "tabset-dropdown" or other supported by RMarkdown.
#'
#' @return A list of character vectors. Each element of the list corresponds to a plot or a
#' group of plots (if `plt_list` is nested). Each character vector within the list represents the markdown
#' code to be processed by `knitr::knit()`. For nested lists, each element will contain a "header" element
#' and an "items" element. The "header" element is the markdown for the tabset header, and the "items" element
#' is a list of markdown chunks for each tab.
#'
#' @examples
#' \dontrun{
#' # Simple list of plots
#' plotlist <- lapply(unique(iris$Species), function(iris_name) {
#'   ggplot2::ggplot(iris[iris$Species == iris_name, c("Sepal.Length", "Sepal.Width")]) +
#'   ggplot2::geom_point(ggplot2::aes(x = Sepal.Length, y = Sepal.Width))
#' })
#' names(plotlist) <- unique(iris$Species)
#'
#' prep_plot_chunk(plotlist, "iris")
#'
#' # Nested list of plots for tabbed output
#' nested_plotlist <- list()
#' for (species in unique(iris$Species)) {
#'   nested_plotlist[[species]] <- list()
#'   nested_plotlist[[species]][["Sepal"]] <- ggplot2::ggplot(iris[iris$Species == species, ],
#'     ggplot2::aes(x = Sepal.Length, y = Sepal.Width)) + ggplot2::geom_point()
#'   nested_plotlist[[species]][["Petal"]] <- ggplot2::ggplot(iris[iris$Species == species, ],
#'     ggplot2::aes(x = Petal.Length, y = Petal.Width)) + ggplot2::geom_point()
#' }
#'
#' prep_plot_chunk(nested_plotlist, "iris_nested", tabset_options = c("tabset", "unnumbered"))
#' }
#' @keywords internal
#' @seealso \code{\link[knitr]{knit}}
#' @export
prep_plot_chunk <- function(plt_list,
                            chunk_name,
                            header_level = 3,
                            tabset_options = c("unnumbered", "tabset", "tabset-dropdown")) {

  checkmate::assert_list(plt_list)
  checkmate::assert_string(chunk_name)
  checkmate::assert_int(header_level, lower = 1)
  checkmate::assert_character(tabset_options, null.ok = TRUE)

  plt_list_name <- deparse(substitute(plt_list))
  lvl <- paste0(rep("#", header_level), collapse = "")

  lapply(seq_along(plt_list), function(nm) {
    group_name <- names(plt_list)[nm]

    if (inherits(plt_list[[nm]], "list") && !is.null(names(plt_list[[nm]]))) {
      # nested list - use tabset options
      header <- if (is.null(tabset_options)) {
        sprintf("%s %s\n\n", lvl, group_name)
      } else {
        tabset_string <- paste0("{.", paste(tabset_options, collapse = " ."), "}")
        sprintf("%s %s %s\n\n", lvl, group_name, tabset_string)
      }

      item_chunks <- lapply(names(plt_list[[nm]]), function(item_name) {
        chunk <- sprintf(
          "%s# %s\n```{r %s_%s_%s, echo = FALSE}\n%s[[\"%s\"]][[\"%s\"]] \n```\n\n",
          lvl, item_name, chunk_name, group_name, item_name, plt_list_name, group_name, item_name
        )
        knitr::knit_expand(text = chunk)
      })

      list(header = knitr::knit_expand(text = header), items = item_chunks)

    } else {
      # not nested - no tabset, access element by index
      chunk <- c(
        sprintf("%s %s\n", lvl, group_name),
        sprintf("```{r %s_%s, echo = FALSE}\n", chunk_name, group_name),
        sprintf("%s[[%d]] \n", plt_list_name, nm),  # Use %d and nm directly
        "```\n",
        "\n"
      )
      knitr::knit_expand(text = chunk)
    }
  })
}


#' Prepare markdown chunk based on the nested plots list
#'
#' Function output should be generated with \code{knitr::knit(text = unlist(<result>))}
#'
#' @param plt_list named list with generated plots to be shown in tabs; list of plots in nested
#'   hierarchy, where last 4th level is plot and 3rd level is \code{normalization_type} described by
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
#' @return list of character vectors - input for \code{knitr::knit}
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

                      chunk <- c(
                        sprintf("%s {{nm_vis}} \n\n", lvl_4),
                        sprintf("```{r %s {{nm_vis}}, echo = FALSE}\n", chunk_name),
                        sprintf('%s[["{{nm_vis}}"]] \n', plt_list_name),
                        "```\n",
                        "\n"
                      )
                      knitr::knit_expand(text = chunk)
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
#' escape_special_characters("AD/12")
#'
#' @return Original string with \code{:}s and \code{#}s and \code{/}s escaped
#' @keywords internal
#'
#' @export
escape_special_characters <- function(x) {
  checkmate::assert_string(x)
  if (grepl("\\:", x)) x <- gsub(pattern = "\\:", replacement = "[colon]", x = x)
  if (grepl("\\/", x)) x <- gsub(pattern = "\\/", replacement = "[slash]", x = x)
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

#' Extract path of the executed R file
#'
#' @param test_mode logical flag whether the function be run in the test mode
#' @export
#' @return string with the path to the executed Rscript file
#' @keywords internal
get_r_file_path <-  function(test_mode = FALSE) {
  checkmate::assert_flag(test_mode)

  # on Rstudio
  fpath <- if (.Platform$GUI == "RStudio" && !test_mode) {
    rstudioapi::getActiveDocumentContext()$path
  } else {
    # in terminal/test mode
    ca <- commandArgs()
    fpath <- strsplit(ca[grepl("^--file=", ca)], "=")[[1]][2]
    tools::file_path_as_absolute(fpath)
  }
  checkmate::assert_file_exists(fpath)
  fpath
}

#' Prepare markdown chunk based on a doubly nested list of tables
#'
#' Generates markdown code for displaying tables in a document using `knitr::knit()`.
#' Handles doubly nested lists, allowing for tabbed sections for cell lines and then metrics.
#' The inner header level (for metrics) is automatically set one level greater than the outer header level.
#'
#' @param tbl_list A doubly nested named list of tables. The outer list represents cell lines,
#'   and the inner lists represent metrics.  Names are used as headings.
#' @param chunk_name Base name for generated code chunks. Avoid spaces.
#' @param header_level Markdown header level for the outer tabset (cell lines).
#' @param tabset_options Options for the tabset. Can be "unnumbered", "tabset", "tabset-dropdown".
#'
#' @return A list of character vectors. Each element corresponds to a cell line. Each character vector
#'   represents markdown code for the cell line's tabset.
#'
#' @examples
#' \dontrun{
#' nested_tables <- list(
#'   CellLine1 = list(MetricA = mtcars[1:5, ], MetricB = mtcars[6:10, ]),
#'   CellLine2 = list(MetricC = iris[1:5, ], MetricD = iris[6:10, ])
#' )
#' prep_double_table_chunk(nested_tables, "nested_tables", header_level = 2, tabset_options = "tabset")
#' 
#' # Example using DT::datatable
#' prep_double_table_chunk(nested_tables, "dt_tables", header_level = 2, tabset_options = "tabset")
#' }
#' @export
prep_double_table_chunk <- function(tbl_list,
                                    chunk_name,
                                    header_level = 3,
                                    tabset_options = c("unnumbered", "tabset", "tabset-dropdown")) {
  
  checkmate::assert_list(tbl_list, min.len = 1)
  checkmate::assert_string(chunk_name)
  checkmate::assert_int(header_level, lower = 1)
  checkmate::assert_character(tabset_options, null.ok = TRUE)
  
  tbl_list_name <- deparse(substitute(tbl_list))
  lvl <- paste0(rep("#", header_level), collapse = "")
  inner_lvl <- paste0(rep("#", header_level), collapse = "") # Inner level is one greater
  
  lapply(names(tbl_list), function(cell_line) {
    tabset_string <- if (is.null(tabset_options)) {
      ""
    } else {
      paste0("{.", paste(tabset_options, collapse = " ."), "}")
    }
    
    header <- sprintf("%s %s %s\n\n", lvl, cell_line, tabset_string)
    
    item_chunks <- lapply(names(tbl_list[[cell_line]]), function(metric) {
      chunk <- sprintf(
        "%s# %s\n```{r %s_%s_%s, echo = FALSE}\n%s \n```\n\n",
        inner_lvl, 
        metric, 
        chunk_name, 
        cell_line, 
        metric, 
        paste0(
          "DT::formatRound(",
          "DT::datatable(", 
          tbl_list_name, 
          "[[\"", cell_line, "\"]][[\"", metric, "\"]]), ",
          "columns = names(Filter(is.numeric, ", 
          tbl_list_name, 
          "[[\"", cell_line, "\"]][[\"", metric, "\"]])), ", 
          "digits = 5)"
        )
      )
      knitr::knit_expand(text = chunk)
    })
    
    list(header = knitr::knit_expand(text = header), items = item_chunks)
  })
}
