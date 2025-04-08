#' Prepare markdown chunk based on a list of plots
#'
#' Generates markdown code for displaying plots in a document using \code{knitr::knit()}.
#' The function handles both simple lists of plots and nested lists, allowing for the creation of tabbed
#' sections for grouped plots.
#'
#' @param plt_list A named list of plots. Names will be used as headings for plots/tab groups.
#' If unnamed, ordinal numbers will be used.  Can be nested lists for tabbed output.
#' If a nested list is provided, the inner lists should also be named.
#' @param chunk_name A character string specifying the base name for the generated code chunks. Avoid spaces.
#' @param link_list A named list of links to the location (relative paths) where plots are saved, 
#' which when clicked, will be displayed in a new browser tab. It must have the same structure as \code{plt_list}.
#' @param dwn_list A named list of links to location (relative paths) where table or plots are saved, 
#' which when clocked, will be downloaded. It must have the same structure as \code{plt_list}.
#' @param header_level An integer specifying the markdown header level to use (e.g., 1 for `#`, 2 for `##`, etc.).
#' @param tabset_options A character vector of options for the tabset. This is only used 
#' when \code{plt_list} is a nested list.
#' Possible values are "unnumbered", "tabset", and "tabset-dropdown".
#'
#' @return A list of character vectors. Each element of the list corresponds to a plot or a
#' group of plots (if \code{plt_list} is nested). Each character vector within the list represents the markdown
#' code to be processed by \code{knitr::knit()}. For nested lists, each element will contain a "header" element
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
#' 
#' @export
prep_plot_chunk <- function(plt_list,
                            chunk_name,
                            link_list = NULL,
                            dwn_list = NULL,
                            header_level = 3,
                            tabset_options = c("tabset", "tabset-dropdown")) {
  
  checkmate::assert_list(plt_list)
  checkmate::assert_list(link_list, null.ok = TRUE)
  checkmate::assert_list(dwn_list, null.ok = TRUE)
  checkmate::assert_string(chunk_name)
  checkmate::assert_int(header_level, lower = 1)
  checkmate::assert_character(tabset_options, null.ok = TRUE, any.missing = FALSE,
                              pattern = "unnumbered|tabset|tabset-dropdown")
  
  plt_list_name <- deparse(substitute(plt_list))
  lvl <- paste0(rep("#", header_level), collapse = "")
  
  # checking if the structure of plt_list and link_list is identical
  link_structure_condition <- if (inherits(plt_list[[1]], "list")) {
    all(unlist(lapply(seq_along(plt_list), function(i) NROW(plt_list[[i]]))) == 
          unlist(lapply(seq_along(link_list), function(i) NROW(link_list[[i]]))))
  } else {
    NROW(plt_list) == NROW(link_list)
  }
  if (!link_structure_condition) link_list <- NULL
  
  # checking if the structure of plt_list and dwn_list is identical
  dwn_structure_condition <- if (inherits(plt_list[[1]], "list")) {
    all(unlist(lapply(seq_along(plt_list), function(i) NROW(plt_list[[i]]))) == 
          unlist(lapply(seq_along(dwn_list), function(i) NROW(dwn_list[[i]]))))
  } else {
    NROW(plt_list) == NROW(dwn_list)
  }
  if (!dwn_structure_condition) dwn_list <- NULL
  
  lapply(seq_along(plt_list), function(nm) {
    group_name <- ifelse(is.null(names(plt_list)[nm]), nm, names(plt_list)[nm]) # number on name
    
    if (inherits(plt_list[[nm]], "list")) {
      # nested list - use tabset options
      header <- if (is.null(tabset_options)) {
        sprintf("%s %s\n\n", lvl, group_name)
      } else {
        tabset_string <- paste0("{.", paste(tabset_options, collapse = " ."), "}")
        sprintf("%s %s %s\n\n", lvl, group_name, tabset_string)
      }
      
      item_chunks <- lapply(seq_along(plt_list[[nm]]), function(i_nm) {
        item_name <- 
          ifelse(is.null(names(plt_list[[nm]])[i_nm]), i_nm, names(plt_list[[nm]])[i_nm]) # number on name
        
        chunk <- c(
          sprintf("%s# %s\n", lvl, item_name),
          if (!is.null(link_list)) c(create_zoom_link(link_list[[nm]][[i_nm]]), "\n"),
          if (!is.null(dwn_list)) c(create_download_link(dwn_list[[nm]][[i_nm]]), "\n"),
          sprintf("```{r %s_%s_%s, echo = FALSE}\n%s[[%d]][[%d]] \n```\n\n",
                  chunk_name, group_name, item_name, plt_list_name, nm, i_nm)
        )
        knitr::knit_expand(text = chunk)
      })
      
      c(knitr::knit_expand(text = header), unlist(item_chunks))
      
    } else {
      # not nested - no tabset, access element by index
      chunk <- c(
        sprintf("%s %s\n", lvl, group_name),
        if (!is.null(link_list)) c(create_zoom_link(link_list[[nm]]), "\n"),
        if (!is.null(dwn_list)) c(create_download_link(dwn_list[[nm]]), "\n"),
        sprintf("```{r %s_%s, echo = FALSE}\n%s[[%d]] \n```\n\n",
                chunk_name, group_name, plt_list_name, nm)  # Use %d and nm directly
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
#' @inheritParams prep_plot_chunk
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
#' 
#' @keywords internal
#' @seealso \code{\link[knitr]{knit}}
#' 
#' @export
prep_nested_plot_chunk <- function(plt_list,
                                   chunk_name,
                                   link_list = NULL,
                                   header_level = 2) {
  checkmate::assert_list(plt_list)
  checkmate::assert_named(plt_list)
  checkmate::assert_string(chunk_name)
  checkmate::assert_list(link_list, null.ok = TRUE)
  checkmate::assert_int(header_level, lower = 1)
  
  lvl_1 <- paste0(rep("#", header_level), collapse =  "")
  lvl_2 <- paste0(rep("#", header_level + 1), collapse =  "")
  lvl_3 <- paste0(rep("#", header_level + 2), collapse =  "")
  lvl_4 <- paste0(rep("#", header_level + 3), collapse =  "")
  plt_list_name <- deparse(substitute(plt_list))
  
  # checking if the structure of plt_list and link_list is identical
  link_structure_condition <-
    all(
      names(plt_list) %in% names(link_list),
      vapply(names(plt_list), function(nm) {
        all(names(plt_list[[nm]]) == names(link_list[[nm]]),
            vapply(names(plt_list[[nm]]), function(nm_2) {
              all(names(plt_list[[nm]][[nm_2]]) == names(link_list[[nm]][[nm_2]]), 
                  vapply(names(plt_list[[nm]][[nm_2]]), function(nm_3) {
                    all(names(plt_list[[nm]][[nm_2]][[nm_3]]) == names(link_list[[nm]][[nm_2]][[nm_3]]))
                  }, logical(1))
              )
            }, logical(1))
        )
      }, logical(1))
    )
  if (!link_structure_condition) link_list <- NULL
  
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
                      
                      chunk_name <- sprintf("%s__%s_%s_%s_%s",
                                            chunk_name, nm_1, nm_2, nm_norm, nm_vis)
                      
                      plt_list_name_i <- sprintf('%s[["%s"]][["%s"]][["%s"]][["%s"]]',
                                                 plt_list_name, nm_1, nm_2, nm_norm, nm_vis)
                      
                      link_vis <- if (!is.null(link_list)) {
                        c(create_zoom_link(link_list[[nm_1]][[nm_2]][[nm_norm]][[nm_vis]]), "\n")
                      } else {
                        NULL
                      }
                      
                      chunk <- c(sprintf("%s %s \n\n", lvl_4, nm_vis),
                                 link_vis,
                                 sprintf("```{r %s, echo = FALSE}\n%s\n```\n\n",
                                         chunk_name, plt_list_name))
                      
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
#' 
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
#' 
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
#'
#' @examples
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) + ggplot2::geom_point()
#' estimate_plot_size(p)
#'
#' @keywords internal
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
#' 
#' @examples
#' tmp_dir <- file.path(tempdir(), "plot_dir")
#' dir.create(tmp_dir, showWarnings = FALSE)
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) + ggplot2::geom_point()
#' save_plot(plt = p, path = paste(tmp_dir, "mtcars_scatter", sep = "/"), format = "png")
#'
#' @keywords internal
#' @seealso \code{\link[ggplot2]{ggsave}}
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
#' 
#' @return string with the path to the executed Rscript file
#' 
#' @keywords internal
#' 
#' @export
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
#' The inner header level (for metrics) is automatically set to one level greater than the outer header level.
#'
#' @inheritParams prep_plot_chunk
#' @param tbl_list A doubly nested named list of tables. The outer list represents cell lines,
#'   and the inner lists represent metrics.  Names are used as headings.
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
#' }
#' 
#' @keywords internal
#' 
#' @export
prep_double_table_chunk <- function(tbl_list,
                                    chunk_name,
                                    dwn_list = NULL,
                                    header_level = 3,
                                    tabset_options = c("tabset", "tabset-dropdown")) {
  
  checkmate::assert_list(tbl_list, min.len = 1)
  checkmate::assert_string(chunk_name)
  checkmate::assert_list(dwn_list, null.ok = TRUE)
  checkmate::assert_int(header_level, lower = 1)
  checkmate::assert_character(tabset_options, null.ok = TRUE, any.missing = FALSE,
                              pattern = "unnumbered|tabset|tabset-dropdown")
  
  tbl_list_name <- deparse(substitute(tbl_list))
  lvl <- paste0(rep("#", header_level), collapse = "")
  inner_lvl <- paste0(rep("#", header_level), collapse = "") # Inner level is one greater
  
  # checking if the structure of plt_list and dwn_list is identical
  dwn_structure_condition <- all(names(tbl_list) %in% names(dwn_list))
  if (!dwn_structure_condition) dwn_list <- NULL
  
  lapply(names(tbl_list), function(cell_line) {
    tabset_string <- if (is.null(tabset_options)) {
      ""
    } else {
      paste0("{.", paste(tabset_options, collapse = " ."), "}")
    }
    
    header <- c(
      sprintf("%s %s %s\n\n", lvl, cell_line, tabset_string),
      if (!is.null(dwn_list)) c(create_download_link(dwn_list[[cell_line]]), "\n")
    )
    
    item_chunks <- lapply(names(tbl_list[[cell_line]]), function(metric) {
      chunk <- c(
        sprintf("%s# %s\n", inner_lvl, metric),
        sprintf("```{r %s_%s_%s, echo = FALSE}\n%s \n```\n\n",
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
                  "digits = 5)")
        )
      )
      knitr::knit_expand(text = chunk)
    })
    
    c(knitr::knit_expand(text = header), unlist(item_chunks))
  })
}

#' Prepare markdown chunk with zoom link
#' 
#' Generate markdown code for html link item which when clicked opens plot in new browser.
#' The function output should be wrapped in \code{knitr::knit()}.
#'
#' @param img_path string with relative path to file with plot to be shown
#' @param link_txt string with text describing link
#'
#' @return string with html link code
#'
#' @keywords internal
#' @seealso \code{\link[knitr]{knit}}
#' 
#' @export
create_zoom_link <- function(img_path,
                             link_txt = "Zoom In for Details") {
  checkmate::assert_string(img_path, na.ok = TRUE)
  checkmate::assert_string(link_txt)
  
  if (is.na(img_path)) {
    ""
  } else {
    sprintf("<a href=\"%s\" target=\"_blank\">\U1F50D %s</a>",
            img_path, link_txt)
  }
}

#' Prepare markdown chunk with download link
#' 
#' Generate markdown code with a html link item that, when clicked, downloads a file.
#' The function output should be wrapped in \code{knitr::knit()}.
#'
#' @param dwn_path string with relative path to file with plot to be downloaded
#' @param link_txt string with text describing link
#'
#' @return string with html download code
#'
#' @keywords internal
#' @seealso \code{\link[knitr]{knit}}
#' 
#' @export
create_download_link <- function(dwn_path,
                                 link_txt = "Download Table") {
  checkmate::assert_string(dwn_path, na.ok = TRUE)
  checkmate::assert_string(link_txt)
  
  if (is.na(dwn_path)) {
    ""
  } else {
    sprintf("<a href=\"%s\" download>\U0001F4BE %s</a>",
            dwn_path, link_txt)
  }
}


#' Prepare list of names or paths
#' 
#' This function can creates lists that can be used as an input to \code{\link[gDRplots]{prep_plot_chunk}} -
#' param \code{link_list} and \code{dwn_list}. 
#' Depending on the inputs, outputted string will have a format:
#' \emph{<path_file><prefix><name of item on the list>}\strong{.}\emph{<file_format>}
#' Note: An unnamed item in the list will be numbered.
#'
#' @inheritParams prep_plot_chunk
#' @param prefix string to be added as prefix to the name of file
#' @param path_file string path to directory, where data will be read from 
#' @param file_format string specifying the format of file
#'
#' @return list of strings describing names or paths depending on the input; the list is structured 
#' like the input \code{plt_list} and retains the same names as \code{plt_list}.
#'
#' @examples
#' \dontrun{
#' # Simple list of plots
#' plotlist <- lapply(unique(iris$Species), function(iris_name) {
#'   ggplot2::ggplot(iris[iris$Species == iris_name, c("Sepal.Length", "Sepal.Width")]) +
#'     ggplot2::geom_point(ggplot2::aes(x = Sepal.Length, y = Sepal.Width))
#' })
#' 
#' prep_filename_path(plt_list = plotlist,
#'                    prefix = "iris__",
#'                    path_file = file.path(".", "plots"))
#'
#' # Nested list of plots for tabbed output
#' nested_plotlist <- list()
#' for (species in unique(iris$Species)) {
#'   nested_plotlist[[species]] <- list()
#'   nested_plotlist[[species]][["Sepal"]] <- 
#'     ggplot2::ggplot(iris[iris$Species == species, ],
#'                     ggplot2::aes(x = Sepal.Length, y = Sepal.Width)) + ggplot2::geom_point()
#'   nested_plotlist[[species]][["Petal"]] <- 
#'     ggplot2::ggplot(iris[iris$Species == species, ],
#'                     ggplot2::aes(x = Petal.Length, y = Petal.Width)) + ggplot2::geom_point()
#' }
#' 
#' prep_filename_path(plt_list = nested_plotlist, 
#'                    file_format = "png")
#' }
#'
#' @keywords internal
#' @seealso \code{\link[gDRplots]{prep_plot_chunk}}
#' 
#' @export
prep_filename_path <- function(plt_list,
                               prefix = NULL,
                               path_file = NULL,
                               file_format = NULL) {
  
  checkmate::assert_list(plt_list)
  checkmate::assert_string(prefix, null.ok = TRUE)
  checkmate::assert_string(path_file, null.ok = TRUE)
  checkmate::assert_string(file_format, null.ok = TRUE)
  
  ls_file_name <- lapply(seq_along(plt_list), function(nm) {
    lvl1_name <- ifelse(is.null(names(plt_list)[nm]), nm, names(plt_list)[nm]) # number on name
    
    if (inherits(plt_list[[nm]], "list")) {
      
      ls_nested <- lapply(seq_along(plt_list[[nm]]), function(i_nm) {
        lvl2_name <-
          ifelse(is.null(names(plt_list[[nm]])[i_nm]), i_nm, names(plt_list[[nm]])[i_nm]) # number on name
        # file name
        file_name <- paste0(prefix,
                            neutralize_spaces(as.character(
                              ifelse(is.null(file_format), lvl2_name,
                                     paste(lvl2_name, file_format, sep = ".")))))
        # path
        ifelse(is.null(path_file), file_name, file.path(path_file, file_name))
      })
      
      if (is.null(names(plt_list[[nm]]))) {
        names(ls_nested) <- seq_along(plt_list[[nm]])
      } else {
        names(ls_nested) <- names(plt_list[[nm]])
      }
      
      ls_nested
    } else {
      # file name
      file_name <- paste0(prefix,
                          neutralize_spaces(as.character(
                            ifelse(is.null(file_format), lvl1_name,
                                   paste(lvl1_name, file_format, sep = ".")))))
      # path
      ifelse(is.null(path_file), file_name, file.path(path_file, file_name))
    }
  })
  
  if (is.null(names(plt_list))) {
    names(ls_file_name) <- seq_along(plt_list)
  } else {
    names(ls_file_name) <- names(plt_list)
  }
  # final
  ls_file_name
} 
