# Prepare markdown chunk based on the nested plots list

Function output should be generated with
`knitr::knit(text = unlist(<result>))`

## Usage

``` r
prep_nested_plot_chunk(
  plt_list,
  chunk_name,
  link_list = NULL,
  dwn_list = NULL,
  header_level = 2
)
```

## Arguments

- plt_list:

  named list with generated plots to be shown in tabs; list of plots in
  nested hierarchy, where last 4th level is plot and 3rd level is
  `normalization_type` described by one of: "GR" ("GR Value") or "RV"
  ("Relative Viability")#'

- chunk_name:

  A character string specifying the base name for the generated code
  chunks. Avoid spaces.

- link_list:

  A named list of links to the location (relative paths) where plots are
  saved, which when clicked, will be displayed in a new browser tab. It
  must have the same structure as `plt_list`.

- dwn_list:

  A named list of links to location (relative paths) where table or
  plots are saved, which when clocked, will be downloaded. It must have
  the same structure as `plt_list`.

- header_level:

  An integer specifying the markdown header level to use (e.g., 1 for
  \`#\`, 2 for \`##\`, etc.).

## Value

list of character vectors - input for
[`knitr::knit`](https://rdrr.io/pkg/knitr/man/knit.html)

## See also

[`knitr::knit`](https://rdrr.io/pkg/knitr/man/knit.html)

## Author

Janina Smoła <janina.smola@contractors.roche.com>

## Examples

``` r
if (FALSE) { # \dontrun{
mae <- gDRutils::get_synthetic_data("small")
se <- mae[[gDRutils::get_supported_experiments("sa")]]
dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")

# help function
plot_col <- function(tab_plt, norm_type, col = "red") {
  tab_plt <- data.table::melt(
    data = tab_plt[normalization_type == norm_type][, c("rId", "xc50", "x_mean", "x_max")],
    id = "rId")
  plt <- ggplot2::ggplot(tab_plt, ggplot2::aes(x = variable, y = value)) +
    ggplot2::geom_col(fill = col)
  return(plt)
}

# creating nested list with plots
plotlist <- list()
ls_color <- c("darkred", "orange", "darkcyan")
for (drug in unique(dt_metrics$DrugName)) {
  for (cl in unique(dt_metrics$CellLineName)) {
    tab_plot <- dt_metrics[DrugName == drug & CellLineName == cl]

    plt_GR <- lapply(ls_color, function(col) plot_col(tab_plot, "RV", col))
    names(plt_GR) <- sprintf("%s_%s", "GR", ls_color)
    plt_RV <- lapply(ls_color, function(col) plot_col(tab_plot, "RV", col))
    names(plt_RV) <- sprintf("%s_%s", "RV", ls_color)

    plotlist[[drug]][[cl]][["RV"]] <- plt_RV
    plotlist[[drug]][[cl]][["GR"]] <- plt_GR
  }
}

prep_nested_plot_chunk(plotlist, "metric_value")

prep_nested_plot_chunk(plotlist, "metric_value")

} # }
```
