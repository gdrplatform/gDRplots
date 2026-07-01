# Prepare markdown chunk based on a list of plots

Generates markdown code for displaying plots in a document using
[`knitr::knit()`](https://rdrr.io/pkg/knitr/man/knit.html). The function
handles both simple lists of plots and nested lists, allowing for the
creation of tabbed sections for grouped plots.

## Usage

``` r
prep_plot_chunk(
  plt_list,
  chunk_name,
  link_list = NULL,
  dwn_list = NULL,
  saved_plot_list = NULL,
  header_level = 3,
  tabset_options = c("tabset", "tabset-dropdown")
)
```

## Arguments

- plt_list:

  A named list of plots. Names will be used as headings for plots/tab
  groups. If unnamed, ordinal numbers will be used. Can be nested lists
  for tabbed output. If a nested list is provided, the inner lists
  should also be named.

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

- saved_plot_list:

  A character vector or list of absolute paths to pre-saved plot files
  (e.g. SVG). When provided, plots are embedded from these files instead
  of re-rendering the ggplot objects, which avoids expensive
  double-rendering and can dramatically speed up report generation. Must
  have the same structure as `plt_list`.

- header_level:

  An integer specifying the markdown header level to use (e.g., 1 for
  \`#\`, 2 for \`##\`, etc.).

- tabset_options:

  A character vector of options for the tabset. This is only used when
  `plt_list` is a nested list. Possible values are "unnumbered",
  "tabset", and "tabset-dropdown".

## Value

A list of character vectors. Each element of the list corresponds to a
plot or a group of plots (if `plt_list` is nested). Each character
vector within the list represents the markdown code to be processed by
[`knitr::knit()`](https://rdrr.io/pkg/knitr/man/knit.html). For nested
lists, each element will contain a "header" element and an "items"
element. The "header" element is the markdown for the tabset header, and
the "items" element is a list of markdown chunks for each tab.

## See also

[`knitr::knit`](https://rdrr.io/pkg/knitr/man/knit.html)

## Author

Janina Smoła <janina.smola@contractors.roche.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# Simple list of plots
plotlist <- lapply(unique(iris$Species), function(iris_name) {
  ggplot2::ggplot(iris[iris$Species == iris_name, c("Sepal.Length", "Sepal.Width")]) +
  ggplot2::geom_point(ggplot2::aes(x = Sepal.Length, y = Sepal.Width))
})
names(plotlist) <- unique(iris$Species)

prep_plot_chunk(plotlist, "iris")

# Nested list of plots for tabbed output
nested_plotlist <- list()
for (species in unique(iris$Species)) {
  nested_plotlist[[species]] <- list()
  nested_plotlist[[species]][["Sepal"]] <- ggplot2::ggplot(iris[iris$Species == species, ],
    ggplot2::aes(x = Sepal.Length, y = Sepal.Width)) + ggplot2::geom_point()
  nested_plotlist[[species]][["Petal"]] <- ggplot2::ggplot(iris[iris$Species == species, ],
    ggplot2::aes(x = Petal.Length, y = Petal.Width)) + ggplot2::geom_point()
}

prep_plot_chunk(nested_plotlist, "iris_nested", tabset_options = c("tabset", "unnumbered"))
} # }
```
