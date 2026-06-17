# Prepare list of names or paths

This function can creates lists that can be used as an input to
[`prep_plot_chunk`](https://gdrplatform.github.io/gDRplots/reference/prep_plot_chunk.md) -
param `link_list` and `dwn_list`. Depending on the inputs, outputted
string will have a format: *\<path_file\>\<prefix\>\<name of item on the
list\>***.***\<file_format\>* Note: An unnamed item in the list will be
numbered.

## Usage

``` r
prep_filename_path(
  plt_list,
  prefix = NULL,
  path_file = NULL,
  file_format = NULL
)
```

## Arguments

- plt_list:

  A named list of plots. Names will be used as headings for plots/tab
  groups. If unnamed, ordinal numbers will be used. Can be nested lists
  for tabbed output. If a nested list is provided, the inner lists
  should also be named.

- prefix:

  string to be added as prefix to the name of file

- path_file:

  string path to directory, where data will be read from

- file_format:

  string specifying the format of file

## Value

list of strings describing names or paths depending on the input; the
list is structured like the input `plt_list` and retains the same names
as `plt_list`.

## See also

[`prep_plot_chunk`](https://gdrplatform.github.io/gDRplots/reference/prep_plot_chunk.md)

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

prep_filename_path(plt_list = plotlist,
                   prefix = "iris__",
                   path_file = file.path(".", "plots"))

# Nested list of plots for tabbed output
nested_plotlist <- list()
for (species in unique(iris$Species)) {
  nested_plotlist[[species]] <- list()
  nested_plotlist[[species]][["Sepal"]] <-
    ggplot2::ggplot(iris[iris$Species == species, ],
                    ggplot2::aes(x = Sepal.Length, y = Sepal.Width)) + ggplot2::geom_point()
  nested_plotlist[[species]][["Petal"]] <-
    ggplot2::ggplot(iris[iris$Species == species, ],
                    ggplot2::aes(x = Petal.Length, y = Petal.Width)) + ggplot2::geom_point()
}

prep_filename_path(plt_list = nested_plotlist,
                   file_format = "png")
} # }
```
