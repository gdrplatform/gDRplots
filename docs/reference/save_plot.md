# Save gDR plots to a specified path

Save gDR plots to a specified path

## Usage

``` r
save_plot(plt, path, format = "svg")
```

## Arguments

- plt:

  A plot object; either a ggplot2 or pheatmap object.

- path:

  A string specifying the path where the plot should be saved.

- format:

  A string specifying the format for saving the plot; either "svg",
  "png", or "pdf". Default is "svg".

## Value

`NULL`

## See also

[`ggplot2::ggsave`](https://ggplot2.tidyverse.org/reference/ggsave.html)

## Examples

``` r
tmp_dir <- file.path(tempdir(), "plot_dir")
dir.create(tmp_dir, showWarnings = FALSE)
p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) + ggplot2::geom_point()
save_plot(plt = p, path = paste(tmp_dir, "mtcars_scatter", sep = "/"), format = "png")
```
