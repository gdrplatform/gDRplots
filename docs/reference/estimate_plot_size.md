# Estimate the optimal plot size (either ggplot or pheatmap) for saving plots

Estimate the optimal plot size (either ggplot or pheatmap) for saving
plots

## Usage

``` r
estimate_plot_size(plt, base_width = 10, base_height = 6, scale_factor = 0.5)
```

## Arguments

- plt:

  a ggplot or pheatmap object

- base_width:

  an integer with default base_width

- base_height:

  an integer with default base_height

- scale_factor:

  an integer with default scale_factor

## Value

named vector with optimal width and height used in the
[`ggplot2::ggsave`](https://ggplot2.tidyverse.org/reference/ggsave.html)
function

## Author

Bartosz Czech <czech.bartosz@external.gene.com>

## Examples

``` r
p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) + ggplot2::geom_point()
estimate_plot_size(p)
#>  width height 
#>   10.5    6.5 
```
