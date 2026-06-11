# Compute color for number font in pheatmap::pheatmap based on given color palette and breaks

Compute color for number font in pheatmap::pheatmap based on given color
palette and breaks

## Usage

``` r
.get_pheatmap_number_color(
  mat_with_metric,
  colors_vec,
  breaks,
  dark_color_font = "white",
  light_color_font = "black"
)
```

## Arguments

- mat_with_metric:

  numeric matrix with metric values; must have named rows and columns

- colors_vec:

  character vector of colors (valid name or hex) used in heatmap; must
  to be one item shorter than `no_breaks`

- breaks:

  numeric vector of breaks on scale used for mapping values to colors

- dark_color_font:

  string with valid color name of font for field with dark background

- light_color_font:

  string with valid color name of font for field without dark background

## Value

named `matrix` with number color

## Author

Janina Smoła <janina.smola@contractors.roche.com>

## Examples

``` r
if (FALSE) { # \dontrun{
mat <- matrix(-14:30, ncol = 5,
              dimnames = list(letters[1:9], LETTERS[1:5]))
no_breaks <- 15
breaks <- seq(from = min(mat), to = max(mat), length.out = no_breaks + 1)
ls_colors <- c("limegreen", "darkblue", "orange")
hm_colors <- grDevices::colorRampPalette(ls_colors)(no_breaks)

number_color <- .get_pheatmap_number_color(mat, hm_colors, breaks)

pheatmap::pheatmap(mat,
                   breaks = breaks,
                   color = hm_colors,
                   display_numbers = TRUE,
                   number_color = number_color,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE)
} # }
```
