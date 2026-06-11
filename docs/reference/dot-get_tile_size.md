# Calculate the size of tiles based on pos_x/pos_y values

Since
[`ggplot2::geom_tile`](https://ggplot2.tidyverse.org/reference/geom_tile.html)
uses the center of the tile and its size (x, y, width, height), x and y
are given as pos_x and pos_y, it is required to calculate the width and
height.

## Usage

``` r
.get_tile_size(pos_vec)
```

## Arguments

- pos_vec:

  numeric vector with pos_x or pos_y values

## Value

size of tile in
[`ggplot2::geom_tile`](https://ggplot2.tidyverse.org/reference/geom_tile.html)
