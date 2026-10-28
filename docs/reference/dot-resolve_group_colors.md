# Resolve a color mapping for group names

Builds a named color vector for `group_names`, resolving `colors_vec` by
name wherever possible so each color maps to its intended group
regardless of the order or length of the inputs.

## Usage

``` r
.resolve_group_colors(group_names_requested, group_names, colors_vec = NULL)
```

## Arguments

- group_names_requested:

  character vector of group names in their originally requested order,
  used to pair an unnamed `colors_vec` by name; `NULL` if no
  `group_names` were requested. Duplicates are ignored (the first
  occurrence wins).

- group_names:

  character vector (no duplicates) of group names to produce a color
  for, in their required order

- colors_vec:

  character vector of colors (name or hex value). If fully named (e.g.
  `c(cell_DC = "blue")`), matched to `group_names` by name, independent
  of order; otherwise matched positionally against
  `group_names_requested` (partial naming is treated as unnamed).
  Invalid, missing, or extra entries fall back to the default palette
  ([`get_qual_colors`](https://gdrplatform.github.io/gDRplots/reference/get_qual_colors.md));
  nothing is recycled.

## Value

named character vector of colors, one per entry of `group_names`

## Author

Janina Smoła <janina.smola@external.roche.com>

## Examples

``` r
# \donttest{
gDRplots:::.resolve_group_colors(
 group_names_requested = c("cell_DC", "cell_XY", "cell_AB"),
 group_names = c("cell_AB", "cell_DC", "cell_XY"),
 colors_vec = c("blue", "red", "grey"))
#> cell_AB cell_DC cell_XY 
#>  "grey"  "blue"   "red" 

gDRplots:::.resolve_group_colors(
 group_names_requested = c("cell_DC", "cell_XY", "cell_AB"),
 group_names = c("cell_AB", "cell_DC", "cell_XY"),
 colors_vec = c(cell_DC = "blue", cell_XY = "red", cell_AB = "grey"))
#> cell_AB cell_DC cell_XY 
#>  "grey"  "blue"   "red" 

# missing entries fall back to the default palette
gDRplots:::.resolve_group_colors(
 group_names_requested = c("cell_DC", "cell_XY", "cell_AB"),
 group_names = c("cell_AB", "cell_DC", "cell_XY"),
 colors_vec = c("blue", "red"))
#>   cell_AB   cell_DC   cell_XY 
#> "#1B9E77"    "blue"     "red" 

# partial naming is treated as unnamed
gDRplots:::.resolve_group_colors(
 group_names_requested = c("cell_DC", "cell_XY", "cell_AB"),
 group_names = c("cell_AB", "cell_DC", "cell_XY"),
 colors_vec = c(cell_AB = "grey", "blue", "red"))
#> Warning: 'colors_vec' has some, but not all, elements named - names are ignored and colors are matched to requested group_names by position instead.
#> cell_AB cell_DC cell_XY 
#>   "red"  "grey"  "blue" 
# }
```
