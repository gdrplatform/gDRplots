# Create list of qualitative colors

Create list of qualitative colors

## Usage

``` r
get_qual_colors(n = NULL)
```

## Arguments

- n:

  number of required colors

## Value

vector with hex colors from qualitative palettes

## Author

Janina Smoła <janina.smola@contractors.roche.com>

## Examples

``` r
get_qual_colors()
#>  [1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D"
#>  [8] "#666666" "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C"
#> [15] "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928" "#66C2A5"
#> [22] "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494" "#B3B3B3"
get_qual_colors(0)
#> [1] "#000000"
get_qual_colors(5)
#> [1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E"
get_qual_colors(35)
#>  [1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D"
#>  [8] "#666666" "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C"
#> [15] "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928" "#66C2A5"
#> [22] "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494" "#B3B3B3"
#> [29] "#52C098" "#FF8551" "#9C97DA" "#FF71AC" "#85C54F" "#FDC140" "#CD9B53"
```
