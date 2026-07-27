# Get color palette for the smooth values

Get color palette for the smooth values

## Usage

``` r
.get_smooth_palette(no_breaks)
```

## Arguments

- no_breaks:

  numeric number of breaks on scale

## Value

gDR palette for smooth values with given `no_breaks`

## Examples

``` r
# \donttest{
gDRplots:::.get_smooth_palette(25)
#>  [1] "#251739" "#2D1C46" "#352153" "#3D2760" "#462C6D" "#4E3179" "#563787"
#>  [8] "#5E3C93" "#6742A1" "#704DA7" "#7A59AE" "#8464B5" "#8E70BC" "#977CC2"
#> [15] "#A187C9" "#AB93D0" "#B59FD7" "#BCA9DA" "#C4B3DD" "#CBBEE1" "#D3C8E4"
#> [22] "#DBD2E7" "#E2DDEB" "#EAE7EE" "#F2F2F2"
# }
```
