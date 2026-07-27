# Get color palette for the excess values

Get color palette for the excess values

## Usage

``` r
.get_excess_palette(no_breaks)
```

## Arguments

- no_breaks:

  numeric number of breaks on scale

## Value

gDR palette for excess values with given `no_breaks`

## Examples

``` r
# \donttest{
gDRplots:::.get_excess_palette(20)
#>  [1] "#3A5FCD" "#3D65DA" "#416BE7" "#4571F4" "#507CFE" "#7D9DFA" "#AABDF7"
#>  [8] "#D7DEF4" "#F2F2F2" "#F2F2F2" "#F2F2F2" "#F2F2F2" "#F4D3D3" "#F7A0A0"
#> [15] "#FA6D6D" "#FE3A3A" "#F42D2D" "#E72B2B" "#DA2828" "#CD2626"
# }
```
