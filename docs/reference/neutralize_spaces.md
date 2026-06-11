# Replace spaces with another character

Replace spaces with another character

## Usage

``` r
neutralize_spaces(x, replacement = "_")
```

## Arguments

- x:

  String where matches are sought

- replacement:

  String replacement for spaces

## Value

String with spaces replaced by the specified character

## Examples

``` r
neutralize_spaces("GDC-123|Abc x G01234")
#> [1] "GDC-123|Abc_x_G01234"
neutralize_spaces("MNO-321P 789R YY#1 ")
#> [1] "MNO-321P_789R_YY#1"
neutralize_spaces("drug_001 x drug_002", ".")
#> [1] "drug_001.x.drug_002"
```
