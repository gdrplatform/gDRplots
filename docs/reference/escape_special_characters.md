# Escape colon and hash

Escape colon and hash

## Usage

``` r
escape_special_characters(x)
```

## Arguments

- x:

  String

## Value

Original string with `:`s and `#`s and `/`s escaped

## Examples

``` r
escape_special_characters("ABC:123")
#> [1] "ABC[colon]123"
escape_special_characters("AD_12")
#> [1] "AD_12"
escape_special_characters("AD#12")
#> [1] "AD[hash]12"
escape_special_characters("AD/12")
#> [1] "AD[slash]12"
```
