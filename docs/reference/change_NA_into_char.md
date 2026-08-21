# Change NA into given string

Change NA into given string

## Usage

``` r
change_NA_into_char(x, lbl_NA = "NA")
```

## Arguments

- x:

  vector with items suspected of being NA

- lbl_NA:

  string - replacement for NA - as default "NA"

## Value

character (for NA -\> given string)

## Author

Janina Smoła <janina.smola@external.roche.com>

## Examples

``` r
change_NA_into_char(c(1, NA, 3))
#> [1] "1"  "NA" "3" 
change_NA_into_char(c("a", NA, "b"), lbl_NA = "missing")
#> [1] "a"       "missing" "b"      
```
