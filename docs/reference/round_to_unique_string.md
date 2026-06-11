# Round numbers to unique string

Rounds a numeric vector to the minimum precision needed to ensure
uniqueness

## Usage

``` r
round_to_unique_string(num_vec, initial_digits = 4)
```

## Arguments

- num_vec:

  a numeric vector to be rounded and change into character

- initial_digits:

  numeric value for number of decimal places to start rounding with

## Value

a character vector of unique numeric strings.

## Author

Janina Smoła <janina.smola@contractors.roche.com>

## Examples

``` r
vec <- c(0.00000000, 0.00000256, 0.00001280, 0.00006400, 0.00032000,
         0.00160000, 0.00800000, 0.04000000, 0.20000000, 1.00000000)

round_to_unique_string(vec)
#>  [1] "0.000000" "0.000003" "0.000013" "0.0001"   "0.0003"   "0.0016"  
#>  [7] "0.0080"   "0.0400"   "0.2000"   "1.0000"  
```
