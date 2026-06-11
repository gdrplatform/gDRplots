# Trim labels to the required number of characters

Trim labels to the required number of characters

## Usage

``` r
.trim_labels(lbls_vec, max_lbl_length = Inf)
```

## Arguments

- lbls_vec:

  character vector with labels to be trimmed

- max_lbl_length:

  numeric value for the maximum number of characters in the label; if
  set to Inf, no trimming will be performed

## Value

character vectors with trimmed labels

## Author

Janina Smoła <janina.smola@contractors.roche.com>

## Examples

``` r
if (FALSE) { # \dontrun{
ls_lbls <- c(
  "short_lbl", "short_lbl", "veryveryverylong",
  "long_duplicates|lbl_1AB", "long_duplicates|lbl_1AB", "long_duplicates|lbl_123"
)

.trim_labels(lbls_vec = ls_lbls)
.trim_labels(lbls_vec = ls_lbls, max_lbl_length = 15)
} # }
```
