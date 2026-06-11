# Check data type

Check data type

## Usage

``` r
.get_data_type(dt_, desc_col = NULL)
```

## Arguments

- dt\_:

  `data.table` with dependent variables data in the wide format, where
  rows are samples, columns are feature levels

- desc_col:

  a character vector with column names describing the data and which do
  not contain data itself

## Value

a string describing type of data - "numeric" or "categorical"

## Author

Janina Smoła <janina.smola@contractors.roche.com>

## Examples

``` r
if (FALSE) { # \dontrun{
tab_cat <- data.table::data.table(
  ID = sprintf("ID_%s", seq_len(5)),
  brown = c(0, 1, 1, 0, 0),
  blue = c(1, 0, NA, 0, 1),
  green = c(0, 0, 0, 1, 0)
)
.get_data_type(dt_ = tab_cat, desc_col = "ID")

tab_feat <- data.table::data.table(
  ID = sprintf("ID_%s", seq_len(5)),
  grp = LETTERS[seq_len(5)],
  low = c(0, 1, 1, NA, 0),
  med = c(1, 1, NA, 0, 1),
  high = c(0, 1, 0, 1, 0)
)
.get_data_type(dt_ = tab_feat, desc_col = c("ID", "grp"))
} # }
```
