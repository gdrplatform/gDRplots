# Format Dose List for Legend

Helper function to create a formatted string of dose ranks and values
for the plot legend, breaking lines based on \`n_per_line\`.

## Usage

``` r
.format_dose_list(dose_vec, n_per_line)
```

## Arguments

- dose_vec:

  `numeric` vector of unique positive doses.

- n_per_line:

  `integer` number of items to display per line.

## Value

A single `character` string with HTML line breaks.

## Author

Bartosz Czech <czech.bartosz@external.gene.com>
