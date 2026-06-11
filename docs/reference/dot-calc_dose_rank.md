# Calculate Dose Ranks

Helper function to rank unique positive concentration values. Zero
values or NAs are assigned a rank of "-".

## Usage

``` r
.calc_dose_rank(vals)
```

## Arguments

- vals:

  `numeric` vector of concentration values.

## Value

A `character` vector of ranks.

## Author

Bartosz Czech <czech.bartosz@external.gene.com>
