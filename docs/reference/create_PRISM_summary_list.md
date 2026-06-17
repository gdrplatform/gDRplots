# Create a list of PRISM association table

Create a list of PRISM association table

## Usage

``` r
create_PRISM_summary_list(assoc_summary_RV, assoc_summary_GR = NULL)
```

## Arguments

- assoc_summary_RV:

  A `data.table` with associations for normalization type of "Relative
  Viability", outputted by
  [`gDRplots::prep_assoc_summary()`](https://gdrplatform.github.io/gDRplots/reference/prep_assoc_summary.md)

- assoc_summary_GR:

  A `data.table` with associations for normalization type of "GR Value",
  outputted by
  [`gDRplots::prep_assoc_summary()`](https://gdrplatform.github.io/gDRplots/reference/prep_assoc_summary.md)

## Value

A list of table split by drug name and normalization type

## Author

Janina Smoła <janina.smola@contractors.roche.com>
