# Render a data.table as a text-based ggplot table

Render a data.table as a text-based ggplot table

## Usage

``` r
.table_to_ggplot(dt, base_size = 10, digits = 3)
```

## Arguments

- dt:

  data.table to render

- base_size:

  numeric base font size (default 10)

## Value

ggplot object displaying the table

## Examples

``` r
# \donttest{
dt <- data.table::data.table(drug = c("A", "B"), r2 = c(0.9, 0.8))
gDRplots:::.table_to_ggplot(dt)

# }
```
