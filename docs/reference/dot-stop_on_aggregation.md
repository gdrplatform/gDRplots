# stop wrapper for \`data.table::dcast\` to handle unexpected aggregation

idea from: https://github.com/Rdatatable/data.table/issues/5386

## Usage

``` r
.stop_on_aggregation(fname, formula)
```

## Arguments

- fname:

  string with the name of the function that failed to the
  [`data.table::dcast`](https://rdrr.io/pkg/data.table/man/dcast.data.table.html)
  aggregation

- formula:

  string with the formula used in
  [`data.table::dcast`](https://rdrr.io/pkg/data.table/man/dcast.data.table.html)

## Value

`NULL`

## Author

Arkadiusz Gladki <arkadiusz.gladki@contractors.roche.com>
