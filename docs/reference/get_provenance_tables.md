# Tables describing how the results of a gDR object were produced

Builds the tables shown in the "How These Numbers Were Produced" section
of the processing report from the provenance that gDR records in each
experiment (see
[`get_SE_provenance`](https://gdrplatform.github.io/gDRstyle/reference/fit_provenance.html)):
one row per fit, with the function, its version, its settings and the
model, and one row per processing step, described in a sentence.

## Usage

``` r
get_provenance_tables(x)
```

## Arguments

- x:

  `MultiAssayExperiment` or `SummarizedExperiment`.

## Value

A list with two `data.table`s:

- `fits`: one row per fit, with columns `assay`, `fit_source`, `fn`,
  `settings` and `model`.

- `steps`: one row per processing step, with columns `step` and
  `description`.

For a `MultiAssayExperiment`, both tables start with an `experiment`
column naming the experiment of each row; a `SummarizedExperiment` is a
single experiment, so it has no such column.

Both tables are empty for objects without provenance.

## Author

Bartosz Czech <czech.bartosz@external.gene.com>

## Examples

``` r
mae <- gDRutils::get_synthetic_data("small")
tables <- get_provenance_tables(mae)
tables$fits
#> Empty data.table (0 rows and 6 cols): experiment,assay,fit_source,fn,settings,model
```
