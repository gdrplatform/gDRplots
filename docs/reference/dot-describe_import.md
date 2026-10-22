# Describe the import step recorded by gDRimport

Describe the import step recorded by gDRimport

## Usage

``` r
.describe_import(settings)
```

## Arguments

- settings:

  list; the record that gDRimport attaches to the loaded data (see
  `gDRimport::get_import_provenance`): `gDRimport` (version),
  `instrument` or, for PRISM data, `source$format`, `files` (a table
  with a `role` column) and `corrections`.

## Value

A character vector of sentences, or `NULL` when `settings` is `NULL`.
