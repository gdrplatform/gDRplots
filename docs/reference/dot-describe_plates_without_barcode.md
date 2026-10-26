# Describe the plates read without a barcode

Describe the plates read without a barcode

## Usage

``` r
.describe_plates_without_barcode(plates)
```

## Arguments

- plates:

  table recorded by gDRimport in `corrections$plates_without_barcode`:
  one row per plate, with `file`, `plate` and `used`, the barcode the
  plate was imported under (`NA` when it was skipped); or `NULL`.

## Value

A character vector with a sentence for the plates imported under another
barcode and one for the skipped plates, or `NULL` when `plates` is
empty.
