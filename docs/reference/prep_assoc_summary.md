# Prepare summary table with statistically significant associations

Prepare summary table with statistically significant associations

## Usage

``` r
prep_assoc_summary(
  dir_path,
  ls_file,
  alpha = 0.05,
  n_stat_sig_row = 10,
  read_file_fun = readxl::read_excel,
  as_list = FALSE
)
```

## Arguments

- dir_path:

  A string path to the directory containing files with associations
  data.

- ls_file:

  A character vector with names of files containing associations data.

- alpha:

  A numeric cutoff to identify statistically significant correlations

- n_stat_sig_row:

  A numeric value specifying the maximum number of statistically
  significant associations (rows) to include from each file.

- read_file_fun:

  A function to read the data from file; default is
  [`readxl::read_excel`](https://readxl.tidyverse.org/reference/read_excel.html)

- as_list:

  A logical flag indicating whether the result should be returned as a
  list or as a table.

## Value

A [`DT::datatable`](https://rdrr.io/pkg/DT/man/datatable.html) object.

## Author

Janina Smoła <janina.smola@contractors.roche.com>
