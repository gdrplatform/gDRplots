# Prepare markdown chunk with download link

Generate markdown code with a html link item that, when clicked,
downloads a file. The function output should be wrapped in
[`knitr::knit()`](https://rdrr.io/pkg/knitr/man/knit.html).

## Usage

``` r
create_download_link(dwn_path, link_txt = "Download Table")
```

## Arguments

- dwn_path:

  string with relative path to file with plot to be downloaded

- link_txt:

  string with text describing link

## Value

string with html download code

## See also

[`knitr::knit`](https://rdrr.io/pkg/knitr/man/knit.html)

## Author

Janina Smoła <janina.smola@contractors.roche.com>
