# Prepare markdown chunk with zoom link

Generate markdown code for html link item which when clicked opens plot
in new browser. The function output should be wrapped in
[`knitr::knit()`](https://rdrr.io/pkg/knitr/man/knit.html).

## Usage

``` r
create_zoom_link(img_path, link_txt = "Zoom In for Details")
```

## Arguments

- img_path:

  string with relative path to file with plot to be shown

- link_txt:

  string with text describing link

## Value

string with html link code

## See also

[`knitr::knit`](https://rdrr.io/pkg/knitr/man/knit.html)

## Author

Janina Smoła <janina.smola@contractors.roche.com>
