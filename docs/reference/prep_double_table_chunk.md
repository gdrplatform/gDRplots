# Prepare markdown chunk based on a doubly nested list of tables

Generates markdown code for displaying tables in a document using
\`knitr::knit()\`. Handles doubly nested lists, allowing for tabbed
sections for cell lines and then metrics. The inner header level (for
metrics) is automatically set to one level greater than the outer header
level.

## Usage

``` r
prep_double_table_chunk(
  tbl_list,
  chunk_name,
  dwn_list = NULL,
  header_level = 3,
  tabset_options = c("tabset", "tabset-dropdown"),
  sorting_opts = NULL
)
```

## Arguments

- tbl_list:

  A doubly nested named list of tables. The outer list represents cell
  lines, and the inner lists represent metrics. Names are used as
  headings.

- chunk_name:

  A character string specifying the base name for the generated code
  chunks. Avoid spaces.

- dwn_list:

  A named list of links to location (relative paths) where table or
  plots are saved, which when clocked, will be downloaded. It must have
  the same structure as `plt_list`.

- header_level:

  An integer specifying the markdown header level to use (e.g., 1 for
  \`#\`, 2 for \`##\`, etc.).

- tabset_options:

  A character vector of options for the tabset. This is only used when
  `plt_list` is a nested list. Possible values are "unnumbered",
  "tabset", and "tabset-dropdown".

- sorting_opts:

  A vector specifying global sorting options for all tables. Column
  names can be preceded by "-" to indicate descending order.

## Value

A list of character vectors. Each element corresponds to a cell line.
Each character vector represents markdown code for the cell line's
tabset.

## Author

Bartosz Czech <czech.bartosz@external.gene.com>

## Examples

``` r
nested_tables <- list(
  CellLine1 = list(MetricA = mtcars[1:5, ], MetricB = mtcars[6:10, ]),
  CellLine2 = list(MetricC = iris[1:5, ], MetricD = iris[6:10, ])
)
sorting_options <- c("cyl", "-hp")
prep_double_table_chunk(nested_tables, "nested_tables", header_level = 2,
  tabset_options = "tabset", sorting_opts = sorting_options)
#> [[1]]
#> [1] "## CellLine1 {.tabset}\n\n"                                                                                                                                                                                                                                                                                                                                         
#> [2] "### MetricA\n\n```{r nested_tables_CellLine1_MetricA, echo = FALSE}\nDT::formatRound(generate_datatable(nested_tables[[\"CellLine1\"]][[\"MetricA\"]], options = list(scrollX = TRUE, dom = \"t\", order = list(list(2L, \"asc\"), list(4L, \"desc\")))), columns = names(Filter(is.numeric, nested_tables[[\"CellLine1\"]][[\"MetricA\"]])), digits = 5) \n```\n\n"
#> [3] "### MetricB\n\n```{r nested_tables_CellLine1_MetricB, echo = FALSE}\nDT::formatRound(generate_datatable(nested_tables[[\"CellLine1\"]][[\"MetricB\"]], options = list(scrollX = TRUE, dom = \"t\", order = list(list(2L, \"asc\"), list(4L, \"desc\")))), columns = names(Filter(is.numeric, nested_tables[[\"CellLine1\"]][[\"MetricB\"]])), digits = 5) \n```\n\n"
#> 
#> [[2]]
#> [1] "## CellLine2 {.tabset}\n\n"                                                                                                                                                                                                                                                                                    
#> [2] "### MetricC\n\n```{r nested_tables_CellLine2_MetricC, echo = FALSE}\nDT::formatRound(generate_datatable(nested_tables[[\"CellLine2\"]][[\"MetricC\"]], options = list(scrollX = TRUE, dom = \"t\")), columns = names(Filter(is.numeric, nested_tables[[\"CellLine2\"]][[\"MetricC\"]])), digits = 5) \n```\n\n"
#> [3] "### MetricD\n\n```{r nested_tables_CellLine2_MetricD, echo = FALSE}\nDT::formatRound(generate_datatable(nested_tables[[\"CellLine2\"]][[\"MetricD\"]], options = list(scrollX = TRUE, dom = \"t\")), columns = names(Filter(is.numeric, nested_tables[[\"CellLine2\"]][[\"MetricD\"]])), digits = 5) \n```\n\n"
#> 
```
