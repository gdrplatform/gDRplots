# Compare two screening conditions of a chemical genomics screen and perform GSEA

Compares the response to each library drug between two conditions, each
a cell line at one co-treatment concentration (`0` for the drug alone),
and runs GSEA over the drug classes (`drug_moa`). The conditions may use
different cell lines and different co-treatment concentrations, for
example a resistant line kept in the anchor drug against its parental
line in standard medium. The difference is `condition1` minus
`condition2`.

## Usage

``` r
analyze_cgs_conditions(
  dt_metrics,
  condition1,
  condition2,
  metrics = "x_mean",
  normalization_type = "RV"
)
```

## Arguments

- dt_metrics:

  A data.table of combination metrics, as for
  [`analyze_cgs`](https://gdrplatform.github.io/gDRplots/reference/analyze_cgs.md),
  with `cotrt_value`.

- condition1, condition2:

  Named lists with `cell_line` (a value of `CellLineName`) and
  `cotrt_value` (a co-treatment concentration of that cell line, `0` for
  the drug alone).

- metrics:

  A character vector of metrics to compare: "x_mean", "x_AOC_range",
  "xc50", "ec50", "x_max". `xc50` is compared on the log10 scale.

- normalization_type:

  A string, "RV" or "GR".

## Value

A list of one element, named `"<condition1> vs. <condition2>"`, with
`fgsea`, `metrics_diff` and `moa_list` as returned by
[`analyze_cgs`](https://gdrplatform.github.io/gDRplots/reference/analyze_cgs.md),
so
[`plot_cgs_ranking`](https://gdrplatform.github.io/gDRplots/reference/plot_cgs_ranking.md)
can draw it.

## Examples

``` r
dt_metrics <- qs2::qs_read(system.file("testdata/cgs_data.qs2", package = "gDRplots"))
analyze_cgs_conditions(dt_metrics,
                       condition1 = list(cell_line = "CellLineName_1", cotrt_value = 1),
                       condition2 = list(cell_line = "CellLineName_2", cotrt_value = 0))
#> $`CellLineName_1 (co-treatment 1) vs. CellLineName_2 (single agent)`
#> $`CellLineName_1 (co-treatment 1) vs. CellLineName_2 (single agent)`$fgsea
#> $`CellLineName_1 (co-treatment 1) vs. CellLineName_2 (single agent)`$fgsea$x_mean
#> Empty data.table (0 rows and 9 cols): pathway,pval,padj,log2err,ES,NES...
#> 
#> 
#> $`CellLineName_1 (co-treatment 1) vs. CellLineName_2 (single agent)`$metrics_diff
#> Key: <DrugName, drug_moa>
#>         DrugName     drug_moa x_mean_c1 x_mean_c2      x_mean
#>           <char>       <char>     <num>     <num>       <num>
#>  1: DrugName_120 drug_moa_308 0.9769556 1.0823125 -0.10535694
#>  2: DrugName_137 drug_moa_308 0.8933290 0.9787796 -0.08545061
#>  3: DrugName_154 drug_moa_308 0.9236561 0.9668918 -0.04323572
#>  4:  DrugName_28 drug_moa_308 0.6991135 0.7101604 -0.01104683
#>  5: DrugName_281 drug_moa_308 0.7974026 0.8092212 -0.01181856
#>  6: DrugName_303 drug_moa_308 0.5696739 0.6981507 -0.12847688
#>  7: DrugName_381 drug_moa_308 0.9178351 0.9828750 -0.06503987
#>  8: DrugName_404 drug_moa_308 0.4606015 0.6073132 -0.14671165
#>  9:  DrugName_41 drug_moa_308 0.8491840 0.9568987 -0.10771467
#> 10: DrugName_626 drug_moa_308 0.9311594 0.9116004  0.01955895
#> 11: DrugName_667 drug_moa_308 0.9760000 1.0029111 -0.02691114
#>     normalization_type
#>                 <char>
#>  1:                 RV
#>  2:                 RV
#>  3:                 RV
#>  4:                 RV
#>  5:                 RV
#>  6:                 RV
#>  7:                 RV
#>  8:                 RV
#>  9:                 RV
#> 10:                 RV
#> 11:                 RV
#> 
#> $`CellLineName_1 (co-treatment 1) vs. CellLineName_2 (single agent)`$moa_list
#> $`CellLineName_1 (co-treatment 1) vs. CellLineName_2 (single agent)`$moa_list$drug_moa_308
#>  [1] "DrugName_120" "DrugName_137" "DrugName_154" "DrugName_28"  "DrugName_281"
#>  [6] "DrugName_303" "DrugName_381" "DrugName_404" "DrugName_41"  "DrugName_626"
#> [11] "DrugName_667"
#> 
#> 
#> 
```
