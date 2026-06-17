# Analyze Chemical Genomics (CGS) Data and Perform GSEA

This function analyzes chemical genomics (CGS) data, filters by compound
mechanism of action (MOA), prepares the data for Gene Set Enrichment
Analysis (GSEA), and performs GSEA for specified cell lines and metrics.

## Usage

``` r
analyze_cgs(
  dt_metrics,
  metrics,
  cl_name = NULL,
  cellline1 = NULL,
  cellline2 = NULL,
  normalization_type = "RV"
)
```

## Arguments

- dt_metrics:

  A data.table containing screening data. Requires columns:
  \`drug_moa\`, \`DrugName\`, \`CellLineName\`, and columns specified in
  the \`metrics\` argument.

- metrics:

  A character vector specifying the response metrics to analyze:
  "x_mean", "x_AOC_range", "xc50", "ec50", "x_max".

- cl_name:

  An optional string specifying a single cell line to analyze. If NULL
  (default), all cell lines in the data are analyzed. Should be NULL if
  `cellline1` and `cellline2` are provided.

- cellline1:

  An optional string representing the first cell line name. Should be
  specified alongside `cellline2`.

- cellline2:

  An optional string representing the second cell line name. Should be
  specified alongside `cellline1`.

- normalization_type:

  A string with normalization types to be selected, one of: "GR"
  ("GRvalue") or "RV" ("RelativeViability"). Passed to
  [`prep_dt_response_metric_diff`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_response_metric_diff.md).

## Value

A list of results, where each element corresponds to a cell line or cell
line difference. Each result contains:

- `fgsea`: GSEA results for the specified metrics,

- `metrics_diff`: The prepped data.table used for the GSEA analysis,

- `moa_list`: A list of DrugNames grouped by drug_moa.

## Author

Bartosz Czech <czech.bartosz@external.gene.com>

## Examples

``` r
dt_metrics <- qs2::qs_read(system.file("testdata/cgs_data.qs2", package = "gDRplots"))
analyze_cgs(dt_metrics, metrics = "xc50", cl_name = "CellLineName_1")
#> In `dt_metrics` some xc50 values are infinite.
#> $CellLineName_1
#> $CellLineName_1$fgsea
#> $CellLineName_1$fgsea$xc50
#>         pathway       pval       padj    log2err         ES       NES  size
#>          <char>      <num>      <num>      <num>      <num>     <num> <int>
#> 1:  drug_moa_69 0.19299520 0.19299520 0.02411837  0.7000000  1.354846     4
#> 2: drug_moa_308 0.04650697 0.09301394 0.02184993 -0.9247266 -1.511269    10
#>                                                                           leadingEdge
#>                                                                                <list>
#> 1:                                 DrugName_352,DrugName_488,DrugName_43,DrugName_197
#> 2: DrugName_518,DrugName_28,DrugName_41,DrugName_452,DrugName_404,DrugName_220,...[7]
#>         median
#>          <num>
#> 1: -0.05227476
#> 2: -0.37216031
#> 
#> 
#> $CellLineName_1$metrics_diff
#>         rId    cId   CellLineName     DrugName DrugName_2     drug_moa
#>      <char> <char>         <char>       <char>     <char>       <char>
#>  1: rId_179  cId_1 CellLineName_1 DrugName_220 DrugName_1 drug_moa_308
#>  2: rId_199  cId_1 CellLineName_1 DrugName_675 DrugName_1 drug_moa_308
#>  3: rId_228  cId_1 CellLineName_1 DrugName_452 DrugName_1 drug_moa_308
#>  4: rId_313  cId_1 CellLineName_1 DrugName_352 DrugName_1  drug_moa_69
#>  5: rId_379  cId_1 CellLineName_1  DrugName_41 DrugName_1 drug_moa_308
#>  6: rId_471  cId_1 CellLineName_1 DrugName_488 DrugName_1  drug_moa_69
#>  7:  rId_53  cId_1 CellLineName_1 DrugName_197 DrugName_1  drug_moa_69
#>  8: rId_586  cId_1 CellLineName_1  DrugName_28 DrugName_1 drug_moa_308
#>  9:  rId_62  cId_1 CellLineName_1 DrugName_112 DrugName_1 drug_moa_308
#> 10: rId_624  cId_1 CellLineName_1 DrugName_518 DrugName_1 drug_moa_308
#> 11: rId_680  cId_1 CellLineName_1 DrugName_137 DrugName_1 drug_moa_308
#> 12: rId_682  cId_1 CellLineName_1 DrugName_404 DrugName_1 drug_moa_308
#> 13: rId_698  cId_1 CellLineName_1 DrugName_541 DrugName_1 drug_moa_308
#> 14: rId_740  cId_1 CellLineName_1  DrugName_43 DrugName_1  drug_moa_69
#>     RV_gDR_log10_xc50_cotrt_zero_1_drug_1 RV_gDR_log10_xc50_cotrt_1_drug_1
#>                                     <num>                            <num>
#>  1:                            -0.2819911                       -0.6363582
#>  2:                             2.2111867                        2.0686472
#>  3:                             0.7543934                        0.3360619
#>  4:                             1.3645538                        1.3590174
#>  5:                             0.9619740                        0.5435939
#>  6:                             0.9576023                        0.9232680
#>  7:                             0.8618681                        0.6784762
#>  8:                             1.1041723                        0.5346264
#>  9:                             0.9205277                        0.8163443
#> 10:                             1.9542957                        1.1662197
#> 11:                             1.0527963                        0.8604461
#> 12:                            -0.7559844                       -1.1459378
#> 13:                            -0.6727913                       -0.6809344
#> 14:                             1.3479973                        1.2777822
#>             xc50 normalization_type
#>            <num>             <char>
#>  1: -0.354367161                 RV
#>  2: -0.142539517                 RV
#>  3: -0.418331488                 RV
#>  4: -0.005536396                 RV
#>  5: -0.418380126                 RV
#>  6: -0.034334371                 RV
#>  7: -0.183391876                 RV
#>  8: -0.569545915                 RV
#>  9: -0.104183456                 RV
#> 10: -0.788075927                 RV
#> 11: -0.192350146                 RV
#> 12: -0.389953458                 RV
#> 13: -0.008143045                 RV
#> 14: -0.070215159                 RV
#> 
#> $CellLineName_1$moa_list
#> $CellLineName_1$moa_list$drug_moa_308
#>  [1] "DrugName_220" "DrugName_167" "DrugName_675" "DrugName_452" "DrugName_41" 
#>  [6] "DrugName_28"  "DrugName_112" "DrugName_518" "DrugName_137" "DrugName_404"
#> [11] "DrugName_541"
#> 
#> $CellLineName_1$moa_list$drug_moa_69
#> [1] "DrugName_352" "DrugName_488" "DrugName_197" "DrugName_43" 
#> 
#> 
#> 
cellline1 <- "CellLineName_1"
cellline2 <- "CellLineName_2"
analyze_cgs(dt_metrics, "xc50", cl_name = NULL, cellline1, cellline2)
#> In `dt_metrics` some xc50 values are infinite.
#> $`CellLineName_1 vs. CellLineName_2`
#> $`CellLineName_1 vs. CellLineName_2`$fgsea
#> $`CellLineName_1 vs. CellLineName_2`$fgsea$xc50
#> Empty data.table (0 rows and 9 cols): pathway,pval,padj,log2err,ES,NES...
#> 
#> 
#> $`CellLineName_1 vs. CellLineName_2`$metrics_diff
#> Key: <drug_moa_2, DrugName_2, drug_moa, DrugName>
#>        DrugName     drug_moa DrugName_2 drug_moa_2 CellLineName_c1
#>          <char>       <char>     <char>     <char>          <char>
#> 1: DrugName_137 drug_moa_308 DrugName_1 drug_moa_1  CellLineName_1
#> 2: DrugName_220 drug_moa_308 DrugName_1 drug_moa_1  CellLineName_1
#> 3: DrugName_255 drug_moa_308 DrugName_1 drug_moa_1  CellLineName_1
#> 4: DrugName_381 drug_moa_308 DrugName_1 drug_moa_1  CellLineName_1
#> 5: DrugName_452 drug_moa_308 DrugName_1 drug_moa_1  CellLineName_1
#> 6: DrugName_518 drug_moa_308 DrugName_1 drug_moa_1  CellLineName_1
#> 7: DrugName_675 drug_moa_308 DrugName_1 drug_moa_1  CellLineName_1
#>    CellLineName_c2  xc50
#>             <char> <num>
#> 1:  CellLineName_2     0
#> 2:  CellLineName_2     0
#> 3:  CellLineName_2     0
#> 4:  CellLineName_2     0
#> 5:  CellLineName_2     0
#> 6:  CellLineName_2     0
#> 7:  CellLineName_2     0
#> 
#> $`CellLineName_1 vs. CellLineName_2`$moa_list
#> $`CellLineName_1 vs. CellLineName_2`$moa_list$drug_moa_308
#> [1] "DrugName_137" "DrugName_220" "DrugName_255" "DrugName_381" "DrugName_452"
#> [6] "DrugName_518" "DrugName_675"
#> 
#> 
#> 
```
