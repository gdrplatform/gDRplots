# Create a list of PRISM association table

Create a list of PRISM association table

## Usage

``` r
create_PRISM_summary_list(assoc_summary_RV, assoc_summary_GR = NULL)
```

## Arguments

- assoc_summary_RV:

  A `data.table` with associations for normalization type of "Relative
  Viability", outputted by
  [`gDRplots::prep_assoc_summary()`](https://gdrplatform.github.io/gDRplots/reference/prep_assoc_summary.md)

- assoc_summary_GR:

  A `data.table` with associations for normalization type of "GR Value",
  outputted by
  [`gDRplots::prep_assoc_summary()`](https://gdrplatform.github.io/gDRplots/reference/prep_assoc_summary.md)

## Value

A list of table split by drug name and normalization type

## Author

Janina Smoła <janina.smola@external.roche.com>

## Examples

``` r
dir_path <- system.file("testdata", package = "gDRplots")
ls_file <- list.files(dir_path, pattern = "[.]xlsx$")
assoc_tab <- prep_assoc_summary(dir_path = dir_path, ls_file = ls_file)
create_PRISM_summary_list(
  assoc_summary_RV = assoc_tab[grepl("RV", assoc_tab[["src"]]), ],
  assoc_summary_GR = assoc_tab[grepl("GR", assoc_tab[["src"]]), ]
)
#> $drug_001
#> $drug_001$RV
#>     feat_meta    feature      response   rho q_value neglog_q_value drug_grid
#>        <char>     <char>        <char> <num>   <num>          <num>    <char>
#>  1:   featNUX NU_006_X1F  RV_gDR_x_max  0.05  0.0035       2.455932  drug_001
#>  2:   featNUX NU_021_X1U  RV_gDR_x_max  0.05  0.0035       2.455932  drug_001
#>  3:   featNUX NU_002_X1B  RV_gDR_x_max  0.40  0.0060       2.221849  drug_001
#>  4:   featNUX NU_025_X1Y  RV_gDR_x_max  0.05  0.0110       1.958607  drug_001
#>  5:   featNUX NU_009_X1I  RV_gDR_x_max -1.00  0.0135       1.869666  drug_001
#>  6:   featNUX NU_020_X1T  RV_gDR_x_max -0.30  0.0235       1.628932  drug_001
#>  7:   featNUX NU_010_X1J  RV_gDR_x_max  0.05  0.0285       1.545155  drug_001
#>  8:   featNUX NU_008_X1H  RV_gDR_x_max -0.65  0.0310       1.508638  drug_001
#>  9:   featNUX NU_001_X1A  RV_gDR_x_max -1.00  0.0335       1.474955  drug_001
#> 10:   featNUX NU_017_X1Q  RV_gDR_x_max -0.30  0.0460       1.337242  drug_001
#> 11:   featNUX NU_023_X1W RV_gDR_x_mean -0.15  0.0070       2.154902  drug_001
#> 12:   featNUX NU_017_X1Q RV_gDR_x_mean -0.50  0.0295       1.530178  drug_001
#> 13:   featNUX NU_011_X1K RV_gDR_x_mean  1.25  0.0355       1.449772  drug_001
#> 
#> $drug_001$GR
#>     feat_meta    feature          response   rho q_value neglog_q_value
#>        <char>     <char>            <char> <num>   <num>          <num>
#>  1:   metaGRP GRP_007_XC GR_gDR_log10_xc50  0.75  0.0190       1.721246
#>  2:   metaGRP GRP_006_XC      GR_gDR_x_max  0.05  0.0035       2.455932
#>  3:   metaGRP GRP_021_XC      GR_gDR_x_max  0.05  0.0035       2.455932
#>  4:   metaGRP GRP_002_XC      GR_gDR_x_max  0.40  0.0060       2.221849
#>  5:   metaGRP GRP_025_XC      GR_gDR_x_max  0.05  0.0110       1.958607
#>  6:   metaGRP GRP_009_XC      GR_gDR_x_max -1.00  0.0135       1.869666
#>  7:   metaGRP GRP_020_XC      GR_gDR_x_max -0.30  0.0235       1.628932
#>  8:   metaGRP GRP_010_XC      GR_gDR_x_max  0.05  0.0285       1.545155
#>  9:   metaGRP GRP_008_XC      GR_gDR_x_max -0.65  0.0310       1.508638
#> 10:   metaGRP GRP_001_XC      GR_gDR_x_max -1.00  0.0335       1.474955
#> 11:   metaGRP GRP_017_XC      GR_gDR_x_max -0.30  0.0460       1.337242
#>     drug_grid
#>        <char>
#>  1:  drug_001
#>  2:  drug_001
#>  3:  drug_001
#>  4:  drug_001
#>  5:  drug_001
#>  6:  drug_001
#>  7:  drug_001
#>  8:  drug_001
#>  9:  drug_001
#> 10:  drug_001
#> 11:  drug_001
#> 
#> 
#> $drug_002
#> $drug_002$RV
#>    feat_meta    feature      response   rho q_value neglog_q_value drug_grid
#>       <char>     <char>        <char> <num>   <num>          <num>    <char>
#> 1:   metaGRP GRP_016_XC RV_gDR_x_mean  0.55  0.0485       1.314258  drug_002
#> 
#> 
#> $drug_002_RV_gDR_x_mean.xlsx
#> $drug_002_RV_gDR_x_mean.xlsx$GR
#>    feat_meta    feature      response   rho q_value neglog_q_value
#>       <char>     <char>        <char> <num>   <num>          <num>
#> 1:   metaGRP GRP_016_XC RV_gDR_x_mean  0.55  0.0485       1.314258
#>                      drug_grid
#>                         <char>
#> 1: drug_002_RV_gDR_x_mean.xlsx
#> 
#> 
#> $drug_003
#> $drug_003$GR
#>    feat_meta    feature          response   rho q_value neglog_q_value
#>       <char>     <char>            <char> <num>   <num>          <num>
#> 1:   featNUX NU_023_X1W GR_gDR_log10_xc50 -0.15  0.0070       2.154902
#> 2:   featNUX NU_017_X1Q GR_gDR_log10_xc50 -0.50  0.0295       1.530178
#> 3:   featNUX NU_011_X1K GR_gDR_log10_xc50  1.25  0.0355       1.449772
#>    drug_grid
#>       <char>
#> 1:  drug_003
#> 2:  drug_003
#> 3:  drug_003
#> 
#> 
```
