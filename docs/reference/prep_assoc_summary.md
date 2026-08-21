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

Janina Smoła <janina.smola@external.roche.com>

## Examples

``` r
dir_path <- system.file("testdata", package = "gDRplots")
ls_file <- list.files(dir_path, pattern = "[.]xlsx$")
prep_assoc_summary(dir_path = dir_path, ls_file = ls_file)
#>        feature          response   rho q_value neglog_q_value
#>         <char>            <char> <num>   <num>          <num>
#>  1: NU_023_X1W GR_gDR_log10_xc50 -0.15  0.0070       2.154902
#>  2: NU_017_X1Q GR_gDR_log10_xc50 -0.50  0.0295       1.530178
#>  3: NU_011_X1K GR_gDR_log10_xc50  1.25  0.0355       1.449772
#>  4: GRP_007_XC GR_gDR_log10_xc50  0.75  0.0190       1.721246
#>  5: GRP_006_XC      GR_gDR_x_max  0.05  0.0035       2.455932
#>  6: GRP_021_XC      GR_gDR_x_max  0.05  0.0035       2.455932
#>  7: GRP_002_XC      GR_gDR_x_max  0.40  0.0060       2.221849
#>  8: GRP_025_XC      GR_gDR_x_max  0.05  0.0110       1.958607
#>  9: GRP_009_XC      GR_gDR_x_max -1.00  0.0135       1.869666
#> 10: GRP_020_XC      GR_gDR_x_max -0.30  0.0235       1.628932
#> 11: GRP_010_XC      GR_gDR_x_max  0.05  0.0285       1.545155
#> 12: GRP_008_XC      GR_gDR_x_max -0.65  0.0310       1.508638
#> 13: GRP_001_XC      GR_gDR_x_max -1.00  0.0335       1.474955
#> 14: GRP_017_XC      GR_gDR_x_max -0.30  0.0460       1.337242
#> 15: NU_006_X1F      RV_gDR_x_max  0.05  0.0035       2.455932
#> 16: NU_021_X1U      RV_gDR_x_max  0.05  0.0035       2.455932
#> 17: NU_002_X1B      RV_gDR_x_max  0.40  0.0060       2.221849
#> 18: NU_025_X1Y      RV_gDR_x_max  0.05  0.0110       1.958607
#> 19: NU_009_X1I      RV_gDR_x_max -1.00  0.0135       1.869666
#> 20: NU_020_X1T      RV_gDR_x_max -0.30  0.0235       1.628932
#> 21: NU_010_X1J      RV_gDR_x_max  0.05  0.0285       1.545155
#> 22: NU_008_X1H      RV_gDR_x_max -0.65  0.0310       1.508638
#> 23: NU_001_X1A      RV_gDR_x_max -1.00  0.0335       1.474955
#> 24: NU_017_X1Q      RV_gDR_x_max -0.30  0.0460       1.337242
#> 25: NU_023_X1W     RV_gDR_x_mean -0.15  0.0070       2.154902
#> 26: NU_017_X1Q     RV_gDR_x_mean -0.50  0.0295       1.530178
#> 27: NU_011_X1K     RV_gDR_x_mean  1.25  0.0355       1.449772
#> 28: GRP_016_XC     RV_gDR_x_mean  0.55  0.0485       1.314258
#>        feature          response   rho q_value neglog_q_value
#>         <char>            <char> <num>   <num>          <num>
#>                                                       src
#>                                                    <char>
#>  1: tab_assoc_GR__featNUX_drug_003_GR_gDR_log10_xc50.xlsx
#>  2: tab_assoc_GR__featNUX_drug_003_GR_gDR_log10_xc50.xlsx
#>  3: tab_assoc_GR__featNUX_drug_003_GR_gDR_log10_xc50.xlsx
#>  4: tab_assoc_GR__metaGRP_drug_001_GR_gDR_log10_xc50.xlsx
#>  5:      tab_assoc_GR__metaGRP_drug_001_GR_gDR_x_max.xlsx
#>  6:      tab_assoc_GR__metaGRP_drug_001_GR_gDR_x_max.xlsx
#>  7:      tab_assoc_GR__metaGRP_drug_001_GR_gDR_x_max.xlsx
#>  8:      tab_assoc_GR__metaGRP_drug_001_GR_gDR_x_max.xlsx
#>  9:      tab_assoc_GR__metaGRP_drug_001_GR_gDR_x_max.xlsx
#> 10:      tab_assoc_GR__metaGRP_drug_001_GR_gDR_x_max.xlsx
#> 11:      tab_assoc_GR__metaGRP_drug_001_GR_gDR_x_max.xlsx
#> 12:      tab_assoc_GR__metaGRP_drug_001_GR_gDR_x_max.xlsx
#> 13:      tab_assoc_GR__metaGRP_drug_001_GR_gDR_x_max.xlsx
#> 14:      tab_assoc_GR__metaGRP_drug_001_GR_gDR_x_max.xlsx
#> 15:      tab_assoc_RV__featNUX_drug_001_RV_gDR_x_max.xlsx
#> 16:      tab_assoc_RV__featNUX_drug_001_RV_gDR_x_max.xlsx
#> 17:      tab_assoc_RV__featNUX_drug_001_RV_gDR_x_max.xlsx
#> 18:      tab_assoc_RV__featNUX_drug_001_RV_gDR_x_max.xlsx
#> 19:      tab_assoc_RV__featNUX_drug_001_RV_gDR_x_max.xlsx
#> 20:      tab_assoc_RV__featNUX_drug_001_RV_gDR_x_max.xlsx
#> 21:      tab_assoc_RV__featNUX_drug_001_RV_gDR_x_max.xlsx
#> 22:      tab_assoc_RV__featNUX_drug_001_RV_gDR_x_max.xlsx
#> 23:      tab_assoc_RV__featNUX_drug_001_RV_gDR_x_max.xlsx
#> 24:      tab_assoc_RV__featNUX_drug_001_RV_gDR_x_max.xlsx
#> 25:     tab_assoc_RV__featNUX_drug_001_RV_gDR_x_mean.xlsx
#> 26:     tab_assoc_RV__featNUX_drug_001_RV_gDR_x_mean.xlsx
#> 27:     tab_assoc_RV__featNUX_drug_001_RV_gDR_x_mean.xlsx
#> 28:     tab_assoc_RV__metaGRP_drug_002_RV_gDR_x_mean.xlsx
#>                                                       src
#>                                                    <char>
```
