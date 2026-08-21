# Prep table with calculated linear associations

Prep table with calculated linear associations

## Usage

``` r
prep_dt_assoc(dt_response, dt_depmap, selected_feat_meta_col = NULL)
```

## Arguments

- dt_response:

  `data.table` with experimental response data (rows are samples) for
  one metric

- dt_depmap:

  `data.table` with dependent variables data load from DepMap. (rows are
  samples, columns are features or meta); outputted by one of
  [`prep_dt_depmap_feat`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_depmap_feat.md)
  or
  [`prep_dt_depmap_meta`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_depmap_meta.md)

- selected_feat_meta_col:

  string name of feature/meta column in DepMap

## Value

A named list with elements, that may be input to
[`plot_volcano_assoc`](https://gdrplatform.github.io/gDRplots/reference/plot_volcano_assoc.md)

- `dt_assoc` `data.table` with calculated association values between
  feature/meta of DepMap and selected metric,

- `condition_info` string describing experiment condition (drugs),

- `selected_feat_meta_col` string name of feature/meta.

## Author

Janina Smoła <janina.smola@external.roche.com>

## Examples

``` r
mae_prism <- gDRutils::get_synthetic_data("prism")
se_prism <- mae_prism[[gDRutils::get_supported_experiments("sa")]]
dt_metrics <- gDRutils::convert_se_assay_to_dt(se_prism, "Metrics")
dt_response <- prep_dt_response_metric_sa(
  dt_metrics = dt_metrics,
  d_name = unique(dt_metrics[["DrugName"]])[1],
  normalization_type = "RV", metric = "x_mean")
feat_data_path <- system.file("depmap_data", package = "gDRtestData")
meta_data_path <- system.file("depmap_data/Model.csv.gz", package = "gDRtestData")
obj_depmap_feat <- prep_dt_depmap_feat(feat_data_path, meta_data_path,
                                       feature_set = "OmicsCNGene")
prep_dt_assoc(
  dt_response = dt_response,
  dt_depmap = obj_depmap_feat[["dt_depmap"]],
  selected_feat_meta_col = obj_depmap_feat[["selected_feat_meta_col"]])
#> $dt_assoc
#>              feature      response           rho     q_value
#>               <char>        <char>         <num>       <num>
#>   1:      ABCA1 (19) RV_gDR_x_mean -0.0005082346 0.667729711
#>   2:   ABCC9 (10060) RV_gDR_x_mean  0.0748785160 0.386959685
#>   3:      ACACB (32) RV_gDR_x_mean  0.0506714309 0.440011930
#>   4:  ACTR10 (55860) RV_gDR_x_mean  0.0215449751 0.700357070
#>   5:     ACVR1B (91) RV_gDR_x_mean  0.0373464528 0.501222920
#>  ---                                                        
#> 146: ZNF534 (147658) RV_gDR_x_mean -0.1482499643 0.002867780
#> 147: ZNF615 (284370) RV_gDR_x_mean -0.1424800139 0.003925022
#> 148:     ZNF7 (7553) RV_gDR_x_mean  0.0390147465 0.722562541
#> 149: ZNF732 (654254) RV_gDR_x_mean -0.0230505747 0.683037687
#> 150:  ZSCAN20 (7579) RV_gDR_x_mean -0.0645321004 0.250615186
#> 
#> $condition_info
#> [1] "G00106_GDC-8025_TEAD1|TEAD2|TEAD3|TEAD4_120"
#> 
#> $selected_metric
#> [1] "RV_gDR_x_mean"
#> 
#> $selected_feat_meta_col
#> [1] "OmicsCNGene"
#> 
```
