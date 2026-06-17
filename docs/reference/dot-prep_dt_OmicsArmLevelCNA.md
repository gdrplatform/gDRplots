# Encode OmicsArmLevelCNA as not mutated and mutated

OmicsArmLevelCNA is arm-level copy number alteration inferred using
absolute copy number data from PureCN, method from the Ben-David et al.
(2021) paper (https://www.nature.com/articles/s41586-020-03114-6).
Chromosome arms: *1* indicates arm-level gain, *-1* indicates arm-level
loss, and *0* indicates copy-neutral.

## Usage

``` r
.prep_dt_OmicsArmLevelCNA(dt_depmap)
```

## Arguments

- dt_depmap:

  `data.table` with dependent variables data load from DepMap. (rows are
  samples, columns are features or meta); outputted by one of
  [`prep_dt_depmap_feat`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_depmap_feat.md)
  for OmicsArmLevelCNA

## Value

`data.table` with OmicsArmLevelCNA decoded as mutated - not mutated

## Details

This function transform each chromosome column (e.g., `3p`) into two new
binary columns: `3p_loss` and `3p_gain`. `3p_loss` is *1* for values of
*-1* in the original column and *0* otherwise. `3p_gain` is *1* for
values of *1* in the original column and *0* otherwise. The original
chromosome column is then removed.

## Author

Janina Smoła <janina.smola@contractors.roche.com>
