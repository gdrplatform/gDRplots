# Binarize somatic mutations in OmicsSomaticMutationsMatrixHotspot and OmicsSomaticMutationsMatrixDamaging

OmicsSomaticMutationsMatrixHotspot is genotyped matrix determining for
each cell line whether each gene has at least one hot spot mutation. A
variant is considered a hot spot if it's present in one of the
following: Hess et al. (2019) paper, OncoKB hotspot, COSMIC mutation
significance tier 1.

## Usage

``` r
.prep_dt_OmicsSomaticMutationsMatrix(dt_depmap)
```

## Arguments

- dt_depmap:

  `data.table` with dependent variables data load from DepMap. (rows are
  samples, columns are features or meta); outputted by one of
  [`prep_dt_depmap_feat`](https://gdrplatform.github.io/gDRplots/reference/prep_dt_depmap_feat.md)
  for OmicsSomaticMutationsMatrixHotspot or
  OmicsSomaticMutationsMatrixDamaging

## Value

`data.table` with OmicsSomaticMutationsMatrixHotspot or
OmicsSomaticMutationsMatrixDamaging decoded as not mutated and mutated

## Details

OmicsSomaticMutationsMatrixDamaging is genotyped matrix determining for
each cell line whether each gene has at least one damaging mutation. A
variant is considered a damaging mutation if LikelyLoF is True

*0* means no mutation; if there is one or more hot spot mutations or
damaging mutations respectively, in the same gene for the same cell
line, the allele frequencies are summed, and if the sum is greater than
0.95, a value of *2* is assigned (representing a likely homozygous
mutation), otherwise a value of *1* is assigned (likely heterozygous).

This function transforms each gene column into binary columns: *0*
indicates no mutation, and *1* indicates mutation, regardless of
zygosity (the original value was *1* or *2*).

## Author

Janina Smoła <janina.smola@contractors.roche.com>
