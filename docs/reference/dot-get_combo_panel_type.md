# Check type of concentrations set for combination of drug per cell line

This function checks if the concentration vectors per cell line have
common part or not. It is required to decide which function to use for
plotting the heatmaps panel (the set of heatmaps with combo metrics plot
for combination of selected drug with selected codrug, each heatmap per
one cell line):
[`heatmap_combo_with_isoref_panel_common`](heatmap_combo_with_isoref_panel_common.md)
that plot heatmaps with shared axes or
[`heatmap_combo_with_isoref_panel_independent`](heatmap_combo_with_isoref_panel_independent.md)
that plot heatmaps independently.

## Usage

``` r
.get_combo_panel_type(ls_vec_conc)
```

## Arguments

- ls_vec_conc:

  a list with vectors with concentration per cell line

## Value

a string decribing type of concentration list - one of:

- `common` vectors in the input list have common part; heatmaps in the
  panel should be created jointly with function
  [`heatmap_combo_with_isoref_panel_common`](heatmap_combo_with_isoref_panel_common.md)

- `independent` vectors in the input list do not have common part;
  heatmaps in the panel should be created independently with function
  [`heatmap_combo_with_isoref_panel_independent`](heatmap_combo_with_isoref_panel_independent.md)

## Details

Possible combinations (0 is not taken into account as it is always
present):

- `common`

  - all vectors have common part and start and end conc are the same  
    0, 0.03, 0.1, 0.3, 1  

  - all vectors have common part and start and end conc are the same,
    but there are some gap  
    0, 0.003, 0.01, 0.03, 0.1, 0.3, 1  
    0, 0.003, 0.01, 0.03, \_\_\_, 0.3, 1  

- `independent`

  - all vectors have common part and start conc are the same, but end
    conc is different (no gap)  
    0, 0.003, 0.01, 0.03, 0.1, 0.3, \_\_  
    0, 0.003, 0.01, 0.03, 0.1, 0.3, 1  

  - no common part between vectors  
    0, 0.003, 0.01, 0.03, \_\_, \_\_, \_\_  
    0, \_\_\_\_, \_\_\_\_, \_\_\_\_, 0.1, 0.3, 1  

  - all vectors have common part but start and end conc are different
    (shifted range)  
    0, 0.003, 0.01, 0.03, 0.1, 0.3, 1, \_\_, \_\_  
    0, \_\_\_\_, \_\_\_\_, 0.03, 0.1, 0.3, 1, 3, 10  

## Author

Janina Smoła <janina.smola@contractors.roche.com>

## Examples

``` r
if (FALSE) { # \dontrun{
ls_conc <- list(c(0, 0.003, 0.01, 0.03), c(0, 0.003, 0.01, 0.03, 0.1))
.get_combo_panel_type(ls_conc)
} # }
```
