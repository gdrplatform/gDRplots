# Generate Gradient Colors for Dose Legends

Internal helper to create a named vector of colors. Each drug is
assigned a base color, and individual concentrations are assigned a
gradient intensity of that color.

## Usage

``` r
.get_gradient_colors(
  col_names,
  col_concs,
  col_labels,
  base_map,
  untreated_tag = gDRutils::get_env_identifiers("untreated_tag")
)
```

## Arguments

- col_names:

  `character` vector of drug names corresponding to the data rows.

- col_concs:

  `numeric` or `character` vector of concentrations.

- col_labels:

  `character` vector of labels to be used as names for the returned
  colors (e.g., "Drug 0.5").

- base_map:

  `character` named vector mapping drug names to their base (darkest)
  color.

- untreated_tag:

  `character` vector of identifiers for untreated/vehicle controls.
  Defaults to `gDRutils::get_env_identifiers("untreated_tag")`.

## Value

A named `character` vector of hex color codes.
