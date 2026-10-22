# Describe one recorded processing step in a sentence

Steps recorded by gDRimport and gDRcore (`"import"`, `"create"`,
`"normalize"`, `"average"`, `"growth_rates"`) have their own
description; any other step lists its settings.

## Usage

``` r
.describe_provenance_step(step, entry)
```

## Arguments

- step:

  string; name of the step.

- entry:

  list; the step entry from `provenance$steps`, with `settings` and, for
  some steps, `functions`.

## Value

A string.
