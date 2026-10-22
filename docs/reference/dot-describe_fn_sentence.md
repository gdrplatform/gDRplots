# Describe a recorded function in a sentence

Uses, in this order, the description set with
[`gDRutils::set_fit_fn_info`](https://gdrplatform.github.io/gDRstyle/reference/set_fit_fn_info.html),
the one-line source of a function written inline, or the package, name
and version of the function.

## Usage

``` r
.describe_fn_sentence(lead, fn)
```

## Arguments

- lead:

  string; the start of the sentence.

- fn:

  list; a function record from `provenance$steps`, or `NULL`.

## Value

A string, or `NULL` when `fn` is `NULL`.
