# gDRplots

The `gDRplots` package is intended to store functions for gDR static
visualizations and their dedicated helpers.

# Website

A package website is available under [this
link](https://gdrplatform.github.io/gDRplots/).

# AI use

gDRplots has been developed by the gDR team since May 2024. Since April
2026, parts of the code, tests and documentation have been written with
the assistance of Claude Code (Anthropic) under the maintainers’
direction. Plot design, package scope and methods are the team’s; every
change is reviewed by another team member and covered by the test suite.
Third-party code: the internal
[`.lin_associations()`](https://gdrplatform.github.io/gDRplots/reference/dot-lin_associations.md)
used by
[`calc_assoc()`](https://gdrplatform.github.io/gDRplots/reference/calc_assoc.md)
is adapted from `lin_associations()` in
[cdsrmodels](https://github.com/broadinstitute/cdsr_models) (MIT
License; authors William Colgan and Andrew Boghossian) to avoid a
GitHub-only dependency; the license text and copyright notice are in
[`inst/licenses/cdsrmodels-MIT.txt`](https://gdrplatform.github.io/gDRplots/inst/licenses/cdsrmodels-MIT.txt).
