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
No third-party code is included;
[`calc_assoc()`](https://gdrplatform.github.io/gDRplots/reference/calc_assoc.md)
reimplements
[`cdsrmodels::lin_associations`](https://rdrr.io/pkg/cdsrmodels/man/lin_associations.html)
to avoid a GitHub-only dependency.
