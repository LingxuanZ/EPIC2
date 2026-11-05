# EPIC2

**Integrating single-cell chromatin accessibility with human genetics to identify trait-relevant cell types**

EPIC2 is an R package for prioritizing trait-relevant cell types by integrating GWAS summary statistics with single-cell or single-nucleus ATAC-seq data. It accounts for linkage disequilibrium (LD) in four association models at SNP and regulatory-peak resolutions, using library-size-adjusted accessibility counts or cell-type-specific Wilcoxon rank-sum (WRS) test Z-scores. For each trait–cell-type pair, it combines the four nominal one-sided P-values with equal weights using the Cauchy combination test to obtain an omnibus P-value.

[Installation](#installation) · [Vignette](#vignette)  · [Methods](METHODS.md)

## Installation

EPIC2 requires R 4.3 or later. Install from GitHub using `devtools`:

```r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github(
  "LingxuanZ/EPIC2", dependencies = TRUE,
  upgrade = "never", build_vignettes = FALSE
)
```

`dependencies = TRUE` includes the analysis packages listed under `Suggests` in [DESCRIPTION](DESCRIPTION). Some dependencies require system libraries or compilers.

## Vignette

[Read the HTML vignette](https://LingxuanZ.github.io/EPIC2/EPIC2_vignette.html) · [View the Markdown version on GitHub](docs/EPIC2_vignette.md)

The vignette provides a step-by-step guide to a total cholesterol (TC) analysis using a sampled subset of real single-cell ATAC-seq data and all eligible TC SNPs. It includes downloadable example data, explains required inputs, supported file formats and cell annotations, and provides configuration guidance and code for running the workflow as a whole or stage by stage. It also includes core-function examples without YAML and describes the intermediate files, result tables and figures generated at each stage. Use the supplied `config.yml` and data preparation steps to reproduce the example.

## Support and citation

Report issues through [GitHub Issues](https://github.com/LingxuanZ/EPIC2/issues), including the package version, configuration, error message, and a minimal reproducible example. Use `citation("EPIC2")` for the software citation and record the version or commit used in your analysis.

> Zhu, L. and Jiang, Y. *EPIC2: integrating single-cell chromatin accessibility with human genetics to identify trait-relevant cell types.* [Software repository](https://github.com/LingxuanZ/EPIC2).

EPIC2 is distributed under the [MIT License](LICENSE).
