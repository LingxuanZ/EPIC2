# EPIC2

**LD-aware integration of GWAS and single-cell chromatin accessibility**

EPIC2 integrates GWAS summary statistics, reference genotypes, and single-cell or single-nucleus ATAC-seq data to test associations between traits and annotated cell groups. It fits four LD-adjusted models using SNP- and peak-level accessibility counts and Wilcoxon rank-sum (WRS) scores, then combines their P-values with the aggregated Cauchy association test (ACAT).

[Installation](#installation) · [Vignette](#vignette) · [Run your data](#run-your-data) · [Outputs](#outputs) · [Configuration](CONFIGURATION.md) · [Methods](METHODS.md)

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

[Read the HTML vignette](https://LingxuanZ.github.io/EPIC2/EPIC2_vignette_knit.html) · [View the Markdown version on GitHub](docs/EPIC2_vignette_knit.md)

The vignette demonstrates a total cholesterol (TC) analysis using sampled real single-cell ATAC-seq data and all eligible TC SNPs. It provides input downloads, analysis code, result tables, and figures. Use its supplied `config_knit.yml` and data preparation steps to reproduce the example.

## Required inputs

| Input | Requirements |
| --- | --- |
| ATAC object | A Seurat/Signac object with a chromatin count assay, peak ranges, and grouping annotations. Defaults are assay `ATAC` and metadata columns `tissue` and `tissue_cell_type`. |
| Fragments | Accessible fragment files with tabix indexes and cell-barcode mappings, attached to the assay or supplied through configuration. SNP accessibility requires these files. |
| GWAS | Variant IDs, chromosome, position, effect and other alleles, MAF, effect estimates, standard errors, and P-values. Column names are configurable. |
| Reference genotypes | A variant-by-sample matrix or supported genotype representation with unique variant IDs matching the GWAS. Allele dosage must count the GWAS effect allele. |

All coordinates must use the declared common assembly: `hg19`/`GRCh37` or `hg38`/`GRCh38`. EPIC2 checks identifiers, assembly metadata, and coordinate bounds; inputs must be harmonized before analysis. The pipeline does not perform coordinate conversion or infer strand alignment from variant IDs.

GWAS and genotype inputs may be `.rds`, `.rda`/`.RData`, `.csv[.gz]`, `.tsv[.gz]`, or `.txt[.gz]`. Export PLINK, VCF, BGEN, or GDS data to a supported matrix or table first. Genotype encodings, allele orientation, and missing-value conventions are described in the [configuration reference](CONFIGURATION.md#gwas-and-reference-genotypes).

## ATAC preprocessing

- **Inputs requiring cell and peak filtering:** run the full workflow below. With `atac.compute_qc_metrics: false`, the object must already contain `TSS.enrichment`, `FRIP`, `nCount_<assay>`, `BlacklistRatio`, and `nucleosome_signal`; filtering still runs. To calculate fragment-based QC metrics, set `atac.compute_qc_metrics: true` and provide complete, deduplicated fragment files.
- **Post-QC inputs:** save the prepared object, with valid fragment attachments and cell-barcode mappings, as `01_atac/atac.qc.rds` under `project.output_dir`, then use the command below that skips `atac`. Retain the full post-QC autosomal peak background for normalization, including peaks that do not overlap the selected GWAS. The vignette demonstrates this approach.

See [ATAC configuration](CONFIGURATION.md#atac-input-and-qc) for input settings and [Methods](METHODS.md#quality-control) for QC definitions.

## Run your data

### 1. Configure the analysis

Copy the installed template:

```r
library(EPIC2)
file.copy(
  system.file("extdata", "config.example.yml", package = "EPIC2"),
  "config.yml"
)
```

Edit `config.yml` to set:

- ATAC, GWAS, and reference-genotype input paths.
- `project.output_dir` and the required `genome.build`.
- GWAS column mappings, genotype encoding, grouping columns, and QC settings.
- The datasets to analyze under `gwas.traits`; remove unused template entries.

Paths are resolved relative to the configuration file. The [configuration reference](CONFIGURATION.md) documents all analysis settings.

### 2. Validate and run

For inputs requiring cell and peak filtering:

```r
config <- read_EPIC2_config("config.yml", check_files = TRUE)
run_EPIC2(config, gwas_name = "HDL")
```

For a post-QC checkpoint prepared as described above, use this instead:

```r
config <- read_EPIC2_config("config.yml", check_files = TRUE)
run_EPIC2(
  config,
  stages = c("gwas", "accessibility", "ld", "wrs", "model", "plot"),
  gwas_name = "HDL"
)
```

Replace `HDL` with a dataset name in your configuration. Omit `gwas_name` to process all configured datasets.

The full workflow runs `atac → gwas → accessibility → ld → wrs → model → plot`. Both commands include figure generation.

### 3. Inspect results

```r
combined_results <- read.delim(file.path(
  config$output_dir, "05_models", "tissue_cell_type.cauchy_results.tsv"
))
head(combined_results[order(combined_results$p_cauchy), ])
```

ACAT combines all available valid component P-values by default. The combined table reports the number of tests used; set `model.require_all_components: true` to require all four. BH and Bonferroni adjustments are applied within each trait across tested groups, separately for each configured grouping. Review component diagnostic statuses and the number of available tests alongside the combined P-values.

## Outputs

All outputs are written beneath `project.output_dir`:

```text
results/
├── 01_atac/             # Filtered ATAC object and QC summaries
├── 02_gwas/             # Matched GWAS, genotypes, and QC records
├── 02_accessibility/    # SNP/peak matrices and overlap tables
├── 03_ld/               # LD covariance models and inverse square roots
├── 04_wrs/              # WRS scores and optional stability summaries
├── 05_models/           # Component and ACAT results
├── figures/             # Figures, supporting tables, figure_manifest.tsv
├── run_manifest.txt
├── config.used.yml
├── config.used.rds
└── run_history/         # Per-invocation settings and session manifests
```

| Main file | Contents |
| --- | --- |
| `05_models/<grouping>.component_results.tsv` | Focal coefficient estimates, standard errors, P-values, and model statuses. |
| `05_models/<grouping>.cauchy_results.tsv` | Component P-values, `p_cauchy`, number of tests combined, and adjusted P-values. |
| `figures/figure_manifest.tsv` | Figure paths and `generated`, `skipped`, or `error` status. |
| `config.used.yml` / `config.used.rds` | Effective configuration, including in-memory overrides. |
| `run_manifest.txt` | Software environment, selected stages and traits, assembly, seed, and worker settings. |

The default groupings are `tissue_cell_type` and `tissue`. S-LDSC summaries can optionally be imported for comparison figures through `plots.sldsc_summary`; S-LDSC is not required for EPIC2 analysis.

## Execution and performance

Use a separate output directory for each independent analysis. Requested stages rerun explicitly. Sequential model runs preserve completed results for other traits still listed in the configuration; completion records verify saved artifacts before aggregation. Concurrent jobs require separate output directories. [Methods](METHODS.md#saved-results-and-reproducibility) describes these checks and recovery after an interrupted run.

`parameters.workers` controls SNP-accessibility workers. LD correlation batching is controlled by `ld.correlation_block_size`, while BLAS threading depends on the R environment. SNPRelate pruning uses `gwas.ld_num_threads`. For an LD-only batch job, [run_ld_stage.sh](scripts/run_ld_stage.sh) sets thread limits before starting R.

Dense LD blocks are limited to 5,000 features by default. A local distance window can still produce a large connected component, so this guard is not a total-memory limit. See [LD configuration](CONFIGURATION.md#accessibility-ld-and-models) before increasing it. The [correlation benchmark](scripts/benchmark_hotspots.R) and [dense-inverse benchmark](scripts/benchmark_dense_inverse.R) measure selected operations on your machine.

## Testing

Tests cover fragment counting, genotype normalization, GWAS filtering, covariance calculations, transformed regression, WRS, ACAT, and saved-result validation. The optional integration test runs both supported assemblies using synthetic inputs.

Automated checks run through [GitHub Actions](https://github.com/LingxuanZ/EPIC2/actions/workflows/R-CMD-check.yaml). For a worked analysis, see the [vignette](https://LingxuanZ.github.io/EPIC2/EPIC2_vignette_knit.html).

## Support and citation

Report issues through [GitHub Issues](https://github.com/LingxuanZ/EPIC2/issues), including the package version, configuration, error message, and a minimal reproducible example. Use `citation("EPIC2")` for the software citation and record the version or commit used in your analysis.

> Zhu, L. and Jiang, Y. *EPIC2: integrating single-cell chromatin accessibility with human genetics to identify trait-relevant cell types.* [Software repository](https://github.com/LingxuanZ/EPIC2).

EPIC2 is distributed under the [MIT License](LICENSE).
