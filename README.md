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

[Read the HTML vignette](https://LingxuanZ.github.io/EPIC2/EPIC2_vignette.html) · [View the Markdown version on GitHub](docs/EPIC2_vignette.md)

The vignette demonstrates a total cholesterol (TC) analysis using sampled real single-cell ATAC-seq data and all eligible TC SNPs. It provides input downloads, analysis code, result tables, and figures. Use its supplied `config.yml` and data preparation steps to reproduce the example.

## Required inputs

| Input | Requirements |
| --- | --- |
| ATAC object | A Seurat/Signac object with a chromatin assay (`ATAC` by default), full peak-by-cell counts in the `counts` layer, genomic peak ranges, and per-cell grouping annotations. |
| Fragments | Accessible fragment files with tabix indexes and cell-barcode mappings, attached to the assay or supplied through configuration. SNP accessibility requires these files. |
| GWAS | Variant IDs, chromosome, position, effect and other alleles, MAF, effect estimates, standard errors, and P-values. Column names are configurable. |
| Reference genotypes | A variant-by-sample matrix or supported genotype representation with unique variant IDs matching the GWAS. Allele dosage must count the GWAS effect allele. |
| Configuration | A YAML file defining input paths and analysis settings; see [Configure the analysis](#1-configure-the-analysis). |

Cell annotations belong in `atac[[]]`, with one row per cell and row names matching `colnames(atac)`. Every column listed in `model.groupings` must exist; the default columns are `tissue_cell_type` and `tissue`. Each grouping needs at least two labelled groups; cells with missing labels are excluded from that grouping's analysis.

For cell-type-only analysis, one annotation column is sufficient if it is the only entry in `model.groupings`; a separate tissue column is not required for those tests. If running ATAC QC, also set `atac.group_column` to an available annotation column. See the [vignette's annotation examples](https://LingxuanZ.github.io/EPIC2/EPIC2_vignette.html#choose-annotation-levels-for-your-own-data) and [input requirements](https://LingxuanZ.github.io/EPIC2/EPIC2_vignette.html#example-data-and-required-inputs), including metadata for tissue-mapped fragment attachment and raw-barcode mapping.

All coordinates must use the declared common assembly: `hg19`/`GRCh37` or `hg38`/`GRCh38`. EPIC2 checks identifiers, assembly metadata, and coordinate bounds; inputs must be harmonized before analysis. The pipeline does not perform coordinate conversion or infer strand alignment from variant IDs.

GWAS and genotype inputs may be `.rds`, `.rda`/`.RData`, `.csv[.gz]`, `.tsv[.gz]`, or `.txt[.gz]`. Export PLINK, VCF, BGEN, or GDS data to a supported matrix or table first. Genotype encodings, allele orientation, and missing-value conventions are described in the [configuration reference](CONFIGURATION.md#gwas-and-reference-genotypes).

## ATAC preprocessing

- **Inputs requiring cell and peak filtering:** run the full workflow below. With `atac.compute_qc_metrics: false`, the object must already contain `TSS.enrichment`, `FRIP`, `nCount_<assay>`, `BlacklistRatio`, and `nucleosome_signal`; filtering still runs. To calculate fragment-based QC metrics, set `atac.compute_qc_metrics: true` and provide complete, deduplicated fragment files.
- **Post-QC inputs:** save the prepared object, with valid fragment attachments and cell-barcode mappings, as `01_atac/atac.qc.rds` under `project.output_dir`, then use the command below that skips `atac`. Retain the full post-QC autosomal peak background for normalization, including peaks that do not overlap the selected GWAS. QC metrics are not required for downstream calculations when `atac` is skipped; retain them for QC figures. The [vignette](https://LingxuanZ.github.io/EPIC2/EPIC2_vignette.html#prepare-the-shared-inputs) demonstrates checkpoint preparation.

UMAP is optional for statistical analysis. Retain the `umap` reduction in the Seurat object to reproduce UMAP figures.

See [ATAC configuration](CONFIGURATION.md#atac-input-and-qc) for input settings and [Methods](METHODS.md#quality-control) for QC definitions.

## Run your data

### 1. Configure the analysis

Copy the installed [full-workflow template](inst/extdata/config.example.yml), which covers all stages from ATAC preprocessing through figure generation:

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

The [vignette](https://LingxuanZ.github.io/EPIC2/EPIC2_vignette.html#configuration-and-cell-type-annotations) uses a separate [config.yml](docs/config.yml) for a QC-filtered subset of real scATAC-seq data and skips ATAC preprocessing. It also disables repeated GWAS LD pruning and ambiguous-variant removal because those steps have already been applied.

### 2. Validate and run

The commands below show how to run the overall workflow. For details of each stage, including required inputs, code, and outputs, see the [vignette](https://LingxuanZ.github.io/EPIC2/EPIC2_vignette.html#run-the-analysis).

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

The core result data frames are saved as TSV tables in `05_models/` under `project.output_dir`: `<grouping>.cauchy_results.tsv` contains the combined ACAT results, and `<grouping>.component_results.tsv` contains the individual model results. For example:

```r
combined_results <- read.delim(file.path(
  config$output_dir, "05_models", "tissue_cell_type.cauchy_results.tsv"
))
head(combined_results[order(combined_results$p_cauchy), ])
```

ACAT combines all available valid component P-values by default. The combined table reports the number of tests used; set `model.require_all_components: true` to require all four. BH and Bonferroni adjustments are applied within each trait across tested groups, separately for each configured grouping. Review component diagnostic statuses and the number of available tests alongside the combined P-values.

See the [vignette's association results](https://LingxuanZ.github.io/EPIC2/EPIC2_vignette.html#association-models-and-acat) for example result tables. Figures are saved under `figures/` in PNG and PDF formats by default; see the [vignette's figure section](https://LingxuanZ.github.io/EPIC2/EPIC2_vignette.html#generate-figures) for the main figures, plotting functions, filenames, and output formats.

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

Automated checks run through [GitHub Actions](https://github.com/LingxuanZ/EPIC2/actions/workflows/R-CMD-check.yaml). For a worked analysis, see the [vignette](https://LingxuanZ.github.io/EPIC2/EPIC2_vignette.html).

## Support and citation

Report issues through [GitHub Issues](https://github.com/LingxuanZ/EPIC2/issues), including the package version, configuration, error message, and a minimal reproducible example. Use `citation("EPIC2")` for the software citation and record the version or commit used in your analysis.

> Zhu, L. and Jiang, Y. *EPIC2: integrating single-cell chromatin accessibility with human genetics to identify trait-relevant cell types.* [Software repository](https://github.com/LingxuanZ/EPIC2).

EPIC2 is distributed under the [MIT License](LICENSE).
