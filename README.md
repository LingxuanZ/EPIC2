# EPIC2

**LD-aware integration of GWAS and single-cell chromatin accessibility**

EPIC2 integrates GWAS summary statistics, reference genotypes, and single-cell or single-nucleus ATAC-seq data to test associations between traits and annotated cell groups. It fits four LD-adjusted models using SNP- and peak-level accessibility counts and Wilcoxon rank-sum (WRS) scores, then combines their P-values with the aggregated Cauchy association test (ACAT).

[Installation](#installation) · [Example](#example) · [Run your data](#run-your-data) · [Outputs](#outputs) · [Configuration](CONFIGURATION.md) · [Methods](METHODS.md)

## Installation

EPIC2 requires R 4.3 or later. From a terminal, install the workflow dependencies and package:

```bash
git clone https://github.com/LingxuanZ/EPIC2.git
cd EPIC2
Rscript scripts/install_dependencies.R
R CMD INSTALL .
```

The dependency script installs missing CRAN and Bioconductor packages, including Seurat, Signac, GenomicRanges, SNPRelate, and BiocParallel. It installs a missing `presto` and retains existing package installations. Some dependencies require system libraries or compilers.

Alternatively, install directly from GitHub using either `devtools` or `remotes`. Choose one of the following:

**Using devtools**

```r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github(
  "LingxuanZ/EPIC2", dependencies = NA,
  upgrade = "never", build_vignettes = FALSE
)
```

**Using remotes**

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github(
  "LingxuanZ/EPIC2", dependencies = NA,
  upgrade = "never", build_vignettes = FALSE
)
```

Both commands install required package dependencies. The repository dependency script also installs the optional packages needed for all analysis stages. [DESCRIPTION](DESCRIPTION) lists these under `Suggests`; [environment.yml](environment.yml) provides a conda environment recipe. Record installed dependency versions with your analysis.

## Example

After installation, run the synthetic example from the repository root:

```bash
Rscript scripts/run_synthetic_example.R /absolute/path/to/new-demo-directory hg19
```

The script creates ATAC, indexed fragment, GWAS, and reference-genotype inputs, runs all seven stages, and writes result tables and figures. The destination must be new or empty. Replace `hg19` with `hg38` to use that assembly. The synthetic inputs demonstrate the workflow and do not represent biological findings.

## Required inputs

| Input | Requirements |
| --- | --- |
| ATAC object | A Seurat/Signac object with a chromatin count assay, peak ranges, and grouping annotations. Defaults are assay `ATAC` and metadata columns `tissue` and `tissue_cell_type`. |
| Fragments | Accessible fragment files with tabix indexes and cell-barcode mappings, attached to the assay or supplied through configuration. SNP accessibility requires these files. |
| GWAS | Variant IDs, chromosome, position, effect and other alleles, MAF, effect estimates, standard errors, and P-values. Column names are configurable. |
| Reference genotypes | A variant-by-sample matrix or supported genotype representation with unique variant IDs matching the GWAS. Allele dosage must count the GWAS effect allele. |

All coordinates must use the declared common assembly: `hg19`/`GRCh37` or `hg38`/`GRCh38`. EPIC2 checks identifiers, assembly metadata, and coordinate bounds; inputs must be harmonized before analysis. The pipeline does not perform coordinate conversion or infer strand alignment from variant IDs.

GWAS and genotype inputs may be `.rds`, `.rda`/`.RData`, `.csv[.gz]`, `.tsv[.gz]`, or `.txt[.gz]`. Export PLINK, VCF, BGEN, or GDS data to a supported matrix or table first. Genotype encodings, allele orientation, and missing-value conventions are described in the [configuration reference](CONFIGURATION.md#gwas-and-reference-genotypes).

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

```r
config <- read_EPIC2_config("config.yml", check_files = TRUE)
run_EPIC2(config, gwas_name = "HDL")
```

Replace `HDL` with a configured dataset name. Omit `gwas_name` to process all configured datasets. The full workflow runs these stages in order:

| Stage | Purpose |
| --- | --- |
| `atac` | Filter cells and peaks; optionally compute QC metrics and embeddings. |
| `gwas` | Match reference variants, filter summary statistics and genotypes, and optionally prune LD. |
| `accessibility` | Build SNP-by-cell and GWAS-overlapping peak-by-cell matrices. |
| `ld` | Construct local LD covariance models and inverse square roots. |
| `wrs` | Compute group-specific WRS scores and optional balanced-sampling diagnostics. |
| `model` | Fit the four component models and combine their P-values with ACAT. |
| `plot` | Generate figures and a figure manifest from saved results. |

Individual stages require outputs from their prerequisites:

```r
run_EPIC2(config, stages = c("wrs", "model", "plot"), gwas_name = "HDL")
```

From the source repository, the command-line equivalents are:

```bash
Rscript exec/EPIC2 validate config.yml --check-files
Rscript exec/EPIC2 run config.yml --gwas=HDL
Rscript exec/EPIC2 stage config.yml ld --gwas=HDL
```

### 3. Inspect results

```r
component_results <- read.delim(file.path(
  config$output_dir, "05_models", "tissue_cell_type.component_results.tsv"
))
combined_results <- read.delim(file.path(
  config$output_dir, "05_models", "tissue_cell_type.cauchy_results.tsv"
))
head(combined_results[order(combined_results$p_cauchy), ])
```

ACAT combines all available valid component P-values by default. The combined table reports the number of tests used; set `model.require_all_components: true` to require all four. BH and Bonferroni adjustments are applied within each trait across tested groups, separately for each configured grouping. Review component diagnostic statuses and the number of available tests alongside the combined P-values.

## ATAC preprocessing

With `atac.compute_qc_metrics: false`, the object must already contain `TSS.enrichment`, `FRIP`, `nCount_<assay>`, `BlacklistRatio`, and `nucleosome_signal`. The ATAC stage still applies cell and peak filters. TSS and FRIP filtering precede calculation of the count-quantile cutoff.

To compute fragment-based QC, provide complete, deduplicated fragment files. For renamed cells, supply correct attached Signac mappings or set `atac.fragment_barcode_column` to a metadata column containing raw file barcodes. FRIP and blacklist ratios count each fragment record once if it overlaps any relevant interval; read multiplicity in the fifth column is ignored. See [Methods](METHODS.md#quality-control) for the definitions.

To begin with post-QC inputs without filtering again, place the prepared object at `01_atac/atac.qc.rds` under the output directory and omit the `atac` stage. Retain the full autosomal ATAC peak background for library-size normalization, including peaks that do not overlap the selected GWAS.

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

With the workflow dependencies installed, run from the source root:

```bash
R CMD build --no-build-vignettes .
EPIC2_RUN_COMMON_BUILD_TEST=true R CMD check --no-manual --no-build-vignettes EPIC2_*.tar.gz
```

## Support and citation

Report issues through [GitHub Issues](https://github.com/LingxuanZ/EPIC2/issues), including the package version, configuration, error message, and a minimal reproducible example. Use `citation("EPIC2")` for the software citation and record the version or commit used in your analysis.

> Zhu, L. and Jiang, Y. *EPIC2: LD-aware integration of single-cell chromatin accessibility and GWAS.* [Software repository](https://github.com/LingxuanZ/EPIC2).

EPIC2 is distributed under the [MIT License](LICENSE).
