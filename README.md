# EPIC2

**LD-aware integration of GWAS and single-cell chromatin accessibility**

EPIC2 is an R package that integrates GWAS summary statistics, a reference genotype panel, and single-cell or single-nucleus ATAC-seq data. It constructs SNP- and peak-level accessibility features, models local linkage disequilibrium (LD), and evaluates trait associations with annotated cell groups using four component models and the aggregated Cauchy association test (ACAT).

EPIC2 provides a configurable workflow for ATAC quality control, GWAS/reference matching, accessibility matrix construction, LD calculation, Wilcoxon rank-sum (WRS) scoring, association modeling, and visualization. Individual stages can also be run separately when their required inputs have already been prepared.

[Installation](#installation) · [Quick start](#quick-start) · [Required inputs](#required-inputs) · [Tutorial](#tutorial) · [Outputs](#outputs) · [Configuration reference](https://github.com/LingxuanZ/EPIC2/blob/main/CONFIGURATION.md)

## Workflow

```text
Input: ATAC-seq object + GWAS summary statistics + reference genotypes
                              |
                     ATAC and GWAS QC
                              |
              SNP- and peak-level accessibility
                              |
                  Local LD and WRS scores
                              |
              Four LD-adjusted association models
                              |
                  ACAT results and figures
```

| Stage           | Purpose                                                      |
| --------------- | ------------------------------------------------------------ |
| `atac`          | Filter ATAC cells and peaks; optionally compute QC metrics and embeddings. |
| `gwas`          | Match variants to reference genotypes, apply configured GWAS/reference QC, and optionally perform LD pruning. |
| `accessibility` | Construct SNP-by-cell and GWAS-overlapping peak-by-cell matrices. |
| `ld`            | Construct local SNP- and peak-level LD models and inverse-square-root covariance representations. |
| `wrs`           | Calculate group-specific WRS scores, with optional balanced subsampling and stability summaries. |
| `model`         | Fit SNP accessibility, SNP WRS, peak accessibility, and peak WRS models; combine their P-values with ACAT. |
| `plot`          | Generate figures and a figure manifest from available analysis outputs. |

## Installation

EPIC2 requires **R 4.3 or later**. The full workflow also uses CRAN and Bioconductor packages, including Seurat, Signac, GenomicRanges, SNPRelate, BiocParallel, and `presto`.

### Install the full workflow

Clone the repository, then run the following commands from its root directory, which contains `DESCRIPTION`:

```bash
git clone https://github.com/LingxuanZ/EPIC2.git
cd EPIC2
Rscript scripts/install_dependencies.R
R CMD INSTALL .
```

Alternatively, download and extract the repository ZIP, open a terminal in the extracted package directory, and run the last two commands.

The dependency script installs missing CRAN/Bioconductor packages and installs `presto` from its upstream GitHub repository. It also installs `knitr` and `rmarkdown` for tutorial documents. Additional system libraries or compilers may be required by these dependencies on your operating system.

### Install the package directly from GitHub

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  "LingxuanZ/EPIC2",
  dependencies = NA,
  upgrade = "never",
  build_vignettes = FALSE
)
```

This installs EPIC2 and its required package dependencies. Additional packages used by particular analysis stages are listed under `Suggests` in [DESCRIPTION](https://github.com/LingxuanZ/EPIC2/blob/main/DESCRIPTION); use the repository dependency script to prepare the full workflow.

## Quick start

This example runs one configured GWAS dataset using your own inputs. The example YAML is a configuration template; it does not download or prepare analysis data.

### 1. Create a configuration file

```r
library(EPIC2)

example_config <- system.file(
  "extdata", "config.example.yml", package = "EPIC2"
)
file.copy(example_config, "config.yml")
```

Edit `config.yml` before continuing:

- Set the ATAC object, GWAS summary statistics, and reference genotype paths.
- Choose an output directory for this analysis.
- Set the required common `genome.build` and confirm that the GWAS alleles match the reference panel.
- Set the GWAS column mapping, genotype encoding, grouping columns, and QC options for your data.
- Keep only the entries you intend to configure under `gwas.traits`; the installed template includes more than one example name.

Relative input and output paths are resolved relative to `config.yml`.

### 2. Validate and run one GWAS dataset

```r
config <- read_EPIC2_config("config.yml", check_files = TRUE)
run_EPIC2(config, gwas_name = "HDL")
```

Replace `HDL` with a name defined under `gwas.traits`. The call runs all seven stages for the selected dataset. Configuration validation checks settings and configured file paths; it does not establish allele alignment or prove that a declared genome assembly is correct.

### 3. Inspect the main results

With `project.output_dir: results` and the default `tissue_cell_type` grouping:

```r
results_dir <- config$output_dir

component_results <- read.delim(file.path(
  results_dir, "05_models", "tissue_cell_type.component_results.tsv"
))
combined_results <- read.delim(file.path(
  results_dir, "05_models", "tissue_cell_type.cauchy_results.tsv"
))

head(combined_results[order(combined_results$p_cauchy), ])
```

The combined table contains the component P-values, `p_cauchy`, the number of tests combined, and within-trait BH and Bonferroni adjustments. By default, ACAT requires all four component P-values to be valid. Adjustments are within each trait across tested groups; they are not a single correction across all traits.

### Command-line usage

From the source repository root, after installing EPIC2:

```bash
Rscript exec/EPIC2 validate config.yml --check-files
Rscript exec/EPIC2 run config.yml --gwas=HDL
```

To run a single stage after its prerequisite outputs exist:

```bash
Rscript exec/EPIC2 stage config.yml ld --gwas=HDL
Rscript exec/EPIC2 run config.yml --stages=wrs,model,plot --gwas=HDL
```

The relative `exec/EPIC2` path refers to a file in the source checkout. For R sessions outside the checkout, use the exported R functions shown above.

## Required inputs

| Input                   | Requirements                                                 |
| ----------------------- | ------------------------------------------------------------ |
| ATAC object             | A Seurat/Signac object with a chromatin count assay, genomic peak ranges, and annotation columns used for grouping. The default assay is `ATAC`; default metadata names are `tissue` and `tissue_cell_type`. |
| Fragment files          | Fragment objects attached to the chromatin assay, or configured fragment files that can be attached. SNP coverage requires accessible fragment files and their tabix indexes, with barcode mappings matching the ATAC cells. |
| GWAS summary statistics | Variant IDs, chromosome, position, effect and other alleles, MAF, effect estimate, standard error, and P-value columns. Column names are configurable. Alleles must already be aligned to the reference panel. |
| Reference genotypes     | A user-supplied variant-by-sample genotype/dosage matrix, or another supported encoding, with unique variant identifiers that match the GWAS data. A specific reference panel or hg19-only panel is not required. |

GWAS and genotype files can be supplied as `.rds`, `.rda`/`.RData`, `.csv[.gz]`, `.tsv[.gz]`, or `.txt[.gz]`. PLINK, VCF, BGEN, and GDS containers must first be read or exported to a supported matrix/table representation.

Genotype normalization supports hard calls, dosage, GT strings, biallelic allele labels, and genotype probability representations. Set `genotype_format` and `genotype_orientation` explicitly when automatic detection is ambiguous. Carrier-only presence/absence data do not contain enough information for exact MAF calculation and are rejected. See the [configuration reference](https://github.com/LingxuanZ/EPIC2/blob/main/CONFIGURATION.md) for encoding and ploidy requirements.

Variants are matched by identifier and reordered to the GWAS rows; variants genuinely absent from the reference are skipped. Genotype variant IDs must be supplied as unique row names (column names for transposed input) or through `genotype_variant_id_column`. Matching dimensions or setting `genotype_orientation` does not establish correspondence, and GWAS IDs are never assigned automatically to unnamed genotype rows. Coordinates used for matching ATAC features and defining LD windows come from the GWAS table; genotype dosage values themselves are not lifted between assemblies. Matching identifiers alone does not verify allele coding.

### ATAC QC and preprocessed inputs

When `atac.compute_qc_metrics: false`, the input object must contain the QC metadata required by the configured filters, including `TSS.enrichment`, `FRIP`, `nCount_ATAC`, `BlacklistRatio`, and `nucleosome_signal` under the default settings.

**Disabling QC metric calculation does not disable ATAC filtering.** Running the `atac` stage applies the configured cell and peak filters, even when the object has already been filtered elsewhere.

To start from post-QC data without filtering again, prepare the expected `01_atac/atac.qc.rds` stage input in the analysis output directory before omitting the `atac` stage. Keep the full retained ATAC peak background for the selected cells: library-size normalization uses the autosomal ATAC counts, not only the peaks overlapping a particular GWAS.

### Genome builds

The analysis requires one explicitly declared common assembly. Set the required
`genome.build` to `hg19`/`GRCh37` or `hg38`/`GRCh38`; there is no default genome
build. GWAS, ATAC peaks, and fragment coordinates must already be harmonized.

```yaml
genome:
  build: hg38
```

Legacy `genome.target`, `atac.genome`, `gwas.genome`, and per-dataset `genome`
fields may be omitted. Any explicitly supplied value must match `genome.build`.
Nonempty pipeline `genome.chains` settings are rejected. The analysis pipeline
validates the current coordinates and does not call liftOver. The standalone
`lift_genomic_ranges()` utility remains available for separate preprocessing.

Available build metadata and coordinate bounds are checked, but cannot prove
the true assembly or replace allele alignment. See the [common assembly
configuration reference](https://github.com/LingxuanZ/EPIC2/blob/main/CONFIGURATION.md#common-genome-assembly).

## Tutorial

A small-data vignette is being prepared as an R Markdown document with an HTML report. It walks through GWAS/reference matching, SNP and peak accessibility, LD, WRS, association models, and result inspection using post-QC ATAC inputs. The default example uses one trait, HDL.

**Current status: draft.** The draft defaults to inspecting available outputs without executing the analysis (`run_analysis: false`). Its prepared inputs still depend on the original HPC fragment files; a standalone public example dataset and hosted HTML tutorial are not yet provided here.

The tutorial begins after ATAC QC and therefore omits the `atac` stage. The full workflow in the quick start includes that stage. Tutorial results from a randomly selected subset are intended for demonstrating the software and do not replace a full-data analysis.

## Outputs

All pipeline outputs are written beneath `project.output_dir`.

```text
results/
├── 01_atac/             # ATAC object, QC summaries, group counts
├── 02_gwas/             # Per-dataset GWAS, genotypes, QC records
├── 02_accessibility/    # SNP/peak matrices and overlap tables
├── 03_ld/               # Local LD models and covariance representations
├── 04_wrs/              # WRS scores and optional stability summaries
├── 05_models/           # Component and ACAT result tables
├── figures/             # Figures, plotting tables, figure_manifest.tsv
├── run_manifest.txt
└── config.used.yml      # Written when a YAML configuration file is available
```

| Main output                                  | Contents                                                     |
| -------------------------------------------- | ------------------------------------------------------------ |
| `05_models/<grouping>.component_results.tsv` | Group-specific estimates and P-values for the four component models. |
| `05_models/<grouping>.cauchy_results.tsv`    | Component P-values, ACAT P-values, and within-trait adjusted P-values. |
| `figures/figure_manifest.tsv`                | Figure paths and `generated`, `skipped`, or `error` status.  |
| `run_manifest.txt`                           | Software/environment information, selected stages and traits, genome declarations, seed, and worker setting for the invocation. |

Default grouping names include `tissue_cell_type` and `tissue`; configure `model.groupings` for your annotations. The plotting stage can also import an externally generated S-LDSC summary when `plots.sldsc_summary` is provided.

Review model statistics and figure statuses alongside the main P-value tables. Small or uninformative feature sets can yield unavailable component tests; a successfully completed command does not guarantee that every model or optional figure has a valid result.

## Configuration and performance

See the [example configuration](https://github.com/LingxuanZ/EPIC2/blob/main/inst/extdata/config.example.yml) and [complete configuration reference](https://github.com/LingxuanZ/EPIC2/blob/main/CONFIGURATION.md) for QC thresholds, column mappings, genotype formats, LD settings, balanced WRS, and plotting options.

- If GWAS variants have already been LD-pruned, provide matching inputs and set `gwas.ld_prune: false` to avoid repeating pruning.
- `parameters.workers` controls workers for SNP accessibility feature chunks. It does not automatically run different traits in parallel.
- Local LD calculations use configurable blocks (`ld.correlation_block_size`); LD loops are serial, with BLAS threading determined by the R environment. SNPRelate pruning has a separate `gwas.ld_num_threads` setting.
- Balanced WRS repeats are optional and can substantially increase work. Enable them when the analysis requires these diagnostics.
- Stage selection explicitly reruns the requested stages. It is not an automatic resume or cache-validation mechanism.

For concurrent per-trait jobs, use a separate output directory for each independent run. Model tables and figure manifests use shared filenames within an output directory and must not be written concurrently by separate trait jobs. A single run with multiple selected traits produces combined model tables for those traits.

Runtime and memory depend on the number of cells and features, reference sample size, LD connectivity, WRS repeats, and available resources. The repository includes [benchmark_hotspots.R](https://github.com/LingxuanZ/EPIC2/blob/main/scripts/benchmark_hotspots.R) for comparing numerical results and timings of selected numerical routines on your system.

## Reproducibility and testing

Record the installed EPIC2 version, repository commit or release, input provenance, configuration, and R session information for each analysis. Keep separate output directories when comparing configurations, and retain manifests from runs that need to be compared.

```r
packageVersion("EPIC2")
sessionInfo()
```

The repository contains [core tests](https://github.com/LingxuanZ/EPIC2/blob/main/tests/test-core.R), [LD comparison tests](https://github.com/LingxuanZ/EPIC2/blob/main/tests/test-ld-rmd.R), common-build configuration and accessibility tests, and an optional [same-build integration test](https://github.com/LingxuanZ/EPIC2/blob/main/tests/test-common-build-integration.R). The integration fixture tests hg19 and hg38 separately and fails on accidental pipeline liftOver calls; the [mixed-build test](https://github.com/LingxuanZ/EPIC2/blob/main/tests/test-mixed-build-integration.R) verifies early rejection. Full integration requires the complete workflow dependencies and `EPIC2_RUN_COMMON_BUILD_TEST=true`. Configuration validation alone is not an end-to-end analysis test.

## Support

Please use [GitHub Issues](https://github.com/LingxuanZ/EPIC2/issues) for questions, bug reports, and feature requests. For a reproducible report, include the EPIC2 version and commit, relevant configuration settings, error messages, and a minimal example or input dimensions, while excluding private data and credentials.

## Citation

If you use EPIC2, cite the software repository and report the version or commit used:

> Zhu, L. and Jiang, Y. *EPIC2: LD-aware integration of single-cell chromatin accessibility and GWAS.* [Software repository](https://github.com/LingxuanZ/EPIC2).

## License

EPIC2 is distributed under the [MIT License](https://github.com/LingxuanZ/EPIC2/blob/main/LICENSE.md).
