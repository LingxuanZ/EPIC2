# EPIC2

`EPIC2` is the LD-aware workflow described in`Fetal_Tissue.Rmd` before section 6. It integrates aligned GWAS summary statistics, a reference genotype panel, and a Seurat/Signac single-nucleus ATAC-seq object at SNP and peak resolution.

The package implements these manuscript stages:

1. ATAC cell and peak quality control.
2. GWAS/reference matching, MAF filtering, LD pruning, and ambiguous-SNP removal.
3. SNP-by-cell and GWAS-overlapping peak-by-cell accessibility matrices.
4. Local SNP and peak LD covariance, Wilcoxon rank sum test scores, and balanced repeats.
5. Four LD-adjusted association models and an aggregated Cauchy combination.
6. Deterministic publication and QC figures in PNG/PDF with a figure manifest.

The plot stage can import a merged S-LDSC summary and produce the manuscript comparison panels.

## Important status

The numerical core and the complete workflow are tested with generated ATAC, GWAS, genotype, fragment, and S-LDSC inputs, including a real hg19-to-hg38 UCSC chain. 

`EPIC2` intentionally preserves its fitting convention: the response and complete raw design are transformed with the inverse-square-root covariance, the transformed intercept column is discarded, and the model is fit with an ordinary untransformed intercept plus the transformed focal and baseline predictors.

`EPIC2` calculates the reference MAF as `min(mean(dosage) / 2, 1 - mean(dosage) / 2)`. 

## Installation

```bash
Rscript scripts/install_dependencies.R
R CMD INSTALL .
```

The workflow requires R 4.3 or newer. Bioconductor dependencies are installed by the supplied script; `presto` is installed from its upstream GitHub repository. `rtracklayer` performs hg19/hg38 liftOver. A Conda environment specification is also included.

## Inputs

The ATAC input must be a Seurat/Signac object with:

- an `ATAC` count assay whose rows are `chr-start-end` peaks;
- attached fragment objects when SNP-level coverage or QC metrics are computed;
- `tissue` and `tissue_cell_type` metadata;
- `TSS.enrichment`, `FRIP`, `nCount_ATAC`, `BlacklistRatio`, and `nucleosome_signal` metadata if `compute_qc_metrics: false`.

Each GWAS table must already be allele-aligned to the reference panel and must contain rsID, chromosome, position, effect allele, other allele, MAF, effect estimate, standard error, and P-value columns. Genotypes are normalized to a variant-by-sample numeric dosage matrix. Supported in-memory encodings include:

- numeric hard calls or continuous dosage;
- allele-index GT calls such as `0/0`, `0/1`, and phased `0|1`;
- biallelic labels such as `AA`, `AB`, `BB`, `A/G`, and `A|G`;
- GP, GL, or PL values supplied as comma-delimited cells, a three-dimensional variant-by-sample-by-state array, or a list of state matrices;
- variants in rows or columns, with orientation inferred from GWAS rsIDs.

GWAS tables and genotype matrices can be supplied as `.rds`, `.rda/.RData`, `.csv[.gz]`, `.tsv[.gz]`, or `.txt[.gz]`. For a text genotype matrix, set `genotype_variant_id_column` to the column containing rsIDs or other stable IDs; the remaining columns are samples. Coordinates are taken from the matching GWAS rows, while dosage values are never changed by liftOver.

Set `genotype_format` explicitly whenever `auto` cannot distinguish the encoding. Carrier-only presence/absence data cannot yield exact MAF and are rejected. PLINK, VCF, and GDS containers must first be read into one of the supported R objects; the package does not guess file-level schemas.

Numeric `0/1` data with diploid ploidy are also intentionally ambiguous under `auto`: use `genotype_format: dosage` only when the values are allele counts, or use `genotype_ploidy: 1` for haploid calls.

## Genome builds

Both hg19/GRCh37 and hg38/GRCh38 are supported. Set `atac.genome` to the build of the ATAC peak and fragment coordinates, set `gwas.genome` or a per-GWAS `genome` to each summary-statistics build, and choose the common output build with `genome.target`.

When builds differ, provide the required directional UCSC chain files under `genome.chains`. ATAC-to-target maps retained peaks for downstream LD and models. For a GWAS that began in the ATAC build, its preserved input coordinates are reused for fragment queries; otherwise target-to-ATAC is also required. The original ATAC object and tabix fragment files remain in their input build; all saved GWAS, SNP IDs, peak IDs, LD objects, and model inputs use `genome.target`.

Only unique one-to-one mappings are retained by default. Every stage records unmapped, multimapped, duplicate-target, chain-path, and chain-MD5 information. The conversion changes coordinates only; users must still align GWAS alleles and genotype dosage coding before input.

## Configuration

Start from the installed example:

```r
file.copy(
  system.file("extdata/config.example.yml", package = "EPIC2"),
  "config.yml"
)
```

All relative paths are resolved relative to `config.yml`. Validate it before a
large run:

```bash
Rscript exec/EPIC2 validate config.yml --check-files
```

Cell QC, peak QC, reference missingness, MAF, and every SNPRelate pruning threshold used are configurable. Any GWAS-level setting can be overridden inside `gwas.traits.<name>`, and the effective values are saved with the per-GWAS QC counts. Figure formats, resolution, UMAP sampling, heatmap size, top-feature count, and S-LDSC input are also configurable. See `CONFIGURATION.md` for the complete key list.

`ld_prune_tool` is the backend selector, not a token. This release integrates `snprelate`; its `method`, missing-rate, MAF, window, start-position, and thread arguments are explicit YAML settings. For LD pruning done beforehand with PLINK or another program, provide matching pre-pruned inputs and set
`ld_prune: false`.

Harmonized GWAS files retain `EPIC2_reference_maf` and `EPIC2_reference_missing_rate`, input coordinates/build, and target build for variant-level QC auditing. Genome conversion summaries are written beside the GWAS and accessibility outputs.

## Running

Run all stages through publication figures:

```bash
Rscript exec/EPIC2 run config.yml
```

Run one stage or selected GWAS datasets:

```bash
Rscript exec/EPIC2 stage config.yml gwas --gwas=HDL
Rscript exec/EPIC2 stage config.yml ld --gwas=HDL,LDL
Rscript exec/EPIC2 run config.yml --stages=wrs,model --gwas=HDL
Rscript exec/EPIC2 stage config.yml plot --gwas=HDL,LDL
```

The equivalent R call for one configured GWAS is:

```r
config <- read_EPIC2_config("config.yml")
run_gwas_stage(config, gwas_name = "HDL")
```

Run the generated hg19-ATAC/GWAS to hg38 seven-stage integration test after installing all suggested packages and downloading the directional UCSC chain:

```bash
EPIC2_RUN_MIXED_BUILD_TEST=true \
EPIC2_HG19_TO_HG38_CHAIN=/path/to/hg19ToHg38.over.chain.gz \
Rscript tests/test-mixed-build-integration.R
```

The test creates its own tabix fragment file, Signac object, GWAS table, user-style TSV genotype matrix, balanced WRS runs, and S-LDSC summary. It then verifies target-build GWAS/peak outputs, the final model table, and every publication-figure family.

Every name under `gwas.traits` is processed independently and writes its own GWAS, genotype, and QC files. Running names separately updates rather than erases rows already present in `02_gwas/gwas_qc_summary.tsv`.

For Slurm, submit `ld`, `wrs`, and `model` as per-trait jobs after the shared ATAC and accessibility stages. Every stage writes only beneath the configured output directory.

## Output contract

```text
results/
├── 01_atac/             # filtered object and QC summaries
├── 02_gwas/             # per-GWAS harmonized data, genotype, and QC files
├── 02_accessibility/    # section 2 SNP/peak matrices and overlap tables
├── 03_ld/               # section 3 local covariance and inverse square roots
├── 04_wrs/              # section 4 full/balanced Wilcoxon and stability
├── 05_models/           # section 5 component models and Cauchy tables
└── figures/
    ├── 00_atac_qc/      # ATAC retention, metrics, and annotation diagnostics
    ├── 01_gwas_qc/      # per-GWAS filtering and liftOver diagnostics
    ├── 02_accessibility/ # SNP/peak count and coverage diagnostics
    ├── 03_ld/            # block sizes and bounded LD heatmaps
    ├── 04_wrs/           # WRS distributions, dotplots, and stability
    ├── 05_models/        # component and Cauchy model panels
    ├── 06_sldsc/         # optional S-LDSC and method-comparison panels
    ├── tables/           # data used to draw publication figures
    └── figure_manifest.tsv
```

The primary submission outputs are `05_models/tissue_cell_type.component_results.tsv` and
`05_models/tissue_cell_type.cauchy_results.tsv`. The six manuscript-facing figures are written at the root of `figures/` with the filenames referenced by the LaTeX manuscript. Each is saved as a 300-dpi PNG and a vector PDF by default. The full manuscript wrapper is in `../publication/`.

## Performance notes

Genotype imputation and standardization are vectorized, and local LD correlations are calculated in configurable BLAS blocks (`ld.correlation_block_size`) while retaining the same distance mask and
covariance definition. Group accessibility stays sparse and uses the column-wise aggregation that benchmarks faster for the study's roughly two dozen groups. Plotting uses deterministic UMAP downsampling and bounded LD heatmaps so diagnostics do not materialize unbounded dense objects.

Run `Rscript scripts/benchmark_hotspots.R` after installation to verify numerical equivalence and measure machine-specific speedups for genotype scaling and local LD construction.

## Reproducibility checklist

- Keep the YAML used for the manuscript under version control without private paths.
- Record the Git commit, package version, R `sessionInfo()`, and Slurm resources.
- Run `R CMD check EPIC2_0.1.0.tar.gz` on a clean environment.
- Add a small synthetic example that completes every stage in minutes.
- Tag the verified manuscript release and archive that tag and test data on Zenodo.

## License

MIT. Confirm dependency licenses and institutional approval before the public
release.
