# Configuration reference

All paths may be absolute or relative to the YAML configuration file. Values
shown below are defaults unless marked as required.

## Project and execution

| Key | Default | Purpose |
| --- | --- | --- |
| `project.name` | `EPIC2 analysis` | Human-readable run name. |
| `project.output_dir` | `results` | Root for all generated files. |
| `parameters.seed` | `1000` | Reproducible LD pruning and balanced sampling seed. |
| `parameters.workers` | `1` | BiocParallel workers for fragment overlap calculations. |

## Common genome assembly

| Key | Default | Purpose |
| --- | --- | --- |
| `genome.build` | **required; no default** | Assembly already shared by GWAS, ATAC peaks, and fragment coordinates. Accepts hg19/GRCh37 or hg38/GRCh38. |
| `genome.duplicate_action` | `drop` | `drop` removes all rows at duplicated coordinates; `error` rejects them. |

EPIC2 does not convert coordinates in its analysis pipeline. Harmonize all
coordinate-bearing inputs before running it. Missing `genome.build` is an error,
even when old assembly fields are present. GRCh37/hg19 and GRCh38/hg38 aliases
are normalized to hg19 and hg38, respectively.

Legacy `genome.target`, `atac.genome`, `gwas.genome`, and per-trait `genome`
declarations are accepted only when each matches the common assembly. Missing
legacy values inherit it. A conflicting, empty, or invalid explicit declaration
is an error. Nonempty `genome.chains` is rejected; remove chain settings from
pipeline configurations. The standalone `lift_genomic_ranges()` utility remains
available for separate, explicitly requested preprocessing and is never called
by the pipeline.

Assembly metadata and chromosome bounds are checked where available. Missing
metadata is reported as unavailable; coordinate bounds alone cannot identify a
true assembly. Genotype dosage values do not require coordinate conversion, and
variant identifier matching does not establish allele coding alignment.

## ATAC input and QC

| Key | Default | Purpose |
| --- | --- | --- |
| `atac.input` | required | Seurat/Signac `.rds` or `.rda` object. |
| `atac.object_name` | null | Object name when an `.rda` contains multiple objects. |
| `atac.assay` | `ATAC` | Chromatin assay name. |
| `atac.genome` | inherits `genome.build` | Optional legacy declaration; must match the common assembly. |
| `atac.tissue_column` | `tissue` | Tissue metadata column. |
| `atac.group_column` | `tissue_cell_type` | Column used for minimum group-size QC. |
| `atac.attach_fragments` | `false` | Attach files listed under `atac.fragment_files`. |
| `atac.fragment_files` | null | Named tissue-to-fragment-file mapping. |
| `atac.validate_fragments` | `true` | Ask Signac to validate cell barcodes during attachment. |
| `atac.compute_qc_metrics` | `false` | Compute nucleosome, TSS, FRIP, and blacklist metrics. |
| `atac.recompute_embedding` | `false` | Re-run TF-IDF, LSI, UMAP, neighbors, and clustering after QC. |
| `atac.qc.tss_min` | `2.5` | Minimum TSS enrichment. |
| `atac.qc.frip_min` | `0.3` | Minimum fraction of fragments in peaks. |
| `atac.qc.count_quantile_min` | `0.01` | Remove cells below this ATAC-count quantile, calculated after TSS and FRIP filtering. |
| `atac.qc.blacklist_max` | `0.02` | Exclusive upper blacklist-fragment ratio. |
| `atac.qc.nucleosome_max` | `4` | Maximum nucleosome signal. |
| `atac.qc.min_cells_per_group` | `50` | Minimum retained cells per annotation group. |
| `atac.qc.peak_width_min` | `100` | Minimum peak width in base pairs. |
| `atac.qc.peak_width_max` | `1000` | Maximum peak width in base pairs. |
| `atac.qc.peak_min_total_counts` | `1` | Minimum total counts across cells for a peak. |
| `atac.qc.peak_min_cells` | `1` | Minimum number of cells in which a peak is detected. |
| `atac.qc.peak_autosomes_only` | `true` | Retain peaks on chromosomes 1-22 only. |
| `atac.qc.peak_remove_blacklist` | `true` | Remove peaks using the Signac blacklist matching `genome.build`. |

When `compute_qc_metrics` is false, the five QC metadata fields must already be
present. The count field is derived from `atac.assay`: `nCount_ATAC` for `ATAC`,
`nCount_peaks` for `peaks`, and likewise for other assay names. Fragment objects must be attached before the accessibility stage even
if QC metrics are precomputed. The ATAC stage writes both the configured values
and the data-derived ATAC count cutoff to `01_atac/qc_thresholds.tsv`.

## GWAS and reference genotypes

| Key | Default | Purpose |
| --- | --- | --- |
| `gwas.traits.<name>.summary` | required | Aligned GWAS table for one trait. |
| `gwas.traits.<name>.genotype` | required | Matrix-like genotype object or probability array/list. |
| `summary_object` / `genotype_object` | null | Object names for multi-object `.rda` inputs. |
| `gwas.genome` | inherits `genome.build` | Optional legacy declaration; every per-GWAS declaration must also match. |
| `gwas.genome_column` | null | Optional summary-statistics column containing build metadata to verify. |
| `gwas.maf_min` | `0.05` | Minimum MAF in both GWAS and reference panel. |
| `gwas.genotype_format` | `auto` | `dosage`, `hard_call`, `gt`, `allele`, `gp`, `gl`, `pl`, or automatic detection. |
| `gwas.genotype_orientation` | `auto` | `variants_by_samples`, `samples_by_variants`, or infer from rsIDs. |
| `gwas.genotype_ploidy` | `2` | Allele copy number used to derive reference frequency. |
| `gwas.genotype_missing_values` | `[-9]` | Additional missing-value codes; finite numeric GL states remain likelihood values (see below). |
| `gwas.genotype_variant_id_column` | null | Variant ID column for genotype data-frame input. |
| `gwas.reference_missing_rate_max` | `0.01` | Maximum missing reference genotypes per variant. |
| `gwas.remove_ambiguous` | `true` | Remove A/T, T/A, C/G, and G/C variants after pruning (documented order). |
| `gwas.ld_prune` | `true` | Perform SNPRelate pruning. |
| `gwas.ld_prune_tool` | `snprelate` | LD-pruning backend; only SNPRelate is integrated in this release. |
| `gwas.ld_prune_noninteger` | `error` | Reject continuous dosage during pruning, or `round` to opt into nearest hard calls. |
| `gwas.ld_method` | `composite` | SNPRelate statistic: `composite`, `r`, `dprime`, or `corr`. |
| `gwas.ld_threshold` | `0.8` | SNPRelate LD threshold. |
| `gwas.ld_window_bp` | `500000` | LD pruning window in base pairs. |
| `gwas.ld_slide_max_n` | null | Optional maximum variants in each pruning window. |
| `gwas.ld_prune_maf_min` | `0.005` | Additional SNPRelate MAF filter during pruning. |
| `gwas.ld_remove_monomorphic` | `true` | Ask SNPRelate to remove monomorphic variants. |
| `gwas.ld_start_pos` | `random.f500` | Window start: `random.f500`, `random`, `first`, or `last`. |
| `gwas.ld_num_threads` | `1` | Threads used by `snpgdsLDpruning()`. |

`gwas.columns` maps `rsid`, `chr`, `pos`, `a1`, `a2`, `maf`, `p`, `beta`,
and `se` to input columns. A GWAS may override the column mapping and every
GWAS QC, genotype, and LD-pruning setting under `gwas.traits.<name>`. The
effective settings and variant counts are written to each
`02_gwas/<name>.gwas_qc.tsv` and to the combined QC summary.

GWAS and genotype inputs may be `.rds`, `.rda/.RData`, `.csv[.gz]`,
`.tsv[.gz]`, or `.txt[.gz]`. A text genotype table should place stable variant
IDs in the configured `genotype_variant_id_column` and samples in the remaining
columns. Genotype matrices require unique, nonmissing variant IDs in row names
(or column names for transposed input) or the configured ID column. Anonymous
variants are rejected even when dimensions match the GWAS or orientation is
specified; GWAS IDs are orientation hints, not replacement reference IDs. If an
upstream format omits names, assign the actual reference variant IDs only after
verifying their order. Genotypes are joined by ID, genuinely absent variants are
skipped, and retained rows are reordered together with the GWAS. Coordinates come
from the already harmonized GWAS table, and dosage values are not converted.
Native PLINK/VCF/BGEN/GDS containers still need to be read or exported to a
supported matrix first.

`auto` accepts numeric allele dosage, biallelic GT strings, allele labels, and
well-formed GP probabilities. GL and PL must be declared explicitly because
their comma-delimited representation cannot be distinguished safely from GP.
Logical or declared `carrier` matrices are rejected: carrier status does not
identify heterozygous versus homozygous carriers and therefore cannot recover
exact allele frequency. Numeric `0/1` with diploid ploidy is likewise ambiguous
under `auto`; declare `dosage` only for true allele counts, or set ploidy to `1`
for haploid calls.

For explicit `gl` input, finite numeric states such as `-9` are valid log10
likelihoods and are never replaced by numeric missing codes inside a likelihood
vector, array, or state-matrix list. For example, `0,-9,-18` is a called genotype
and `-9,-9,-9` is an uninformative called vector. Use `NA` or nonnumeric missing
tokens for missing array/list states. A scalar text genotype cell consisting
only of a configured code, such as `"-9"`, is still missing. Other genotype
formats retain their existing missing-code behavior.

Finite effect estimates, finite positive standard errors, and finite derived
Z-scores are required before LD pruning, so an unusable summary-statistic row
cannot displace a usable linked variant.

SNPRelate pruning operates on diploid `0/1/2` hard calls. Haploid dosage is
scaled exactly to `0/2`. Continuous dosage, GP/GL/PL-derived expected dosage, and
polyploid values that do not map exactly require either
`ld_prune_noninteger: round` or `ld_prune: false`. The normalized dosage is
retained unchanged for downstream LD calculations even when pruning uses the
explicitly requested nearest hard calls.

The first argument of `snpgdsLDpruning()` is an open GDS file handle created
internally; it is not a token. The user-facing selector is `ld_prune_tool`, and
the LD statistic is selected with `ld_method`. To prune with PLINK or another
external program in this release, prepare matching pre-pruned GWAS and genotype
inputs first and set `ld_prune: false`; arbitrary external commands are not run
from YAML.

The package does not infer allele flips from a dosage matrix without allele
metadata. Summary statistics must be aligned to the reference panel before
input.

## Accessibility, LD, and tests

| Key | Default | Purpose |
| --- | --- | --- |
| `accessibility.chunks` | `20` | Chunks used for SNP FeatureMatrix calculation. |
| `ld.window_bp` | `5000000` | Same-chromosome local covariance distance. |
| `ld.ridge` | `0.01` | Ridge added before inverse square roots. |
| `ld.eigen_tolerance` | `1e-8` | Eigenvalue floor before ridge regularization. |
| `ld.correlation_block_size` | `256` | Focal variants per vectorized BLAS correlation block; performance-only setting. |
| `ld.max_snps_for_chr_ld` | `15000` | Maximum variants for optional chromosome LD caching; also capped by `ld.max_dense_block_size`. Above the smaller limit, calculate peak-pair LD on demand. |
| `ld.max_dense_block_size` | `5000` | Reject a larger dense inverse-square-root block before allocation; applies to SNP/peak covariance components and within-peak SNP inverses. Positive integer, or `.inf` to deliberately disable the guard. |
| `wrs.truncation_quantiles` | `[0.02, 0.98]` | Quantiles of finite raw enrichment-oriented Z-scores; then cap all scores, including infinities. |
| `wrs.balanced.enabled` | `false` | Run repeated equal-size one-versus-rest sampling. |
| `wrs.balanced.cells_per_side` | `[500, 1000]` | Target and background sample sizes. |
| `wrs.balanced.repeats` | `30` | Repeats at each sample size. |
| `wrs.balanced.auc_cutoff` | `0.55` | Selection threshold for repeat-pair stability. |
| `wrs.balanced.stability_frequency_cutoff` | `0.8` | Minimum selected-repeat proportion for a feature to be marked stable. |
| `wrs.balanced.top_n` | `100` | Top-ranked features compared between repeat pairs. |
| `model.groupings` | `[tissue_cell_type, tissue]` | Metadata resolutions modeled independently. |
| `model.require_all_components` | `false` | Combine available valid P-values by default (default); set true to require all four. |

Changing `ld.ridge`, `ld.window_bp`, Wilcoxon truncation, or group definitions
changes the statistical model and should be treated as a new analysis rather
than an execution-only option.

`ld.correlation_block_size` changes batching only. Block size 1 and larger
blocks use the same standardized dosage values, local-distance mask, and output
schema; numerical-equivalence tests cover this contract.

A local LD distance window does not bound the size of a connected component:
successive local connections can span a chromosome. Each component still uses
a dense eigendecomposition. `ld.max_dense_block_size` checks every component
before any component is converted to dense form, and checks each within-peak
SNP block before allocating its dense LD matrix. A stopped run reports the
dimension and memory estimates. At the default limit of 5000, one double matrix
occupies approximately 0.19 GiB; an illustrative six-matrix working estimate is
1.12 GiB. LAPACK workspaces, genotype data, retained blocks, covariance storage,
and copies can require additional memory, so this is not a total-memory limit
or a guarantee that a job will fit in RAM.

The guard can stop large analyses that earlier versions attempted. After
checking the resources available to the job, increase the limit explicitly or
set `ld.max_dense_block_size: .inf` to disable it (`Inf` in R). The standalone
`compute_inverse_sqrt()` and `build_peak_ld_model()` functions accept the same
`max_dense_block_size` argument. The guard does not split components, remove LD
connections, or approximate accepted inverse square roots. It is not a scalable
replacement for the dense eigensolver. Reducing `ld.window_bp` to bypass a large
block changes the statistical model.

Chromosome LD caching uses the smaller of `ld.max_snps_for_chr_ld` and the dense
block limit. Above it, the existing exact peak-pair calculation is used. This
can increase runtime while avoiding the optional full chromosome LD matrix;
it does not change the covariance definition. Sparse covariance construction
and the accumulated output matrices can still exhaust memory independently of
this guard.

## Publication figures

| Key | Default | Purpose |
| --- | --- | --- |
| `plots.formats` | `[png, pdf]` | One or more of `png` and `pdf`. |
| `plots.dpi` | `300` | Raster resolution for PNG output. |
| `plots.base_size` | `11` | Base ggplot text size. |
| `plots.max_umap_cells` | `100000` | Deterministic maximum nuclei per UMAP panel. |
| `plots.count_bin_max` | `4` | Last exact accessibility-count bin; larger values form one `N+` bin. |
| `plots.ld_heatmap_max_features` | `250` | Maximum SNPs/peaks shown in a dense diagnostic LD heatmap. |
| `plots.wrs_histogram_bins` | `50` | WRS histogram bins. |
| `plots.top_features_per_group` | `3` | WRS-ranked features selected per trait and annotation group. |
| `plots.pvalue_cap` | `12` | Display cap for `-log10(P)`; raw table values are unchanged. |
| `plots.sldsc_summary` | null | Optional merged S-LDSC TSV/RDS/RData used for method-comparison figures. |
| `plots.fail_on_error` | `false` | Stop immediately on a plotting error; use `true` for release runs. |

The `plot` stage reads saved results and does not modify analytical objects.
Every attempted figure is recorded in `figures/figure_manifest.tsv` with its
source path, output paths, status, message, and elapsed time. S-LDSC computation
is external; accepted summaries need `trait`, an annotation/group column, and a
coefficient z-score. Common coefficient/enrichment P-value column names are
standardized automatically.
