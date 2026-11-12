# Methods

EPIC2 tests trait associations with annotated cell groups using SNP- and peak-level GWAS responses and accessibility predictors. This document summarizes the calculations implemented in EPIC2 using the package defaults. Full statistical details are provided in the Supplementary Methods of the accompanying manuscript. [CONFIGURATION.md](CONFIGURATION.md) lists adjustable parameters, and the [manuscript reproducibility guide](manuscript_reproducibility/README.md) describes the study-specific inputs and settings.

The full workflow runs `atac`, `gwas`, `accessibility`, `ld`, `wrs`, `model`, and `plot` in that order. GWAS variants, ATAC peaks, and fragment coordinates must already share the assembly declared in `genome.build`; the pipeline validates coordinates but does not perform liftOver.

## Quality control

Cell filtering first requires TSS enrichment ≥ 2.5 and FRIP ≥ 0.3. The count cutoff is the first percentile of assay counts among those remaining cells. Subsequent filters require counts at or above that cutoff, blacklist ratio < 0.02, nucleosome signal ≤ 4, and at least 50 retained cells per group.

FRIP is the fraction of fragment records overlapping at least one peak. The blacklist ratio uses the same denominator and counts records overlapping at least one blacklist interval. Each record contributes at most once to either numerator, including when it spans multiple disjoint intervals. The fifth-column read multiplicity is ignored. Cells with no observed fragments receive `NA`.

Fragment files must be deduplicated, with raw barcodes correctly mapped to Seurat cell names. Distinct files must represent disjoint fragment collections; repeated exports of the same library require upstream deduplication. Repeated attachments of the same physical file are counted once per cell mapping. When QC calculation is disabled, supplied QC metadata is used as provided.

Peak filters retain nonzero autosomal peaks 100–1,000 bp wide with no blacklist overlap.

## GWAS and reference genotypes

Variants are matched by identifier and reordered together. Reference MAF is `min(f, 1 - f)`, where `f = mean(dosage, na.rm = TRUE) / ploidy`. Both GWAS and reference MAF must be at least 0.05; reference missingness must be at most 0.01. Coordinates, allele labels, and effect statistics are checked before LD pruning. Effect estimates and derived Z-scores must be finite, with finite positive standard errors.

Default pruning uses SNPRelate's composite statistic, a 500 kb window, threshold 0.8, and seed 1000. Strand-ambiguous variants are removed after pruning. Disabling pruning retains the other QC steps.

Numeric dosage must count the GWAS effect allele (`A1` by default). GT strings count allele `1`; GP/GL/PL states refer to the encoded counted allele and must have the same orientation. For allele-label inputs, harmonization verifies the observed labels and orients dosage to A1. Strand alignment must be established upstream. Signed LD orientation matters for peak-level quadratic forms even though squared SNP correlations are sign invariant.

Missing reference genotypes that pass filtering are mean-imputed before LD calculation. GL/PL conversion normalizes likelihood weights without applying population-prior probabilities. See the [genotype configuration](CONFIGURATION.md#gwas-and-reference-genotypes) for accepted encodings and missing-value rules.

## Accessibility predictors

SNP accessibility counts fragments covering each one-based SNP position; retained SNP rows have nonzero coverage and overlap retained peaks. Peak accessibility uses the existing count matrix restricted to GWAS-overlapping peaks, independently of SNP-position coverage.

Count-based accessibility is the mean of library-size-adjusted feature counts across cells within each group. Library size is calculated across all retained autosomal peaks, before selecting GWAS-overlapping features. Counts in each cell are scaled by the median positive library size divided by that cell's library size.

The baseline covariate is the mean library-size-adjusted accessibility across all retained cells and is used in both the count and WRS models. Cells with missing group labels or nonfinite/nonpositive library-size factors are excluded from these means. Cells with zero counts at a feature remain in the denominators.

One-versus-rest Wilcoxon rank-sum (WRS) tests use the unadjusted feature-by-cell count matrices to assess higher accessibility in the focal group than in the remaining cells. The enrichment-oriented one-sided P-values are transformed to standard normal Z-scores. These scores are winsorized at the 2nd and 98th percentiles of finite scores across feature–group pairs; values beyond these cutoffs, including infinities, are capped. WRS scores are not divided by group cell counts.

Optional balanced WRS repeats use equal numbers of target and background cells at each configured sample size. These runs provide stability diagnostics; the association models use the ordinary WRS scores.

## LD-adjusted association models

At SNP resolution, the response is the squared GWAS Z-score. The working covariance contains squared reference-genotype correlations with unit diagonal, restricted to same-chromosome pairs within 5 Mb.

At peak resolution, GWAS Z-scores for the retained SNPs within each peak are combined into a quadratic-form statistic using the regularized inverse of within-peak LD, then divided by the number of SNPs. The working covariance between peak statistics is calculated from cross-peak reference LD. Off-diagonal entries are zero for peaks on different chromosomes or whose genomic intervals are separated by more than 5 Mb.

Inverse square roots use an eigenvalue floor of `1e-8`, followed by ridge `0.01`. Nonzero off-diagonal entries define connected components in each response covariance matrix; a stabilized inverse square root is computed for each component using a dense eigendecomposition.

Before ordinary least-squares fitting, the response, intercept column, focal accessibility predictor, and count-based baseline covariate are transformed by the same ridge-stabilized inverse square root of the working covariance matrix. Standard errors are calculated from the transformed fit, and the focal coefficient is tested against a positive alternative using a one-sided t-test with the residual degrees of freedom. These are nominal P-values under the working regression model.

A model with missing or nonfinite inputs for covariance features returns an unavailable test with a diagnostic status. Fitting a reduced feature set requires a newly constructed covariance and inverse square root.

The four component models combine SNP/peak resolution with count/WRS predictors. The Cauchy combination test (ACAT) combines available valid nominal component P-values with equal weights, or requires all four when configured. BH and Bonferroni adjustments are applied separately to each component test and the combined test across groups within each trait and grouping. The output records how many components contributed.

S-LDSC is an optional external comparison. Imported coefficient Z-scores are converted to one-sided P-values for positive coefficients using the upper tail of the standard normal distribution; EPIC2 does not run S-LDSC.

## Saved results and reproducibility

Use separate output directories when changing inputs, statistical settings, or package versions, and for independent concurrent runs. Selecting stages reruns them; rerun dependent stages after changing their inputs.

GWAS and accessibility stages write completion records containing settings and file checksums. These records are verified when present. Manually prepared inputs without records require the user's own provenance checks. LD and WRS completion records are required before model fitting. Incomplete runs, changed artifacts, incompatible settings, or a changed package version prevent reuse.

Original fragment contents are not rehashed by these checks; preserve raw files and archive their checksums separately. Model runs save per-trait components and aggregate completed traits still listed in the configuration.

A stage lock prevents concurrent jobs in one output directory. After a crash, inspect `.EPIC2-stage-lock/owner.txt` and confirm that the process has stopped before removing a stale lock.

`run_EPIC2()` writes effective settings, including in-memory overrides, to `config.used.yml` and `config.used.rds`, and retains per-invocation settings and session manifests in `run_history/`. Direct calls to individual `run_*_stage()` functions do not write these run-level records. Preserve input provenance and checksums, the source commit or archive checksum, and dependency versions with each analysis. The conda environment and dependency installer are installation recipes, not complete dependency locks.
