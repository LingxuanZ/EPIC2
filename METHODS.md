# Methods

EPIC2 tests trait associations with annotated cell groups using SNP- and peak-level GWAS responses and accessibility predictors. This document describes the default calculations; [CONFIGURATION.md](CONFIGURATION.md) lists adjustable parameters.

## Quality control

Cell filtering first requires TSS enrichment ≥ 2.5 and FRIP ≥ 0.3. The count cutoff is the first percentile of assay counts among those remaining cells. Subsequent filters require counts at or above that cutoff, blacklist ratio < 0.02, nucleosome signal ≤ 4, and at least 50 retained cells per group.

FRIP is the fraction of fragment records overlapping at least one peak. The blacklist ratio uses the same denominator and counts records overlapping at least one blacklist interval. Each record contributes at most once to either numerator, including when it spans multiple disjoint intervals. The fifth-column read multiplicity is ignored. Cells with no observed fragments receive `NA`.

Fragment files must be deduplicated, with raw barcodes correctly mapped to Seurat cell names. Distinct files must represent disjoint fragment collections; repeated exports of the same library require upstream deduplication. Repeated attachments of the same physical file are counted once per cell mapping. When QC calculation is disabled, supplied QC metadata is used as provided.

Peak filters retain nonzero autosomal peaks 100–1,000 bp wide with no blacklist overlap. Library-size adjustment uses counts across the retained autosomal peak background.

## GWAS and reference genotypes

Variants are matched by identifier and reordered together. Reference MAF is `min(f, 1 - f)`, where `f = mean(dosage, na.rm = TRUE) / ploidy`. Both GWAS and reference MAF must be at least 0.05; reference missingness must be at most 0.01. Coordinates, allele labels, and effect statistics are checked before LD pruning. Effect estimates and derived Z-scores must be finite, with finite positive standard errors.

Default pruning uses SNPRelate's composite statistic, a 500 kb window, threshold 0.8, and seed 1000. Strand-ambiguous variants are removed after pruning. Disabling pruning retains the other QC steps.

Numeric dosage must count the GWAS effect allele (`A1` by default). GT strings count allele `1`; GP/GL/PL states refer to the encoded counted allele and must have the same orientation. For allele-label inputs, harmonization verifies the observed labels and orients dosage to A1. Strand alignment must be established upstream. Signed LD orientation matters for peak-level quadratic forms even though squared SNP correlations are sign invariant.

Missing reference genotypes that pass filtering are mean-imputed before LD calculation. GL/PL conversion normalizes likelihood weights without applying population-prior probabilities. See the [genotype configuration](CONFIGURATION.md#gwas-and-reference-genotypes) for accepted encodings and missing-value rules.

## Accessibility predictors

SNP accessibility counts fragments covering each one-based SNP position; retained SNP rows have nonzero coverage and overlap retained peaks. Peak accessibility uses the existing count matrix restricted to GWAS-overlapping peaks, independently of SNP-position coverage.

Count predictors average library-size-adjusted feature counts over cells within each group. For SNP counts, let $x_{ij}$ denote the adjusted count in cell $i$ at SNP $j$. With $n_k$ retained cells in group $k$ and $N=\sum_{k=1}^{K}n_k$, the model uses

$$
A_{kj}=\frac{1}{n_k}\sum_{i\in\mathcal{C}_k}x_{ij},
\qquad
A_j=\frac{1}{N}\sum_{i=1}^{N}x_{ij}
   =\sum_{k=1}^{K}\frac{n_k}{N}A_{kj}.
$$

The same definitions apply to peak counts $x_{ip}$, giving $A_{kp}$ and $A_p$. The baseline is a pooled-cell mean. Cells with missing group labels or nonfinite/nonpositive library-size factors are excluded. Cells with zero counts at a feature remain in the denominators.

WRS compares each group with the remaining cells. A two-sided P-value is converted to an enrichment-oriented one-sided value: `P / 2` when AUC ≥ 0.5, and `1 - P / 2` otherwise. Its score is `qnorm(P_one_sided, lower.tail = FALSE)`. The 2nd and 98th percentiles are calculated from finite raw scores across feature–group pairs; all scores, including infinities, are then capped at those cutoffs. Both count and WRS models use the pooled-cell count baseline. WRS scores are not divided by group cell counts.

Optional balanced WRS repeats use equal numbers of target and background cells at each configured sample size. These runs provide stability diagnostics; the association models use the ordinary WRS scores.

## LD-adjusted association models

At SNP resolution, the response is the squared GWAS Z-score. The working covariance contains squared reference-genotype correlations with unit diagonal, restricted to same-chromosome pairs within 5 Mb.

For a peak containing `n_p` SNPs, let `W_p` be the stabilized inverse square root of its SNP correlation matrix. The response is `||W_p Z_p||² / n_p`. For peaks `p` and `q`, define `M = W_p R_pq W_q`, where `R_pq` is their cross-correlation matrix. The working covariance is `2 tr(M M') / (n_p n_q)` for same-chromosome peak pairs whose interval distance is at most 5 Mb.

Inverse square roots use an eigenvalue floor of `1e-8`, followed by ridge `0.01`. Response covariance matrices are decomposed by connected components; each component uses a dense eigendecomposition.

Each model includes an intercept, the focal group predictor, and a baseline predictor. The response and full design are transformed by the same inverse square root `W`:

```text
response:  W y
design:    [W 1, W focal, W baseline]
```

Coefficients, residual variance, and QR-based standard errors use this design. The focal coefficient is tested against a positive alternative using the residual-degrees-of-freedom t distribution. These are nominal P-values under the working regression model. A model with missing or nonfinite covariance-feature inputs returns an unavailable test with a diagnostic status; a reduced feature set requires a newly constructed covariance and inverse square root.

The four component models combine SNP/peak resolution with count/WRS predictors. ACAT combines available valid component P-values with equal weights, or requires all four when configured. BH and Bonferroni adjustments are applied across groups within each trait and grouping. The output records how many components contributed.

S-LDSC is an optional external comparison. Imported coefficient Z-scores are converted to positive-coefficient P-values using `pnorm(z, lower.tail = FALSE)`; EPIC2 does not run S-LDSC.

## Saved results and reproducibility

Use separate output directories when changing inputs, statistical settings, or package versions. Selecting stages reruns them; rerun dependent stages after changing their inputs.

GWAS and accessibility stages write completion records containing settings and file checksums. These records are verified when present. Manually prepared inputs without records require the user's own provenance checks. LD and WRS completion records are required before model fitting. Incomplete runs, changed artifacts, incompatible settings, or a changed package version prevent reuse.

Checksums stream saved artifacts and add disk I/O. Original fragment contents are not rehashed by these checks; preserve raw files and archive their checksums separately. Model runs save per-trait components and aggregate completed traits still listed in the configuration.

A stage lock prevents concurrent jobs in one output directory. After a crash, inspect `.EPIC2-stage-lock/owner.txt` and confirm that the process has stopped before removing a stale lock. Use separate output directories for independent concurrent jobs.

`config.used.yml` and `config.used.rds` store effective settings, including in-memory overrides. `run_history/` retains per-invocation settings and session manifests. Preserve input provenance and checksums, the source commit or archive checksum, and dependency versions with each analysis. The conda environment and dependency installer are installation recipes, not complete dependency locks.
