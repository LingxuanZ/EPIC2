% Maintained against the exported function signature in R/gwas.R.
\name{harmonize_gwas}
\alias{harmonize_gwas}
\title{Harmonize and prune one GWAS with a reference genotype panel}
\description{
Harmonize and prune one GWAS with a reference genotype panel.
}
\usage{
harmonize_gwas(gwas, genotype, columns = EPIC2_defaults()$gwas$columns, 
    maf_min = 0.05, remove_ambiguous = TRUE, ld_prune = TRUE, 
    ld_threshold = 0.8, ld_window_bp = 500000L, seed = 1000L, 
    genotype_ploidy = 2, genotype_format = "auto", genotype_orientation = "auto", 
    genotype_missing_values = -9, genotype_variant_id_column = NULL, 
    ld_prune_noninteger = "error", gwas_name = NULL, source_genome = NULL, 
    target_genome = NULL, chain_file = NULL, genome_column = NULL, 
    liftover_unmapped_action = "drop", liftover_multimapping_action = "drop", 
    liftover_duplicate_action = "drop", reference_missing_rate_max = 0.01, 
    ld_prune_tool = "snprelate", ld_method = "composite", ld_slide_max_n = NULL, 
    ld_prune_maf_min = 0.005, ld_remove_monomorphic = TRUE, ld_start_pos = "random.f500", 
    ld_num_threads = 1L, genome_build = NULL, duplicate_action = liftover_duplicate_action) 
}
\arguments{
  \item{gwas}{GWAS summary-statistics data frame.}
  \item{genotype}{Variant-by-sample dosage matrix with rsIDs as row names.}
  \item{columns}{Column-name mapping from the configuration.}
  \item{maf_min}{Minimum GWAS and reference-panel minor allele frequency.}
  \item{remove_ambiguous}{Remove A/T and C/G variants.}
  \item{ld_prune}{Perform SNPRelate LD pruning.}
  \item{ld_threshold}{LD pruning threshold.}
  \item{ld_window_bp}{LD pruning window in base pairs.}
  \item{seed}{Random seed passed to SNPRelate.}
  \item{genotype_ploidy}{Maximum genotype dosage.}
  \item{genotype_format}{Genotype encoding accepted by \code{\link{normalize_genotype}}.}
  \item{genotype_orientation}{Location of variants in the input matrix.}
  \item{genotype_missing_values}{Input values representing missing genotypes.}
  \item{genotype_variant_id_column}{Optional variant ID column for data-frame genotypes.}
  \item{ld_prune_noninteger}{How SNPRelate pruning handles non-integer dosage.}
  \item{gwas_name}{Optional GWAS name used in diagnostics and the QC summary.}
  \item{source_genome}{Legacy build declaration; if supplied, must match \code{genome_build}.}
  \item{target_genome}{Legacy build declaration; if supplied, must match \code{genome_build}.}
  \item{chain_file}{Legacy conversion argument. Must be \code{NULL}; the pipeline never lifts coordinates.}
  \item{genome_column}{Optional GWAS column containing build metadata to verify.}
  \item{liftover_unmapped_action}{Legacy argument, unused in validation-only analysis.}
  \item{liftover_multimapping_action}{Legacy argument, unused in validation-only analysis.}
  \item{liftover_duplicate_action}{Legacy default for \code{duplicate_action}.}
  \item{reference_missing_rate_max}{Maximum missing-genotype rate per reference variant.}
  \item{ld_prune_tool}{LD-pruning implementation; currently \code{"snprelate"}.}
  \item{ld_method}{SNPRelate LD statistic: \code{"composite"}, \code{"r"}, \code{"dprime"}, or \code{"corr"}.}
  \item{ld_slide_max_n}{Optional maximum number of variants in a pruning window.}
  \item{ld_prune_maf_min}{SNPRelate's internal MAF threshold during pruning.}
  \item{ld_remove_monomorphic}{Ask SNPRelate to remove monomorphic variants.}
  \item{ld_start_pos}{SNPRelate window-start strategy.}
  \item{ld_num_threads}{Number of SNPRelate pruning threads.}
  \item{genome_build}{Required common input assembly, hg19/GRCh37 or hg38/GRCh38.}
  \item{duplicate_action}{Drop or reject every row at a duplicated coordinate.}
}
\details{
EPIC2 requires one explicitly declared common human genome assembly for GWAS, reference genotypes, ATAC peaks, and fragments. The main workflow validates declarations and coordinate bounds; it does not perform liftOver. Coordinate bounds alone cannot establish that inputs use the same assembly. Users must harmonize assemblies and allele orientation before analysis.

Reference variants are matched and reordered by variant identifier. GWAS identifiers missing from the reference panel are excluded. MAF, reference missingness, ambiguous-allele, coordinate, and optional LD-pruning filters are then applied. No allele-strand harmonization is inferred from the reference matrix. The actual retained SNP coordinates define single-base coordinate identifiers. For duplicated coordinates, the default \code{"drop"} excludes every row at that coordinate.
}
\value{
A list containing aligned \code{gwas}, \code{genotype}, and a QC summary.
}
