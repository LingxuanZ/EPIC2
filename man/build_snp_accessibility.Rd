% Maintained against the exported function signature in R/accessibility.R.
\name{build_snp_accessibility}
\alias{build_snp_accessibility}
\title{Build a trait-specific SNP-by-cell accessibility matrix}
\description{
Build a trait-specific SNP-by-cell accessibility matrix.
}
\usage{
build_snp_accessibility(atac, gwas, columns = EPIC2_defaults()$gwas$columns, 
    assay = "ATAC", workers = 1L, chunks = 20L, gwas_genome = NULL, 
    atac_genome = NULL, chain_file = NULL, unmapped_action = "drop", 
    multimapping_action = "drop", duplicate_action = "drop", 
    genome_build = NULL) 
}
\arguments{
  \item{atac}{Filtered Seurat/Signac object with fragments attached.}
  \item{gwas}{Harmonized GWAS table containing \code{coord_id}.}
  \item{columns}{GWAS column mapping.}
  \item{assay}{ATAC assay name.}
  \item{workers}{Number of BiocParallel workers.}
  \item{chunks}{Number of genomic feature chunks.}
  \item{gwas_genome}{Legacy declaration; if supplied, must match \code{genome_build}.}
  \item{atac_genome}{Legacy declaration; if supplied, must match \code{genome_build}.}
  \item{chain_file}{Legacy conversion argument; must be \code{NULL}.}
  \item{unmapped_action}{Legacy argument, unused in validation-only analysis.}
  \item{multimapping_action}{Legacy argument, unused in validation-only analysis.}
  \item{duplicate_action}{Drop or reject duplicate target ranges.}
  \item{genome_build}{Required common assembly for GWAS, peaks, and fragments.}
}
\details{
EPIC2 requires one explicitly declared common human genome assembly for GWAS, reference genotypes, ATAC peaks, and fragments. The main workflow validates declarations and coordinate bounds; it does not perform liftOver. Coordinate bounds alone cannot establish that inputs use the same assembly. Users must harmonize assemblies and allele orientation before analysis.

Coverage is computed from attached fragments at single-base SNP ranges lying within retained ATAC peaks. Zero-signal SNP rows are omitted. This function does not impose additional detection thresholds of 20 cells per feature or 10 SNPs per cell.
}
\value{
Sparse SNP-by-cell matrix containing nonzero SNPs in ATAC peaks.
}
