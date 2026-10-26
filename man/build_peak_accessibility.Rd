% Maintained against the exported function signature in R/accessibility.R.
\name{build_peak_accessibility}
\alias{build_peak_accessibility}
\title{Build a trait-specific peak-by-cell accessibility matrix}
\description{
Build a trait-specific peak-by-cell accessibility matrix.
}
\usage{
build_peak_accessibility(atac, gwas, columns = EPIC2_defaults()$gwas$columns, 
    assay = "ATAC", gwas_genome = NULL, atac_genome = NULL, gwas_to_atac_chain = NULL, 
    atac_to_gwas_chain = NULL, unmapped_action = "drop", multimapping_action = "drop", 
    duplicate_action = "drop", genome_build = NULL) 
}
\arguments{
  \item{atac}{Filtered Seurat/Signac object.}
  \item{gwas}{Harmonized GWAS table.}
  \item{columns}{GWAS column mapping.}
  \item{assay}{ATAC assay name.}
  \item{gwas_genome}{Legacy declaration; if supplied, must match \code{genome_build}.}
  \item{atac_genome}{Legacy declaration; if supplied, must match \code{genome_build}.}
  \item{gwas_to_atac_chain}{Legacy conversion argument; must be \code{NULL}.}
  \item{atac_to_gwas_chain}{Legacy conversion argument; must be \code{NULL}.}
  \item{unmapped_action}{Legacy argument, unused in validation-only analysis.}
  \item{multimapping_action}{Legacy argument, unused in validation-only analysis.}
  \item{duplicate_action}{Drop or reject duplicate target ranges.}
  \item{genome_build}{Required common assembly for GWAS, peaks, and fragments.}
}
\details{
EPIC2 requires one explicitly declared common human genome assembly for GWAS, reference genotypes, ATAC peaks, and fragments. The main workflow validates declarations and coordinate bounds; it does not perform liftOver. Coordinate bounds alone cannot establish that inputs use the same assembly. Users must harmonize assemblies and allele orientation before analysis.

This function selects existing post-QC ATAC peak-count rows overlapped by eligible GWAS SNPs. It does not recalculate peak counts from fragments and does not impose an additional 20-cells-per-peak or 20-peaks-per-cell filter.
}
\value{
A list with a sparse peak matrix and SNP-to-peak overlap table.
}
