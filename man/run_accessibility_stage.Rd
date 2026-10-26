% Maintained against the exported function signature in R/accessibility.R.
\name{run_accessibility_stage}
\alias{run_accessibility_stage}
\title{Construct SNP and peak accessibility matrices}
\description{
Construct SNP and peak accessibility matrices.
}
\usage{
run_accessibility_stage(config, traits = names(config$gwas$traits)) 
}
\arguments{
  \item{config}{Parsed EPIC2 configuration.}
  \item{traits}{Optional trait subset.}
}
\details{
EPIC2 requires one explicitly declared common human genome assembly for GWAS, reference genotypes, ATAC peaks, and fragments. The main workflow validates declarations and coordinate bounds; it does not perform liftOver. Coordinate bounds alone cannot establish that inputs use the same assembly. Users must harmonize assemblies and allele orientation before analysis.

Requires the completed ATAC and GWAS stage outputs. Writes trait-specific SNP and peak matrices, a genotype matrix aligned to the SNP accessibility rows, SNP-to-peak overlap mappings, and genome-validation QC summaries. The upstream GWAS stage excludes variants truly absent from the reference panel. The retained identifiers are used consistently by downstream response and LD construction.
}
\value{
Named list of generated files.
}
\seealso{\code{\link{read_EPIC2_config}}, \code{\link{validate_EPIC2_config}}, \code{\link{run_EPIC2}}, \code{\link{EPIC2_cli}}.}
