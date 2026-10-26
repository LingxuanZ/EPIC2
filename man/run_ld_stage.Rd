% Maintained against the exported function signature in R/ld.R.
\name{run_ld_stage}
\alias{run_ld_stage}
\title{Build local SNP- and peak-level LD models}
\description{
Build local SNP- and peak-level LD models.
}
\usage{
run_ld_stage(config, traits = names(config$gwas$traits)) 
}
\arguments{
  \item{config}{Parsed EPIC2 configuration.}
  \item{traits}{Optional trait subset.}
}
\details{
EPIC2 requires one explicitly declared common human genome assembly for GWAS, reference genotypes, ATAC peaks, and fragments. The main workflow validates declarations and coordinate bounds; it does not perform liftOver. Coordinate bounds alone cannot establish that inputs use the same assembly. Users must harmonize assemblies and allele orientation before analysis.

Requires accessibility-stage matrices, overlap mappings, and aligned reference genotypes. Writes SNP-level and peak-level LD response covariance objects and their stabilized inverse square roots for the retained feature sets.
}
\value{
Named list of generated files.
}
\seealso{\code{\link{read_EPIC2_config}}, \code{\link{validate_EPIC2_config}}, \code{\link{run_EPIC2}}, \code{\link{EPIC2_cli}}.}
