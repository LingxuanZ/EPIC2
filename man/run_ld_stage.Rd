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

Set \code{config$ld$save_peak_checkpoints} to \code{FALSE} to omit the diagnostic
\code{.peak_response.rds} and \code{.peak_covariance.rds} files; the default is
\code{TRUE}. After a successful trait run, older copies of these two files are
removed when disabled. The complete peak model and public return values are
unchanged. Checkpoints aid diagnosis after a failure and are not automatically reused.
}
\value{
Named list of generated files.
}
\seealso{\code{\link{read_EPIC2_config}}, \code{\link{validate_EPIC2_config}}, \code{\link{run_EPIC2}}, \code{\link{EPIC2_cli}}.}
