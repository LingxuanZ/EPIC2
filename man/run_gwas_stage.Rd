% Maintained against the exported function signature in R/gwas.R.
\name{run_gwas_stage}
\alias{run_gwas_stage}
\title{GWAS harmonization and LD pruning}
\description{
GWAS harmonization and LD pruning.
}
\usage{
run_gwas_stage(config, traits = NULL, gwas_name = NULL) 
}
\arguments{
  \item{config}{Parsed EPIC2 configuration.}
  \item{traits}{Optional GWAS-name subset retained for backward compatibility.}
  \item{gwas_name}{Optional name or names under \code{gwas.traits} to process.}
}
\details{
EPIC2 requires one explicitly declared common human genome assembly for GWAS, reference genotypes, ATAC peaks, and fragments. The main workflow validates declarations and coordinate bounds; it does not perform liftOver. Coordinate bounds alone cannot establish that inputs use the same assembly. Users must harmonize assemblies and allele orientation before analysis.

Reads each requested GWAS/reference pair and writes the retained GWAS table, aligned reference dosage matrix, and filter summary. Prerequisites and file names follow the package example YAML and \code{CONFIGURATION.md}.
}
\value{
Named list of output paths.
}
\seealso{\code{\link{read_EPIC2_config}}, \code{\link{validate_EPIC2_config}}, \code{\link{run_EPIC2}}, \code{\link{EPIC2_cli}}.}
