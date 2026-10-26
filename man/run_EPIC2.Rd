% Maintained against the exported function signature in R/pipeline.R.
\name{run_EPIC2}
\alias{run_EPIC2}
\title{Run the EPIC2 analysis and publication-figure workflow}
\description{
Run the EPIC2 analysis and publication-figure workflow.
}
\usage{
run_EPIC2(config, stages = available_stages(), traits = NULL, 
    gwas_name = NULL) 
}
\arguments{
  \item{config}{YAML path or parsed configuration list.}
  \item{stages}{Character vector selecting workflow stages. Execution follows the canonical pipeline order.}
  \item{traits}{Optional GWAS-name subset retained for backward compatibility.}
  \item{gwas_name}{Optional name or names under \code{gwas.traits} to process.}
}
\details{
EPIC2 requires one explicitly declared common human genome assembly for GWAS, reference genotypes, ATAC peaks, and fragments. The main workflow validates declarations and coordinate bounds; it does not perform liftOver. Coordinate bounds alone cannot establish that inputs use the same assembly. Users must harmonize assemblies and allele orientation before analysis.

Stages execute in the canonical order \code{atac}, \code{gwas}, \code{accessibility}, \code{ld}, \code{wrs}, \code{model}, and \code{plot}, restricted to the requested subset. Selecting a later stage does not automatically run its prerequisites: compatible files from preceding stages must already exist. Outputs are written beneath the configured output directory, including a run manifest. Use either \code{traits} or \code{gwas_name}, not both.
}
\value{
A named list of stage return values, invisibly.
}
\seealso{\code{\link{read_EPIC2_config}}, \code{\link{validate_EPIC2_config}}, \code{\link{EPIC2_cli}}.}
