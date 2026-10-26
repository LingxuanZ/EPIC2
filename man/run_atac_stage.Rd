% Maintained against the exported function signature in R/atac.R.
\name{run_atac_stage}
\alias{run_atac_stage}
\title{ATAC preparation and QC}
\description{
ATAC preparation and QC.
}
\usage{
run_atac_stage(config) 
}
\arguments{
  \item{config}{Parsed EPIC2 configuration.}
}
\details{
EPIC2 requires one explicitly declared common human genome assembly for GWAS, reference genotypes, ATAC peaks, and fragments. The main workflow validates declarations and coordinate bounds; it does not perform liftOver. Coordinate bounds alone cannot establish that inputs use the same assembly. Users must harmonize assemblies and allele orientation before analysis.

This stage reads the configured Seurat/Signac input, attaches configured fragments, applies the configured cell and peak QC filters, optionally recomputes reductions, and writes the QC object and diagnostic tables. For a workflow starting from a prepared post-QC object, omit this stage and provide the expected upstream output for later stages.
}
\value{
The path to the filtered ATAC RDS object.
}
\seealso{\code{\link{read_EPIC2_config}}, \code{\link{validate_EPIC2_config}}, \code{\link{run_EPIC2}}, \code{\link{EPIC2_cli}}.}
