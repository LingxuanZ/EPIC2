% Maintained against the exported function signature in R/model.R.
\name{run_model_stage}
\alias{run_model_stage}
\title{LD-adjusted models and Cauchy combination}
\description{
LD-adjusted models and Cauchy combination.
}
\usage{
run_model_stage(config, traits = names(config$gwas$traits)) 
}
\arguments{
  \item{config}{Parsed EPIC2 configuration.}
  \item{traits}{Optional trait subset.}
}
\details{
Requires aligned accessibility, LD, and WRS outputs. Fits SNP-count, SNP-WRS, peak-count, and peak-WRS models, applies within-trait corrections, and combines the four component P-values with ACAT for each configured grouping.
}
\value{
Named list of output tables by grouping.
}
\seealso{\code{\link{read_EPIC2_config}}, \code{\link{validate_EPIC2_config}}, \code{\link{run_EPIC2}}, \code{\link{EPIC2_cli}}.}
