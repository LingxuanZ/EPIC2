% Maintained against the exported function signature in R/wrs.R.
\name{run_wrs_stage}
\alias{run_wrs_stage}
\title{Run section 4: SNP- and peak-level Wilcoxon analyses}
\description{
Run section 4: SNP- and peak-level Wilcoxon analyses.
}
\usage{
run_wrs_stage(config, traits = names(config$gwas$traits)) 
}
\arguments{
  \item{config}{Parsed EPIC2 configuration.}
  \item{traits}{Optional trait subset.}
}
\details{
Requires accessibility-stage feature-by-cell matrices and aligned cell annotations. The ordinary analysis compares each group with the remaining cells using \pkg{presto}; optional balanced repeats and stability summaries are controlled by configuration. Writes per-feature, per-group result tables in RDS and TSV formats, including one-sided P-values and raw and truncated one-sided Z-score columns. Optional balanced runs additionally write repeated results and summary and stability tables.
}
\value{
Named list of generated files.
}
\seealso{\code{\link{read_EPIC2_config}}, \code{\link{validate_EPIC2_config}}, \code{\link{run_EPIC2}}, \code{\link{EPIC2_cli}}.}
