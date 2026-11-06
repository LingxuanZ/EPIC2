\name{run_wrs_stage}
\alias{run_wrs_stage}
\title{Run SNP- and peak-level Wilcoxon analyses}
\description{
Run SNP- and peak-level Wilcoxon analyses.
}
\usage{
run_wrs_stage(config, traits = names(config$gwas$traits)) 
}
\arguments{
  \item{config}{Parsed EPIC2 configuration.}
  \item{traits}{Optional trait subset.}
}
\details{
Requires accessibility-stage feature-by-cell matrices and aligned cell annotations. The ordinary analysis compares each group with the remaining cells using \pkg{presto}; optional balanced repeats and stability summaries are controlled by configuration. Writes per-feature, per-group result tables as RDS, including one-sided P-values and raw and truncated one-sided Z-score columns. Matching TSV exports are controlled by \code{config$wrs$write_tsv}, which defaults to \code{TRUE}. Set it to \code{FALSE} to omit these exports and remove matching older TSV files after a successful trait run. Optional balanced runs additionally write repeated results and summary and stability tables; these are unaffected by \code{write_tsv}.
}
\value{
Named list of generated files.
}
\seealso{\code{\link{read_EPIC2_config}}, \code{\link{validate_EPIC2_config}}, \code{\link{run_EPIC2}}, \code{\link{EPIC2_cli}}.}
