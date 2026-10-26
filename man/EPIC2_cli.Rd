% Maintained against the exported function signature in R/pipeline.R.
\name{EPIC2_cli}
\alias{EPIC2_cli}
\title{Command-line entry point}
\description{
Command-line entry point.
}
\usage{
EPIC2_cli(args = commandArgs(trailingOnly = TRUE)) 
}
\arguments{
  \item{args}{Command-line arguments, excluding the R executable.}
}
\details{
Commands are \code{validate CONFIG.yml [--check-files]}, \code{run CONFIG.yml [--stages=...] [--gwas=...]}, and \code{stage CONFIG.yml STAGE [--gwas=...]}. Comma-separated names select multiple stages or GWAS traits. The historical \code{--traits} selector is accepted, but cannot be combined with \code{--gwas}. Calling with no arguments, \code{help}, \code{-h}, or \code{--help} prints usage.
}
\value{
Invisibly, the normalized configuration for \code{validate}, the workflow result for \code{run} or \code{stage}, or \code{NULL} after printing help. Invalid commands or configurations raise an R error.
}
\seealso{\code{\link{read_EPIC2_config}}, \code{\link{validate_EPIC2_config}}, \code{\link{run_EPIC2}}.}
\examples{
EPIC2_cli("--help")
}
