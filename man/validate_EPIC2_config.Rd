% Maintained against the exported function signature in R/config.R.
\name{validate_EPIC2_config}
\alias{validate_EPIC2_config}
\title{Validate an EPIC2 configuration}
\description{
Validate an EPIC2 configuration.
}
\usage{
validate_EPIC2_config(config, check_files = FALSE) 
}
\arguments{
  \item{config}{Configuration list returned by \code{\link{read_EPIC2_config}}.}
  \item{check_files}{Whether configured input files must already exist.}
}
\details{
Validation returns a normalized copy. Retain this value, for example \code{config <- validate_EPIC2_config(config)}; the caller's list is not modified in place. The required \code{genome.build} field cannot be replaced by legacy build fields. Nonempty conversion-chain settings and conflicting assembly declarations cause errors.
}
\value{
A normalized configuration copy, invisibly, with one canonical \code{genome.build} and matching compatibility aliases. Callers must retain the returned value; the supplied configuration is not modified in place.
}
\seealso{\code{\link{read_EPIC2_config}}, \code{\link{run_EPIC2}}, \code{\link{EPIC2_cli}}.}
