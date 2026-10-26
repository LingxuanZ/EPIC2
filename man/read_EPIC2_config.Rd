% Maintained against the exported function signature in R/config.R.
\name{read_EPIC2_config}
\alias{read_EPIC2_config}
\title{Read an EPIC2 YAML configuration}
\description{
Read an EPIC2 YAML configuration.
}
\usage{
read_EPIC2_config(path, check_files = FALSE) 
}
\arguments{
  \item{path}{Path to a YAML configuration file.}
  \item{check_files}{Whether to verify configured input files immediately.}
}
\details{
The YAML field \code{genome.build} is required; there is no default assembly. Accepted aliases are hg19/GRCh37 and hg38/GRCh38. Optional legacy build fields must agree with this value. Nonempty chain-file settings are rejected. Relative input and output paths are resolved against the configuration file directory.
}
\value{
A validated configuration list with defaults filled in.
}
\seealso{\code{\link{validate_EPIC2_config}}, \code{\link{run_EPIC2}}, \code{\link{EPIC2_cli}}.}
