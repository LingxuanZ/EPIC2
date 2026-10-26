% Maintained against the exported function signature in R/plots.R.
\name{run_plot_stage}
\alias{run_plot_stage}
\title{Generate publication figures from completed EPIC2 stages}
\description{
Plotting is isolated from statistical computation: this stage reads frozen
pipeline outputs and writes PNG/PDF figures plus a machine-readable manifest.
}
\usage{
run_plot_stage(config, traits = names(config$gwas$traits)) 
}
\arguments{
  \item{config}{Parsed EPIC2 configuration.}
  \item{traits}{Optional GWAS-name subset.}
}
\details{
Requires the stage outputs needed by the requested figures. Statistical results are read rather than refitted. Figure files in the configured PNG/PDF formats, supporting tables, and a manifest are written below the \code{figures} subdirectory of \code{config$output_dir}.
}
\value{
The figure manifest, invisibly.
}
\seealso{\code{\link{read_EPIC2_config}}, \code{\link{validate_EPIC2_config}}, \code{\link{run_EPIC2}}, \code{\link{EPIC2_cli}}.}
