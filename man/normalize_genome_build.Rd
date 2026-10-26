% Maintained against the exported function signature in R/genome.R.
\name{normalize_genome_build}
\alias{normalize_genome_build}
\title{Normalize a human genome-build name}
\description{
Normalize a human genome-build name.
}
\usage{
normalize_genome_build(genome) 
}
\arguments{
  \item{genome}{Genome build name such as \code{hg19}, \code{GRCh37}, \code{hg38}, or \code{GRCh38}.}
}
\value{
Canonical UCSC build name, either \code{hg19} or \code{hg38}.
}
\seealso{\code{\link{validate_genome_coordinates}}, \code{\link{lift_genomic_ranges}}.}
\examples{
normalize_genome_build("GRCh37")
normalize_genome_build("hg38")
}
