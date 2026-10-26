% Maintained against the exported function signature in R/genome.R.
\name{validate_genome_coordinates}
\alias{validate_genome_coordinates}
\title{Validate human genomic coordinates against a declared build}
\description{
Coordinate bounds can detect impossible declarations, but cannot distinguish
hg19 from hg38 for positions valid in both builds. Input builds must therefore
be declared explicitly.
}
\usage{
validate_genome_coordinates(chr, start, end = start, genome, label = "coordinates") 
}
\arguments{
  \item{chr}{Chromosome values.}
  \item{start}{One-based inclusive starts.}
  \item{end}{One-based inclusive ends. Defaults to \code{start}.}
  \item{genome}{Declared input build.}
  \item{label}{Input label used in errors.}
}
\value{
A one-row coordinate QC data frame.
}
\seealso{\code{\link{normalize_genome_build}}, \code{\link{lift_genomic_ranges}}.}
\examples{
validate_genome_coordinates("chr1", 100L, genome = "hg19")
}
