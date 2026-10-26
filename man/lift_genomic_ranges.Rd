% Maintained against the exported function signature in R/genome.R.
\name{lift_genomic_ranges}
\alias{lift_genomic_ranges}
\title{Lift genomic ranges between hg19 and hg38}
\description{
Only ranges with exactly one mapping are retained by default. A directional
UCSC chain file must be supplied whenever source and target builds differ.
}
\usage{
lift_genomic_ranges(ranges, source_genome, target_genome, chain_file = NULL, 
    unmapped_action = "drop", multimapping_action = "drop", duplicate_action = "drop", 
    label = "ranges") 
}
\arguments{
  \item{ranges}{A named \code{GRanges} object.}
  \item{source_genome}{Declared source genome.}
  \item{target_genome}{Requested target genome.}
  \item{chain_file}{Directional UCSC chain file.}
  \item{unmapped_action}{Either \code{drop} or \code{error}.}
  \item{multimapping_action}{Either \code{drop} or \code{error}.}
  \item{duplicate_action}{How duplicate target intervals are handled.}
  \item{label}{Input label used in diagnostics.}
}
\details{
This is a standalone utility for explicit, user-controlled preprocessing. It is not called by the EPIC2 main analysis workflow. The main workflow instead requires \code{genome.build} and inputs already harmonized to that assembly.
}
\value{
A list containing lifted ranges, retained input indices, and QC counts.
}
\seealso{\code{\link{normalize_genome_build}}, \code{\link{validate_genome_coordinates}}.}
