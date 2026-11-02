% Maintained against the exported function signature in R/atac.R.
\name{add_blacklist_ratio}
\alias{add_blacklist_ratio}
\title{Add blacklist-overlap ratios}
\description{
Count each fragment record at most once across all blacklist intervals, using the same cell mappings and total-fragment denominator as \code{add_frip()}. Read multiplicity is ignored. Cells without observed fragments receive an NA ratio.
}
\usage{
add_blacklist_ratio(atac, blacklist = NULL, assay = "ATAC", genome = NULL) 
}
\arguments{
  \item{atac}{A Seurat object with fragment objects.}
  \item{blacklist}{A \code{GRanges} blacklist. Defaults to the Signac blacklist for the explicitly declared genome.}
  \item{assay}{ATAC assay name.}
  \item{genome}{ATAC genome build used to select the default blacklist.}
}
\value{
The input object with blacklist count and ratio metadata columns.
}
