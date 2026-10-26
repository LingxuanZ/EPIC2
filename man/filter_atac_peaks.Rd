% Maintained against the exported function signature in R/atac.R.
\name{filter_atac_peaks}
\alias{filter_atac_peaks}
\title{Filter ATAC peaks using manuscript QC thresholds}
\description{
Filter ATAC peaks using manuscript QC thresholds.
}
\usage{
filter_atac_peaks(atac, thresholds = list(), assay = "ATAC", blacklist = NULL, 
    genome = NULL) 
}
\arguments{
  \item{atac}{A Seurat/Signac object.}
  \item{thresholds}{Named list containing peak width thresholds.}
  \item{assay}{ATAC assay name.}
  \item{blacklist}{A \code{GRanges} blacklist. Defaults to the Signac blacklist for the explicitly declared genome.}
  \item{genome}{ATAC genome build used to select the default blacklist.}
}
\value{
The filtered object.
}
