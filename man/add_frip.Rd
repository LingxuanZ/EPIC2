% Maintained against the exported function signature in R/atac.R.
\name{add_frip}
\alias{add_frip}
\title{Add fraction-of-reads-in-peaks metrics}
\description{
Add fraction-of-reads-in-peaks metrics.
}
\usage{
add_frip(atac, assay = "ATAC", genome_sep = c("-", "-")) 
}
\arguments{
  \item{atac}{A Seurat object with Signac fragment objects attached.}
  \item{assay}{ATAC assay name.}
  \item{genome_sep}{Two separators used in peak row names.}
}
\value{
The input object with fragment totals and \code{FRIP} metadata columns.
}
