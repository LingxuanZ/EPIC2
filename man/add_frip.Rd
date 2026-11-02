\name{add_frip}
\alias{add_frip}
\title{Add fraction-of-fragments-in-peaks metrics}
\description{
Count each fragment-file record once in the denominator and once in the
numerator if it overlaps at least one retained peak. The fifth-column read
multiplicity is not used. Fragment files are streamed in bounded blocks, and
raw barcodes are mapped to Seurat cell names separately for each fragment object.
Fragment coordinates use the zero-based, half-open BED convention; peak ranges
are one-based and inclusive.
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
Cells without fragment records have \code{NA} FRIP.
}
