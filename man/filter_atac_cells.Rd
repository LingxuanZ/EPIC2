% Maintained against the exported function signature in R/atac.R.
\name{filter_atac_cells}
\alias{filter_atac_cells}
\title{Filter nuclei using manuscript QC thresholds}
\description{
Filter nuclei using manuscript QC thresholds.
}
\usage{
filter_atac_cells(atac, thresholds = list(), group_column = "tissue_cell_type") 
}
\arguments{
  \item{atac}{A Seurat/Signac object.}
  \item{thresholds}{Named list of cell-level thresholds.}
  \item{group_column}{Metadata column used for minimum group size filtering.}
}
\value{
The filtered object.
}
