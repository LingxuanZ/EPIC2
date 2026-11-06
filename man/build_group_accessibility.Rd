\name{build_group_accessibility}
\alias{build_group_accessibility}
\title{Aggregate library-size-adjusted accessibility by cell group}
\description{
Aggregate library-size-adjusted accessibility by cell group.
}
\usage{
build_group_accessibility(feature_matrix, groups, cell_scale = NULL)
}
\arguments{
  \item{feature_matrix}{Sparse feature-by-cell matrix.}
  \item{groups}{Named group vector aligned to matrix columns.}
  \item{cell_scale}{Optional named per-cell scaling vector. Defaults to one.}
}
\details{
Each cell's feature counts are multiplied by the corresponding cell-scale value.
Each group total is divided by that group's
retained cell count. The baseline is the total across all retained cells divided
by their total number, equivalently the cell-count-weighted mean of group means.
Cells with missing groups or nonfinite/nonpositive scaling factors are excluded.
Zero feature counts remain in the cell-mean denominators.
When library-size adjustment is wanted, supply scales derived from the complete
background ATAC counts, not from the trait-specific matrix alone.
}
\value{
A list with within-group cell-mean feature vectors and a pooled-cell mean baseline vector.
}
