% Maintained against the exported function signature in R/accessibility.R.
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
Each cell's feature counts are multiplied by the corresponding cell-scale value and summed within each group. The baseline is the unweighted mean of the group-specific feature vectors. When library-size adjustment is wanted, supply scales derived from the complete background ATAC counts, not from the trait-specific matrix alone.
}
\value{
A list with group-specific feature vectors and their unweighted mean baseline.
}
