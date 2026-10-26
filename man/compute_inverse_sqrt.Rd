% Maintained against the exported function signature in R/ld.R.
\name{compute_inverse_sqrt}
\alias{compute_inverse_sqrt}
\title{Compute a ridge-stabilized inverse square root}
\description{
Compute a ridge-stabilized inverse square root.
}
\usage{
compute_inverse_sqrt(matrix, ridge = 0.01, eigen_tolerance = 1e-08) 
}
\arguments{
  \item{matrix}{Symmetric square matrix.}
  \item{ridge}{Nonnegative ridge added to every eigenvalue.}
  \item{eigen_tolerance}{Lower eigenvalue floor.}
}
\details{
The eigendecomposition is used to form the inverse square root after flooring eigenvalues at \code{eigen_tolerance} and adding \code{ridge}. Matrix row and column identifiers are retained.
}
\value{
Dense symmetric inverse square-root matrix.
}
\examples{
x <- diag(2)
dimnames(x) <- list(c("a", "b"), c("a", "b"))
compute_inverse_sqrt(x)
}
