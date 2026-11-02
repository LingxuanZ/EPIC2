% Maintained against the exported function signature in R/ld.R.
\name{compute_inverse_sqrt}
\alias{compute_inverse_sqrt}
\title{Compute a ridge-stabilized inverse square root}
\description{
Compute a ridge-stabilized inverse square root.
}
\usage{
compute_inverse_sqrt(matrix, ridge = 0.01, eigen_tolerance = 1e-08,
    max_dense_block_size = 5000L)
}
\arguments{
  \item{matrix}{Symmetric square matrix.}
  \item{ridge}{Nonnegative ridge added to every eigenvalue.}
  \item{eigen_tolerance}{Lower eigenvalue floor.}
  \item{max_dense_block_size}{Maximum dense block dimension checked before conversion or decomposition. A positive integer, or \code{Inf} to deliberately disable the guard.}
}
\details{
The eigendecomposition is used to form the inverse square root after flooring eigenvalues at \code{eigen_tolerance} and adding \code{ridge}. Matrix row and column identifiers are retained.

The default guard rejects dimensions above 5000. Errors report one-matrix memory and an illustrative six-matrix estimate; actual peak memory can be higher. Raising the limit does not change the calculation. This is an allocation guard, not a total-memory bound or a scalable eigensolver.
}
\value{
Dense symmetric inverse square-root matrix.
}
\examples{
x <- diag(2)
dimnames(x) <- list(c("a", "b"), c("a", "b"))
compute_inverse_sqrt(x)
}
