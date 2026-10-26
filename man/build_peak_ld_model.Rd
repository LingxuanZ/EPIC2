% Maintained against the exported function signature in R/ld.R.
\name{build_peak_ld_model}
\alias{build_peak_ld_model}
\title{Build peak-level GWAS response and local LD covariance}
\description{
Build peak-level GWAS response and local LD covariance.
}
\usage{
build_peak_ld_model(gwas, genotype, peak_matrix, overlap, window_bp = 5000000L, 
    ridge = 0.01, eigen_tolerance = 1e-08, max_snps_for_chr_ld = 15000L) 
}
\arguments{
  \item{gwas}{Harmonized GWAS table containing \code{coord_id} and \code{Z}.}
  \item{genotype}{Variant-by-sample genotype matrix named by \code{coord_id}.}
  \item{peak_matrix}{Peak-by-cell accessibility matrix.}
  \item{overlap}{Data frame with \code{snp} and \code{peak} columns.}
  \item{window_bp}{Local peak-pair distance.}
  \item{ridge}{Ridge used for within-peak and response covariance matrices.}
  \item{eigen_tolerance}{Eigenvalue floor.}
  \item{max_snps_for_chr_ld}{Maximum variants for caching a chromosome LD matrix.}
}
\value{
List containing peak response, covariance, and inverse square root.
}
