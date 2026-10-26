% Maintained against the exported function signature in R/ld.R.
\name{compute_snp_ld_covariance}
\alias{compute_snp_ld_covariance}
\title{Construct local SNP-level LD response covariance}
\description{
Construct local SNP-level LD response covariance.
}
\usage{
compute_snp_ld_covariance(genotype, gwas, columns = EPIC2_defaults()$gwas$columns, 
    window_bp = 5000000L, factor_two = FALSE, correlation_block_size = 256L) 
}
\arguments{
  \item{genotype}{Variant-by-sample dosage matrix named by genomic coordinate.}
  \item{gwas}{Harmonized GWAS table containing \code{coord_id}, chromosome, and position.}
  \item{columns}{GWAS column mapping.}
  \item{window_bp}{Maximum same-chromosome distance represented in the covariance.}
  \item{factor_two}{Use \code{2 * R^2} instead of the response correlation \code{R^2}.}
  \item{correlation_block_size}{Number of focal variants handled in each BLAS cross-product.}
}
\details{
Rows and columns are sorted in genomic order. Off-diagonal entries are retained only for variants on the same chromosome within \code{window_bp}. With \code{factor_two = TRUE}, the diagonal is 2 and off-diagonal entries are \eqn{2r^2}. Otherwise the diagonal is 1 and entries are \eqn{r^2}. This is response covariance for squared GWAS statistics, rather than a signed genotype-correlation matrix.
}
\value{
A sparse symmetric covariance matrix ordered by chromosome and position.
}
