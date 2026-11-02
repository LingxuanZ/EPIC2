% Maintained against the exported function signature in R/model.R.
\name{combine_component_tests}
\alias{combine_component_tests}
\title{Combine SNP/peak count and Wilcoxon component tests with ACAT}
\description{
Combine SNP/peak count and Wilcoxon component tests with ACAT.
}
\usage{
combine_component_tests(snp_accessibility, snp_wrs, peak_accessibility, peak_wrs, 
    require_all = FALSE) 
}
\arguments{
  \item{snp_accessibility}{SNP accessibility coefficient table.}
  \item{snp_wrs}{SNP Wilcoxon coefficient table.}
  \item{peak_accessibility}{Peak accessibility coefficient table.}
  \item{peak_wrs}{Peak Wilcoxon coefficient table.}
  \item{require_all}{Require all four valid component P-values.}
}
\details{
Tables are matched by \code{trait} and \code{group}; each table must have a \code{p} column. By default all available valid P-values in [0, 1] are combined with equal weights. With \code{require_all = TRUE}, all four components must be valid. The reported component count uses this same validity rule.
}
\value{
Trait-by-group component and combined P-value table.
}
