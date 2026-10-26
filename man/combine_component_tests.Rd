% Maintained against the exported function signature in R/model.R.
\name{combine_component_tests}
\alias{combine_component_tests}
\title{Combine SNP/peak count and Wilcoxon component tests with ACAT}
\description{
Combine SNP/peak count and Wilcoxon component tests with ACAT.
}
\usage{
combine_component_tests(snp_accessibility, snp_wrs, peak_accessibility, peak_wrs, 
    require_all = TRUE) 
}
\arguments{
  \item{snp_accessibility}{SNP accessibility coefficient table.}
  \item{snp_wrs}{SNP Wilcoxon coefficient table.}
  \item{peak_accessibility}{Peak accessibility coefficient table.}
  \item{peak_wrs}{Peak Wilcoxon coefficient table.}
  \item{require_all}{Require all four valid component P-values.}
}
\details{
Tables are matched by \code{trait} and \code{group}; each table must have a \code{p} column. With \code{require_all = TRUE}, a combined result requires valid P-values for all four components. Otherwise available valid components can be combined.
}
\value{
Trait-by-group component and combined P-value table.
}
