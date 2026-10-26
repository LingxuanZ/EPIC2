% Maintained against the exported function signature in R/model.R.
\name{adjust_pvalues_within_trait}
\alias{adjust_pvalues_within_trait}
\title{Add within-trait multiple-testing corrections}
\description{
Add within-trait multiple-testing corrections.
}
\usage{
adjust_pvalues_within_trait(table, p_column = "p") 
}
\arguments{
  \item{table}{Result data frame containing trait and P-value columns.}
  \item{p_column}{Name of the unadjusted P-value column.}
}
\value{
Data frame with BH and Bonferroni columns.
}
