\name{make_one_sided_wrs}
\alias{make_one_sided_wrs}
\title{Convert two-sided Wilcoxon results to enrichment-oriented one-sided scores}
\description{
Convert two-sided Wilcoxon results to enrichment-oriented one-sided scores.
}
\usage{
make_one_sided_wrs(result, truncation_quantiles = c(0.02, 0.98)) 
}
\arguments{
  \item{result}{Data frame returned by \code{presto::wilcoxauc()}.}
  \item{truncation_quantiles}{Lower and upper empirical winsorization quantiles. Cutoffs are estimated from finite raw normal scores across all feature--group pairs. Infinite scores from P = 0 or 1 are capped at those cutoffs.}
}
\value{
Input data frame with one-sided P-value and Z-score columns.
}
