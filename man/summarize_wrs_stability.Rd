% Maintained against the exported function signature in R/wrs.R.
\name{summarize_wrs_stability}
\alias{summarize_wrs_stability}
\title{Summarize repeat-to-repeat WRS feature stability}
\description{
Summarize repeat-to-repeat WRS feature stability.
}
\usage{
summarize_wrs_stability(result, top_n = 100L, auc_cutoff = 0.55, stability_frequency_cutoff = 0.8) 
}
\arguments{
  \item{result}{Raw repeated balanced WRS table.}
  \item{top_n}{Number of top one-sided features compared between repeats.}
  \item{auc_cutoff}{AUC threshold defining selected features.}
  \item{stability_frequency_cutoff}{Minimum selected-repeat proportion for a stable feature.}
}
\value{
List containing pairwise repeat and per-feature stability tables.
}
