% Maintained against the exported function signature in R/model.R.
\name{fit_whitened_groups}
\alias{fit_whitened_groups}
\title{Fit group-wise LD-whitened linear models}
\description{
To reproduce the notebook, the response and complete raw design are first
transformed by the inverse square-root covariance matrix. The transformed
intercept column is then discarded, and the fit uses an ordinary intercept
together with the transformed focal and baseline predictors.
}
\usage{
fit_whitened_groups(response, predictors, baseline, inverse_sqrt, trait = NA_character_, 
    predictor_name = "predictor") 
}
\arguments{
  \item{response}{Named response vector.}
  \item{predictors}{Named list of group-specific predictor vectors.}
  \item{baseline}{Named baseline predictor vector.}
  \item{inverse_sqrt}{Named inverse square-root covariance matrix.}
  \item{trait}{Trait label included in output.}
  \item{predictor_name}{Name of the focal model term.}
}
\details{
All vectors and covariance dimensions are matched using feature identifiers. A single coefficient row is returned for the focal predictor in each group, including diagnostic status for a model that cannot be fit.
}
\value{
One coefficient row per group.
}
