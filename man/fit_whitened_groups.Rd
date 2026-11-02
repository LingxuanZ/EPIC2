% Maintained against the exported function signature in R/model.R.
\name{fit_whitened_groups}
\alias{fit_whitened_groups}
\title{Fit group-wise LD-whitened linear models}
\description{
The response and complete raw design, including the intercept, are
transformed by the inverse square-root covariance matrix. The fit retains
all transformed design columns to preserve the original-scale linear model.
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
The transformed intercept is retained. Coefficients, residual variance and QR-based standard errors use the same transformed design. Missing or nonfinite inputs for any covariance feature produce an unavailable test with diagnostic status; construct a new covariance and inverse square root for a reduced feature set rather than slicing an existing whitening matrix.
}
\value{
One coefficient row per group.
}
