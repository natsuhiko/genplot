\name{badp}
\alias{badp}
\title{Beta-adjusted P-values}
\usage{
  badp(p, ntests, badj)
}
\arguments{
  \item{p}{a vector of the minimum P-values}

  \item{ntests}{the numbers of tests for the set of
  P-values}
}
\description{
  Given the parameter estimate by b.adjust, returns
  p-values adjusted.
}
\examples{
# not run
badp(p, 1:100, badj)
}

