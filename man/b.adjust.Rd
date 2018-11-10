\name{b.adjust}
\alias{b.adjust}
\title{Beta-adjust Function}
\usage{
  b.adjust(p, ntests, fixParam = c(0, 0, 0))
}
\arguments{
  \item{p}{a vector of the minimum P-values}

  \item{ntests}{the numbers of tests for the set of
  P-values}

  \item{fixParam}{an indicator vector of 0 or 1 (length =
  3)}
}
\description{
  Given a set of minimum P-values under the different
  numbers of tests, returns p-values adjusted using the
  Beta correction that allows you to control the family
  wise error rates.
}
\examples{
p=NULL
for(i in 1:100){p = c(p, min(runif(i)))}
b.adjust(p, 1:100)
}

