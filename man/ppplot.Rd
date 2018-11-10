\name{ppplot}
\alias{ppplot}
\title{Plotting QQ-plot of P-values}
\usage{
  ppplot(p, col=1, add=F, ylim=c(0,10))
}
\arguments{
  \item{p}{a vector of P-values}

}
\description{
  Visual confirmation of P-values under uniform assumption.
}
\examples{
# not run
ppplot(badp(p, 1:100, badj))
}

