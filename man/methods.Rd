\name{Ising-methods}
\alias{print.IsingFit}
\alias{plot.IsingFit}
\alias{summary.IsingFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for IsingFit objects
}
\description{
Print method prints the IsingFit output , plot method plots the estimated network (with the \code{qgraph} package), and summary method returns density of the network, the value of gamma used, the rule used, and the time the analysis took.
}
\usage{
\method{print}{IsingFit}(x, \dots)
\method{summary}{IsingFit}(object, \dots)
\method{plot}{IsingFit}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
output of \code{\link{IsingFit}}
}
  \item{object}{
output of \code{\link{IsingFit}}
}
  \item{\dots}{
  Arguments sent to qgraph. Only used in plot method.
}
}
\author{
Claudia van Borkulo
}