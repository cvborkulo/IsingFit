\name{Rising-methods}
\alias{print.Risingfit}
\alias{plot.Risingfit}
\alias{summary.Risingfit}
%\alias{plot.lisrel}
%\alias{summary.lisrel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for Risingfit objects
}
\description{
Print method prints the Risingfit output , plot method plots the estimated network (with the \code{qgraph} package), and summary method returns density of the network, the value of gamma used, the rule used, and the time the analysis took.
}
\usage{
\method{print}{Risingfit}(x, \dots)
\method{summary}{Risingfit}(object, \dots)
\method{plot}{Risingfit}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
output of \code{\link{Risingfit}}
}
  \item{object}{
output of \code{\link{Risingfit}}
}
  \item{\dots}{
  Arguments sent to qgraph. Only used in plot method.
}
}
\author{
Claudia van Borkulo
}