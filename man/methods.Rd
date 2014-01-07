\name{Rising-methods}
\alias{print.Rising.fit}
\alias{plot.Rising.fit}
\alias{summary.Rising.fit}
%\alias{plot.lisrel}
%\alias{summary.lisrel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for Rising.fit objects
}
\description{
Print method prints the Rising.fit output , plot method plots the estimated network (with the \code{qgraph} package), and summary method returns density of the network, the value of gamma used, the rule used, and the time the analysis took.
}
\usage{
\method{print}{Rising.fit}(x, \dots)
\method{summary}{Rising.fit}(object, \dots)
\method{plot}{Rising.fit}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
output of \code{\link{Rising.fit}}
}
  \item{object}{
output of \code{\link{Rising.fit}}
}
  \item{\dots}{
  Arguments sent to qgraph. Only used in plot method.
}
}
\author{
Claudia van Borkulo
}