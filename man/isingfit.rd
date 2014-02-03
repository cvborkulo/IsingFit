\name{IsingFit}
\alias{IsingFit}
\title{
Network estimation using the eLasso method
}
\description{
This network estimation procedure eLasso, which is based on the Ising model, combines l1-regularized logistic regression with model selection based on the Extended Bayesian Information Criterion (EBIC). EBIC is a fit measure that identifies relevant relationships between variables. The resulting network consists of variables as nodes and relevant relationships as edges.
}
\usage{
IsingFit(x, AND = TRUE, gamma = 0.25, plot = TRUE, progressbar = TRUE, ...)
}

\arguments{
  \item{x}{
A binary matrix.
}
  \item{AND}{
Logical. Can be TRUE of FALSE to indicate whether the AND-rule or the OR-rule should be used to define the edges in the network. Defaults to TRUE.
}
  \item{gamma}{
A value of hyperparameter gamma in the extended BIC. Can be anything between 0 and 1. Defaults to .25.
}
  \item{plot}{
Logical. Should the resulting network be plotted?
}
  \item{progressbar}{
Logical. Should the progressbar be plotted in order to see the progress of the estimation procedure?
}
\item{\dots}{
Arguments sent to \code{qgraph}.
}
}

\value{
IsingFit returns (invisibly) a 'IsingFit' object that contains the following items:
\item{weiadj }{The weighted adjacency matrix.}
\item{thresholds }{Thresholds of the variables.}
\item{q }{The object that is returned by qgraph (class 'qgraph').}
\item{gamma }{The value of hyperparameter gamma.}
\item{AND }{A logical indicating whether the AND-rule is used or not. If not, the OR-rule is used.}
\item{time }{The time it took to estimate the network.}
}
\references{
Van Borkulo, C.D., Borsboom, D., Epskamp, S., Blanken, T., Bosschloo, L., Schoevers, R. A., & Waldorp, L. J. (2013). Fitting like a glove: A new method for constructing networks for psychometric data. Manuscript submitted for publication.
}
\author{
Claudia van Borkulo <cvborkulo@gmail.com>
}
\note{
See also my website: http://cvborkulo.com
}

\examples{
library("IsingSampler")

### Simulate dataset ###
# Input:
N <- 6 # Number of nodes
nSample <- 1000 # Number of samples

# Ising parameters:
Graph <- matrix(sample(0:1,N^2,TRUE,prob = c(0.8, 0.2)),N,N) * runif(N^2,0.5,2)
Graph <- pmax(Graph,t(Graph))
diag(Graph) <- 0
Thresh <- -rowSums(Graph) / 2

# Siumlate:
Data <- IsingSampler(nSample, Graph, Thresh)

### Fit using IsingFit ###
Res <- IsingFit(Data, plot=FALSE)

# Plot results:
library("qgraph")
layout(t(1:2))
qgraph(Res$weiadj,fade = FALSE)
title("Estimated network")
qgraph(Graph,fade = FALSE)
title("Original network")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
