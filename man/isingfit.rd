\name{IsingFit}
\alias{IsingFit}
\title{
Network estimation using the eLasso method
}
\description{
This network estimation procedure eLasso, which is based on the Ising model, combines l1-regularized logistic regression with model selection based on the Extended Bayesian Information Criterion (EBIC). EBIC is a fit measure that identifies relevant relationships between variables. The resulting network consists of variables as nodes and relevant relationships as edges. Can deal with binary data.
}
\usage{
IsingFit(x, family='binomial', AND = TRUE, gamma = 0.25, 
plot = TRUE, progressbar = TRUE, lowerbound.lambda = NA,...)
}

\arguments{
  \item{x}{
Input matrix. The dimension of the matrix is nobs x nvars; each row is a vector of observations of the variables. Must be cross-sectional data.
}
\item{family}{
The default is 'binomial', treating the data as binary. Currently, this procedure is only supported for binary data.
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
Logical. Should the pbar be plotted in order to see the progress of the estimation procedure?
}
  \item{lowerbound.lambda}{
The minimum value of tuning parameter lambda (regularization parameter). Can be used to compare networks that are based on different sample sizes. The value is based on the number of observations in the smallest group n: lambda >= sqrt(log(p)/n). p is the number of variables. When both networks are estimated with the same lowerbound for lambda (based on the smallest group), the two networks can be directly compared.
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
\item{asymm.weights }{The (asymmetrical) weighted adjacency matrix before applying the AND/OR rule.}
\item{lambda.values }{The values of the tuning parameter per node that ensured the best fitting set of neighbors.}
}

\references{
Chen, J., & Chen, Z. (2008). Extended bayesian information criteria for model selection with large model spaces. Biometrika, 95(3), 759-771.

Foygel, R., & Drton, M. (2011). Bayesian model choice and information criteria in sparse generalized linear models. arXiv preprint arXiv:1112.5635.

Ravikumar, P., Wainwright, M. J., & Lafferty, J. D. (2010). High-dimensional Ising model selection using l1-regularized logistic regression. The Annals of Statistics, 38, 1287 - 1319.

van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F., Boschloo, L., Schoevers, R. A., & Waldorp, L. J. (2014). A new method for constructing networks from binary data. Scientific Reports 4, 5918; DOI:10.1038/srep05918. 
}
\author{
Claudia D. van Borkulo, Sacha Epskamp, with contributions from Alexander Robitzsch

Maintainer: Claudia D. van Borkulo <cvborkulo@gmail.com>
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

# Simulate:
Data <- IsingSampler(nSample, Graph, Thresh)

### Fit using IsingFit ###
Res <- IsingFit(Data, family='binomial', plot=FALSE)

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
