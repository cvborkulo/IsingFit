\name{Rising.fit}
\alias{Rising.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Network estimation using the eLasso method
}
\description{
This network estimation procedure combines l1-regularized logistic regression with model selection based on the Extended Bayesian Information Criterion to identify relevant symptom-symptom relationships. These relationships define connections in a network.
}
\usage{
Rising.fit(x, AND = TRUE, gamma = 0.25, plot = TRUE, progressbar = TRUE, ...)
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
Rising.fit returns (invisibly) a 'Rising.fit' object that contains the following items:
\item{weiadj }{The weighted adjacency matrix.}
\item{thresholds }{Thresholds of the variables.}
\item{q }{The object that is returned by qgraph (class 'qgraph').}
\item{gamma }{The value of hyperparameter gamma.}
\item{AND }{A logical indicating whether the AND-rule is used or not. If not, the OR-rule is used.}
\item{time }{The time it took to estimate the network.}
}
\references{
Van Borkulo, C.D., Borsboom, D., Epskamp, S., Blanken, T., Bosschloo, L., Schoevers, R. A., & Waldorp, L. J. (2013). Fitting like a glove: A new method for constructing networks for psychometric data. Manuscript in preparation.
}
\author{
Claudia van Borkulo <cvborkulo@gmail.com>
}
\note{
I should have an awesome website, but I haven't. Sorry about that!
}


\seealso{
qgraph IsingSampler
}
\examples{



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
