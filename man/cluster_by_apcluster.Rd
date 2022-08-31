\name{cluster_by_apcluster}
\alias{cluster_by_apcluster}
\title{
Cluster similarity matrix by apcluster
}
\description{
Cluster similarity matrix by apcluster
}
\usage{
cluster_by_apcluster(mat, s = apcluster::negDistMat(r = 2), ...)
}
\arguments{

  \item{mat}{The similarity matrix.}
  \item{s}{Passed to the \code{s} argument in \code{\link[apcluster]{apcluster}}.}
  \item{...}{Other arguments passed to \code{\link[apcluster]{apcluster}}.}

}
\value{
A vector of cluster labels (in numeric).
}
\examples{
# There is no example
NULL

}
