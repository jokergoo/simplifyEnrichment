\name{cluster_by_hdbscan}
\alias{cluster_by_hdbscan}
\title{
Cluster similarity matrix by hdbscan
}
\description{
Cluster similarity matrix by hdbscan
}
\usage{
cluster_by_hdbscan(mat, minPts = 5, ...)
}
\arguments{

  \item{mat}{The similarity matrix.}
  \item{minPts}{Passed to the \code{minPts} argument in \code{\link[dbscan]{hdbscan}}.}
  \item{...}{Other arguments passed to \code{\link[dbscan]{hdbscan}}.}

}
\value{
A vector of cluster labels (in numeric).
}
\examples{
# There is no example
NULL

}
