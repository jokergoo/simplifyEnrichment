\name{cluster_by_dynamicTreeCut}
\alias{cluster_by_dynamicTreeCut}
\title{
Cluster GO similarity matrix by dynamicTreeCut
}
\description{
Cluster GO similarity matrix by dynamicTreeCut
}
\usage{
cluster_by_dynamicTreeCut(mat, minClusterSize = 5, ...)
}
\arguments{

  \item{mat}{The GO similarity matrix.}
  \item{minClusterSize}{Minimal number of elements in a cluster. Pass to \code{\link[dynamicTreeCut]{cutreeDynamic}}.}
  \item{...}{Other arguments passed to \code{\link[dynamicTreeCut]{cutreeDynamic}}.}

}
\value{
A vector of cluster labels (in numeric).
}
\examples{
# There is no example
NULL

}
