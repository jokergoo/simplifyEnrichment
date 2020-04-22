\name{cluster_by_kmeans}
\alias{cluster_by_kmeans}
\title{
Cluster GO similarity matrix by k-means clustering
}
\description{
Cluster GO similarity matrix by k-means clustering
}
\usage{
cluster_by_kmeans(mat, ...)
}
\arguments{

  \item{mat}{The GO similarity matrix.}
  \item{...}{Other arguments passed to \code{\link[stats]{kmeans}}.}

}
\details{
The number of clusters are tried from 2 to \code{min(round(nrow(mat)/5), 100)}. The best number
of k for k-means clustering is identified according to the "elbow" or "knee" method on
the distribution of within-cluster sum of squares at each k.
}
\value{
A vector of cluster labels (in numeric).
}
\examples{
# There is no example
NULL

}
