\name{cluster_by_kmeans}
\alias{cluster_by_kmeans}
\title{
Cluster similarity matrix by k-means clustering
}
\description{
Cluster similarity matrix by k-means clustering
}
\usage{
cluster_by_kmeans(mat, max_k = max(2, min(round(nrow(mat)/5), 100)), ...)
}
\arguments{

  \item{mat}{The similarity matrix.}
  \item{max_k}{maximal k for k-means clustering.}
  \item{...}{Other arguments passed to \code{\link[stats]{kmeans}}.}

}
\details{
The best number of k for k-means clustering is identified according to the "elbow" or "knee" method on
the distribution of within-cluster sum of squares at each k.
}
\value{
A vector of cluster labels (in numeric).
}
\examples{
# There is no example
NULL
}
