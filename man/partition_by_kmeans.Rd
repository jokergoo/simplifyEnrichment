\name{partition_by_kmeans}
\alias{partition_by_kmeans}
\title{
Partition by kmeans
}
\description{
Partition by kmeans
}
\usage{
partition_by_kmeans(mat, k, n_repeats = 10)
}
\arguments{

  \item{mat}{The similarity matrix.}
  \item{k}{Number of clusters.}
  \item{n_repeats}{Number of repeated run of k-means.}

}
\details{
Since k-means clustering brings randomness, this function performs
k-means clustering several times and uses the final consensus clustering.
}
\examples{
# There is no example
NULL

}
