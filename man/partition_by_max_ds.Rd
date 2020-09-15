\name{partition_by_max_ds}
\alias{partition_by_max_ds}
\title{
Partition by multiple methods
}
\description{
Partition by multiple methods
}
\usage{
partition_by_max_ds(mat, k)
}
\arguments{

  \item{mat}{The similarity matrix.}
  \item{k}{Number of clusters.}

}
\details{
It runs three partitionning methods: \code{\link{partition_by_pam}}, \code{\link{partition_by_kmeans}}
and \code{\link{partition_by_hclust}} and pick the classification with the highest
difference scores (calculated by \code{\link{difference_score}}).

This function is used to set to the \code{partition_fun} argument in \code{\link{binary_cut}}.
}
\examples{
# There is no example
NULL
}
