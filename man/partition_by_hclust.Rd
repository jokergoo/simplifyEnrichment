\name{partition_by_hclust}
\alias{partition_by_hclust}
\title{
Partition by hclust
}
\description{
Partition by hclust
}
\usage{
partition_by_hclust(mat, k)
}
\arguments{

  \item{mat}{The similarity matrix.}
  \item{k}{Number of clusters.}

}
\details{
The clusering method and the distance method both take the defaults of \code{\link[stats]{hclust}} and \code{\link[stats]{dist}}.

This function is used to set to the \code{partition_fun} argument in \code{\link{binary_cut}}.
}
\examples{
# There is no example
NULL

}
