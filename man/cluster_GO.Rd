\name{cluster_GO}
\alias{cluster_GO}
\title{
Cluster GO terms
}
\description{
Cluster GO terms
}
\usage{
cluster_GO(mat, method = "binary_cut", catch_error = FALSE, verbose = TRUE, ...)
}
\arguments{

  \item{mat}{A GO similarity matrix.}
  \item{method}{Method for clustering the matrix.}
  \item{catch_error}{Internally used.}
  \item{verbose}{Whether print messages.}
  \item{...}{Other arguments passed to the clustering function.}

}
\value{
A numeric vector of cluster labels.

If \code{catch_error} is set to \code{TRUE} and if the clustering produces an error,
the function returns a \code{try-error} object.
}
\examples{
# There is no example
NULL

}
