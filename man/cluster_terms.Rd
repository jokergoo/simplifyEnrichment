\name{cluster_terms}
\alias{cluster_terms}
\title{
Cluster functional terms
}
\description{
Cluster functional terms
}
\usage{
cluster_terms(mat, method = "binary_cut", catch_error = FALSE, verbose = TRUE, ...)
}
\arguments{

  \item{mat}{A similarity matrix.}
  \item{method}{Method for clustering the matrix.}
  \item{catch_error}{Internally used.}
  \item{verbose}{Whether to print messages.}
  \item{...}{Other arguments passed to the clustering function.}

}
\details{
The following methods are supported:

\describe{
  \item{\code{binary_cut}}{see \code{\link{binary_cut}}.}
  \item{\code{kmeans}}{see \code{\link{cluster_by_kmeans}}.}
  \item{\code{dynamicTreeCut}}{see \code{\link{cluster_by_dynamicTreeCut}}.}
  \item{\code{mclust}}{see \code{\link{cluster_by_mclust}}.}
  \item{\code{apcluster}}{see \code{\link{cluster_by_apcluster}}.}
  \item{\code{fast_greedy}}{see \code{\link{cluster_by_igraph}}.}
  \item{\code{leading_eigen}}{see \code{\link{cluster_by_igraph}}.}
  \item{\code{louvain}}{see \code{\link{cluster_by_igraph}}.}
  \item{\code{walktrap}}{see \code{\link{cluster_by_igraph}}.}
}

Note the parametes for each clustering method are passes by \code{...} from \code{\link{cluster_terms}}.
}
\value{
A numeric vector of cluster labels (in numeric).

If \code{catch_error} is set to \code{TRUE} and if the clustering produces an error,
the function returns a \code{try-error} object.
}
\examples{
# There is no example
NULL

}
