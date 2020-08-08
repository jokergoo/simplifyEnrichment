\name{cluster_terms}
\alias{cluster_terms}
\title{
Cluster functional terms
}
\description{
Cluster functional terms
}
\usage{
cluster_terms(mat, method = "binary_cut", control = list(), catch_error = FALSE,
    verbose = TRUE)
}
\arguments{

  \item{mat}{A similarity matrix.}
  \item{method}{Method for clustering the matrix.}
  \item{control}{A list of parameters passed to the corresponding clustering function.}
  \item{catch_error}{Internally used.}
  \item{verbose}{Whether to print messages.}

}
\details{
The following methods are the default:

\describe{
  \item{\code{kmeans}}{see \code{\link{cluster_by_kmeans}}.}
  \item{\code{dynamicTreeCut}}{see \code{\link{cluster_by_dynamicTreeCut}}.}
  \item{\code{mclust}}{see \code{\link{cluster_by_mclust}}.}
  \item{\code{apcluster}}{see \code{\link{cluster_by_apcluster}}.}
  \item{\code{hdbscan}}{see \code{\link{cluster_by_hdbscan}}.}
  \item{\code{fast_greedy}}{see \code{\link{cluster_by_igraph}}.}
  \item{\code{leading_eigen}}{see \code{\link{cluster_by_igraph}}.}
  \item{\code{louvain}}{see \code{\link{cluster_by_igraph}}.}
  \item{\code{walktrap}}{see \code{\link{cluster_by_igraph}}.}
  \item{\code{MCL}}{see \code{\link{cluster_by_MCL}}.}
  \item{\code{binary_cut}}{see \code{\link{binary_cut}}.}
}

Also the user-defined methods in \code{\link{all_clustering_methods}} can be used here.

New clustering methods can be registered by \code{\link{register_clustering_methods}}.

Please note it is better to directly call \code{\link{cluster_terms}} for clustering while not the individual \code{cluster_by_*} functions
because \code{\link{cluster_terms}} does additional cluster label adjustment.
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
