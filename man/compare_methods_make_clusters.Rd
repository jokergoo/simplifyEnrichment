\name{compare_methods_make_clusters}
\alias{compare_methods_make_clusters}
\title{
Apply various clustering methods
}
\description{
Apply various clustering methods
}
\usage{
compare_methods_make_clusters(mat, method = setdiff(ALL_CLUSTERING_METHODS, "mclust"),
    verbose = TRUE)
}
\arguments{

  \item{mat}{The similarity matrix.}
  \item{method}{Which methods to compare. All available methods are in \code{simplifyEnrichment:::ALL_CLUSTERING_METHODS}. A value of \code{all} takes all available methods. By default \code{mclust} is excluded because its long runtime.}
  \item{verbose}{Whether to print messages.}

}
\details{
The function compares following clustering methods:

\describe{
  \item{\code{binary_cut}}{see \code{\link{binary_cut}}.}
  \item{\code{kmeans}}{see \code{\link{cluster_by_kmeans}}.}
  \item{\code{dynamicTreeCut}}{see \code{\link{cluster_by_dynamicTreeCut}}.}
  \item{\code{mclust}}{see \code{\link{cluster_by_mclust}}. By default it is not included.}
  \item{\code{apcluster}}{see \code{\link{cluster_by_apcluster}}.}
  \item{\code{fast_greedy}}{see \code{\link{cluster_by_igraph}}.}
  \item{\code{leading_eigen}}{see \code{\link{cluster_by_igraph}}.}
  \item{\code{louvain}}{see \code{\link{cluster_by_igraph}}.}
  \item{\code{walktrap}}{see \code{\link{cluster_by_igraph}}.}
}
}
\value{
A list of cluster label vectors for different clustering methods.
}
\examples{
# There is no example
NULL
}
