\name{all_clustering_methods}
\alias{all_clustering_methods}
\title{
All clustering methods
}
\description{
All clustering methods
}
\usage{
all_clustering_methods()
}
\details{
The default clustering methods are:

\describe{
  \item{\code{kmeans}}{see \code{\link{cluster_by_kmeans}}.}
  \item{\code{dynamicTreeCut}}{see \code{\link{cluster_by_dynamicTreeCut}}.}
  \item{\code{mclust}}{see \code{\link{cluster_by_mclust}}.}
  \item{\code{apcluster}}{see \code{\link{cluster_by_apcluster}}.}
  \item{\code{hdbscan}}{see \code{\link{cluster_by_hdbscan}}.}
  \item{\code{fast_greedy}}{see \code{\link{cluster_by_igraph}}.}
  \item{\code{louvain}}{see \code{\link{cluster_by_igraph}}.}
  \item{\code{walktrap}}{see \code{\link{cluster_by_igraph}}.}
  \item{\code{MCL}}{see \code{\link{cluster_by_MCL}}.}
  \item{\code{binary_cut}}{see \code{\link{binary_cut}}.}
}
}
\value{
A vector of method names.
}
\seealso{
New methods can be added by \code{\link{register_clustering_methods}}.
}
\examples{
all_clustering_methods()
}
