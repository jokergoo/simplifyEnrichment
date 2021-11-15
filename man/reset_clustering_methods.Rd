\name{reset_clustering_methods}
\alias{reset_clustering_methods}
\title{
Reset to default clustering methods
}
\description{
Reset to default clustering methods
}
\usage{
reset_clustering_methods()
}
\details{
The default methods are:

\describe{
  \item{\code{kmeans}}{see \code{\link{cluster_by_kmeans}}.}
  \item{\code{pam}}{see \code{\link{cluster_by_pam}}.}
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
No value is returned.
}
\examples{
all_clustering_methods()
remove_clustering_methods(c("kmeans", "mclust"))
all_clustering_methods()
reset_clustering_methods()
all_clustering_methods()
}
