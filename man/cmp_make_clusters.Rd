\name{cmp_make_clusters}
\alias{cmp_make_clusters}
\title{
Apply various clustering methods
}
\description{
Apply various clustering methods
}
\usage{
cmp_make_clusters(mat, method = setdiff(all_clustering_methods(), "mclust"),
    verbose = TRUE)
}
\arguments{

  \item{mat}{The similarity matrix.}
  \item{method}{Which methods to compare. All available methods are in \code{\link{all_clustering_methods}}. A value of \code{all} takes all available methods. By default \code{mclust} is excluded because its long runtime.}
  \item{verbose}{Whether to print messages.}

}
\details{
The function compares following default clustering methods by default:

\describe{
  \item{\code{kmeans}}{see \code{\link{cluster_by_kmeans}}.}
  \item{\code{dynamicTreeCut}}{see \code{\link{cluster_by_dynamicTreeCut}}.}
  \item{\code{mclust}}{see \code{\link{cluster_by_mclust}}. By default it is not included.}
  \item{\code{apcluster}}{see \code{\link{cluster_by_apcluster}}.}
  \item{\code{hdbscan}}{see \code{\link{cluster_by_hdbscan}}.}
  \item{\code{fast_greedy}}{see \code{\link{cluster_by_igraph}}.}
  \item{\code{leading_eigen}}{see \code{\link{cluster_by_igraph}}.}
  \item{\code{louvain}}{see \code{\link{cluster_by_igraph}}.}
  \item{\code{walktrap}}{see \code{\link{cluster_by_igraph}}.}
  \item{\code{MCL}}{see \code{\link{cluster_by_MCL}}.}
  \item{\code{binary_cut}}{see \code{\link{binary_cut}}.}
}

Also the user-defined methods in \code{\link{all_clustering_methods}} are also compared.
}
\value{
A list of cluster label vectors for different clustering methods.
}
\examples{
\dontrun{
mat = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds",
    package = "simplifyEnrichment"))
clt = cmp_make_clusters(mat)
}
}
