\name{compare_methods}
\alias{compare_methods}
\title{
Compare clustering methods
}
\description{
Compare clustering methods
}
\usage{
compare_methods(mat, method = "all")
}
\arguments{

  \item{mat}{The similarity matrix.}
  \item{method}{Which methods to compare. All available methods are in \code{simplifyEnrichment:::ALL_CLUSTERING_METHODS}. A value of \code{all} takes all available methods.}

}
\details{
The function compares following clustering methods:

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

The function produces a plot with three panels.

\itemize{
  \item A heatmap of the similarity matrix with different classifications as row annotations.
  \item A barplot of the difference scores for each method, calculated by \code{\link{difference_score}}.
  \item A heatmap of the pair-wise concordance of the classifications of every two clustering methods.
}
}
\value{
No value is returned.
}
\examples{
\dontrun{
mat = readRDS(system.file("extdata", "similarity_mat.rds", package = "simplifyEnrichment"))
# here 'mclust' is removed because it needs very long time to run
compare_methods(mat, method = setdiff(simplifyEnrichment:::ALL_CLUSTERING_METHODS, "mclust"))
}
}
