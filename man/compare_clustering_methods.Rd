\name{compare_clustering_methods}
\alias{compare_clustering_methods}
\title{
Compare clustering methods
}
\description{
Compare clustering methods
}
\usage{
compare_clustering_methods(mat, method = setdiff(all_clustering_methods(), "mclust"),
    plot_type = c("mixed", "heatmap"), verbose = TRUE)
}
\arguments{

  \item{mat}{The similarity matrix.}
  \item{method}{Which methods to compare. All available methods are in \code{\link{all_clustering_methods}}. A value of \code{all} takes all available methods. By default \code{mclust} is excluded because its long runtime.}
  \item{plot_type}{See explanation in \code{\link{cmp_make_plot}}.}
  \item{verbose}{Whether to print messages.}

}
\details{
The function compares following clustering methods:

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

This functon is basically a wrapper function. It calls the following two functions:

\itemize{
  \item \code{\link{cmp_make_clusters}}: applies clustering by different methods.
  \item \code{\link{cmp_make_plot}}: makes the plots.
}
}
\value{
No value is returned.
}
\examples{
\dontrun{
mat = readRDS(system.file("extdata", "similarity_mat.rds", package = "simplifyEnrichment"))
compare_clustering_methods(mat)
compare_clustering_methods(mat, plot_type = "heatmap")
}
}
