\name{simplifyGO}
\alias{simplifyGO}
\title{
Simplify GO enrichment results
}
\description{
Simplify GO enrichment results
}
\usage{
simplifyGO(mat, method = "binary_cut", ...,
    
    # parameters for plotting the heatmap
    plot = TRUE, draw_word_cloud = TRUE, min_term = 5, order_by_size = FALSE,
    exclude_words = character(0), max_words = 10,
    word_cloud_width = 60)
}
\arguments{

  \item{mat}{A GO similarity matrix.}
  \item{method}{Method for clustering the matrix.}
  \item{...}{Other arguments passed to the clustering function, see details.}
  \item{plot}{Whether to make the heatmap.}
  \item{draw_word_cloud}{Whether draw the word clouds.}
  \item{min_term}{Minimal number of GO terms in a cluster.}
  \item{order_by_size}{Whether to reorder GO clusters by their sizes.}
  \item{exclude_words}{Words that are excluded in the word cloud.}
  \item{max_words}{Maximal number of words in the word cloud.}
  \item{word_cloud_width}{The maximal width of the viewport to put the word cloud.  The value should be numeric. It is measured in mm.}

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

Note the parametes for each clustering method are passes by \code{...} from \code{\link{simplifyGO}}.
}
\value{
A data frame with three columns:

\itemize{
  \item GO ID
  \item GO term name
  \item Cluster label
}
}
\examples{
\dontrun{
set.seed(123)
go_id = simplifyGO:::random_GO(500)
mat = GO_similarity(go_id)
df = simplifyGO(mat)
}
}
