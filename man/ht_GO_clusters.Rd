\name{ht_GO_clusters}
\alias{ht_GO_clusters}
\title{
Visualize the GO similarity matrix and the classifications
}
\description{
Visualize the GO similarity matrix and the classifications
}
\usage{
ht_GO_clusters(mat, cl, dend = NULL,
    draw_word_cloud = TRUE, min_term = 5, order_by_size = FALSE,
    exclude_words = character(0), max_words = 10,
    word_cloud_grob_param = list(), fontsize_range = c(4, 16), ...)
}
\arguments{

  \item{mat}{A GO similarity matrix.}
  \item{cl}{Cluster labels inferred from the similarity matrix, e.g. from \code{\link{cluster_GO}} or \code{\link{binary_cut}}.}
  \item{dend}{Used internally.}
  \item{draw_word_cloud}{Whether to draw the word clouds.}
  \item{min_term}{Minimal number of GO terms in a cluster. All the clusters with size less than \code{min_term} are all merged into one single cluster in the heatmap.}
  \item{order_by_size}{Whether to reorder GO clusters by their sizes. The cluster that is merged from small clusters (size < 5) is always put to the bottom of the heatmap.}
  \item{exclude_words}{Words that are excluded in the word cloud.}
  \item{max_words}{Maximal number of words visualized in the word cloud.}
  \item{word_cloud_grob_param}{A list of parameters passed to \code{\link{word_cloud_grob}}.}
  \item{fontsize_range}{The range of the font size. The value should be a numeric vector with length two. The minimal font size is mapped to word frequency value of 1 and the maximal font size is mapped to the maximal word frequency. The font size interlopation is linear.}
  \item{...}{Other arguments passed to \code{\link[ComplexHeatmap:draw-HeatmapList-method]{draw,HeatmapList-method}}.}

}
\examples{
mat = readRDS(system.file("extdata", "similarity_mat.rds", package = "simplifyGO"))
cl = binary_cut(mat)
ht_GO_clusters(mat, cl)
}
