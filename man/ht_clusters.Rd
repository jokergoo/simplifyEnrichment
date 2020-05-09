\name{ht_clusters}
\alias{ht_clusters}
\title{
Visualize the similarity matrix and the clustering
}
\description{
Visualize the similarity matrix and the clustering
}
\usage{
ht_clusters(mat, cl, dend = NULL,
    draw_word_cloud = TRUE, term = NULL, min_term = 5, order_by_size = FALSE,
    exclude_words = character(0), max_words = 10,
    word_cloud_grob_param = list(), fontsize_range = c(4, 16),
    ht_list = NULL, ...)
}
\arguments{

  \item{mat}{A similarity matrix.}
  \item{cl}{Cluster labels inferred from the similarity matrix, e.g. from \code{\link{cluster_terms}} or \code{\link{binary_cut}}.}
  \item{dend}{Used internally.}
  \item{draw_word_cloud}{Whether to draw the word clouds.}
  \item{term}{The full name or the description of the corresponding GO IDs. }
  \item{min_term}{Minimal number of functional terms in a cluster. All the clusters with size less than \code{min_term} are all merged into one separated cluster in the heatmap.}
  \item{order_by_size}{Whether to reorder clusters by their sizes. The cluster that is merged from small clusters (size < \code{min_term}) is always put to the bottom of the heatmap.}
  \item{exclude_words}{Words that are excluded in the word cloud.}
  \item{max_words}{Maximal number of words visualized in the word cloud.}
  \item{word_cloud_grob_param}{A list of graphic parameters passed to \code{\link{word_cloud_grob}}.}
  \item{fontsize_range}{The range of the font size. The value should be a numeric vector with length two. The minimal font size is mapped to word frequency value of 1 and the maximal font size is mapped to the maximal word frequency. The font size interlopation is linear.}
  \item{ht_list}{A list of additional heatmaps added to the left of the similarity heatmap.}
  \item{...}{Other arguments passed to \code{\link[ComplexHeatmap:draw-HeatmapList-method]{draw,HeatmapList-method}}.}

}
\value{
A \code{\link[ComplexHeatmap]{HeatmapList-class}} object.
}
\examples{
mat = readRDS(system.file("extdata", "similarity_mat.rds", package = "simplifyEnrichment"))
cl = binary_cut(mat)
ht_clusters(mat, cl, word_cloud_grob_param = list(max_width = 80))
}
