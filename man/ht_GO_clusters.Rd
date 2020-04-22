\name{ht_GO_clusters}
\alias{ht_GO_clusters}
\title{
Simplify GO enrichment results
}
\description{
Simplify GO enrichment results
}
\usage{
ht_GO_clusters(mat, cl, dend = NULL,
    draw_word_cloud = TRUE, min_term = 5, order_by_size = FALSE,
    exclude_words = character(0), max_words = 10, word_cloud_width = 60, ...)
}
\arguments{

  \item{mat}{A GO similarity matrix.}
  \item{cl}{Cluster labels.}
  \item{dend}{Used internally.}
  \item{draw_word_cloud}{Whether draw the word clouds.}
  \item{min_term}{Minimal number of GO terms in a cluster.}
  \item{order_by_size}{Whether to reorder GO clusters by their sizes.}
  \item{exclude_words}{Words that are excluded in the word cloud.}
  \item{max_words}{Maximal number of words in the word cloud.}
  \item{word_cloud_width}{The maximal width of the viewport to put the word cloud.  The value should be numeric. It is measured in mm.}
  \item{...}{Other arguments}

}
\examples{
# There is no example
NULL

}
