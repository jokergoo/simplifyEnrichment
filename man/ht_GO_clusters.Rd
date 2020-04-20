\name{ht_GO_clusters}
\alias{ht_GO_clusters}
\title{
Simplify GO enrichment results
}
\description{
Simplify GO enrichment results
}
\usage{
ht_GO_clusters(mat, cl, min_term = 5, order_by_size = TRUE,
    exclude_words = character(0), max_words = 10, ...)
}
\arguments{

  \item{mat}{A GO similarity matrix.}
  \item{cl}{Cluster labels.}
  \item{min_term}{Minimal number of GO terms in a cluster.}
  \item{order_by_size}{Whether to reorder GO clusters by their sizes.}
  \item{exclude_words}{Words that are excluded in the word cloud.}
  \item{max_words}{Maximal number of words in the word cloud.}
  \item{...}{Other arguments}

}
\examples{
# There is no example
NULL

}
