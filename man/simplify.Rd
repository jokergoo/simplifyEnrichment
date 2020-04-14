\name{simplify}
\alias{simplify}
\title{
Simplify GO enrichment results
}
\description{
Simplify GO enrichment results
}
\usage{
simplify(mat, value_fun = median, cutoff = 0.8, plot = TRUE, min_term = 5,
    order_by_size = TRUE, exclude_words = character(0), max_words = 10)
}
\arguments{

  \item{mat}{A GO similarity matrix.}
  \item{value_fun}{Value function.}
  \item{cutoff}{Cutoff.}
  \item{plot}{Whether to make the heatmap.}
  \item{min_term}{Minimal number of GO terms in a cluster.}
  \item{order_by_size}{Whether to reorder GO clusters by their sizes.}
  \item{exclude_words}{Words that are excluded in the word cloud.}
  \item{max_words}{Maximal number of words in the word cloud.}

}
\examples{
\dontrun{
set.seed(123)
go_id = simplifyGO:::random_GO(500)
mat = GO_similarity(go_id)
lt = simplify(mat)
}
}
