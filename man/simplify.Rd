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
    order_by_size = TRUE, exclude_words = character(0))
}
\arguments{

  \item{mat}{A GO similarity matrix}
  \item{vavlue_fun}{-vavlue_fun}
  \item{cutoff}{cutoff}
  \item{plot}{Whether to make the plot}
  \item{min_term}{-min_term}
  \item{order_by_sizez}{-order_by_sizez}
  \item{exclude_words}{-exclude_words}

}
\examples{
\dontrun{
set.seed(123)
go_id = simplifyGO:::random_GO(500)
mat = GO_similarity(go_id)
lt = simplify(mat)
}
}
