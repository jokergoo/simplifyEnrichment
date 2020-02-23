\name{simplify}
\alias{simplify}
\title{
Simplify GO enrichment results
}
\description{
Simplify GO enrichment results
}
\usage{
simplify(mat, cutoff = 0.8, plot = TRUE)
}
\arguments{

  \item{mat}{A GO similarity matrix}
  \item{cutoff}{cutoff}
  \item{plot}{Whether to make the plot}

}
\examples{
\dontrun{
set.seed(123)
go_id = simplifyGO:::random_GO(500)
mat = GO_similarity(go_id)
lt = simplify(mat)
}
}
