\name{GO_similarity}
\alias{GO_similarity}
\title{
Calculate GO similarity matrix
}
\description{
Calculate GO similarity matrix
}
\usage{
GO_similarity(go_id, ont, db = 'org.Hs.eg.db', measure = "Rel")
}
\arguments{

  \item{go_id}{A vector of GO IDs.}
  \item{ont}{GO ontology, value should be one of "BP", "CC" or "MF". If it is not specified, the function automatically identifies it by random sampling 10 IDs from \code{go_id} (see \code{\link{guess_ont}}).}
  \item{db}{Annotation database. It should be from \url{https://bioconductor.org/packages/3.10/BiocViews.html#___OrgDb}}
  \item{measure}{Measurement for the GO similarity, pass to \code{\link[GOSemSim]{mgoSim}}.}

}
\details{
This function is basically a wrapper on \code{\link[GOSemSim]{mgoSim}}.
}
\value{
A symmetric matrix.
}
\examples{
# There is no example
NULL

}
