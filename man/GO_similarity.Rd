\name{GO_similarity}
\alias{GO_similarity}
\title{
Calculate Gene Ontology (GO) semantic similarity matrix
}
\description{
Calculate Gene Ontology (GO) semantic similarity matrix
}
\usage{
GO_similarity(go_id, ont, db = 'org.Hs.eg.db', measure = "Rel")
}
\arguments{

  \item{go_id}{A vector of GO IDs.}
  \item{ont}{GO ontology. Value should be one of "BP", "CC" or "MF". If it is not specified, the function automatically identifies it by random sampling 10 IDs from \code{go_id} (see \code{\link{guess_ont}}).}
  \item{db}{Annotation database. It should be from \url{https://bioconductor.org/packages/3.10/BiocViews.html#___OrgDb}}
  \item{measure}{Semantic measure for the GO similarity, pass to \code{\link[GOSemSim]{termSim}}.}

}
\details{
This function is basically a wrapper on \code{\link[GOSemSim]{termSim}}.
}
\value{
A symmetric matrix.
}
\examples{
\donttest{
go_id = random_GO(100)
mat = GO_similarity(go_id)
}
}
