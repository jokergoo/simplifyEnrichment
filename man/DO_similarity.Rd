\name{DO_similarity}
\alias{DO_similarity}
\title{
Calculate Disease Ontology (DO) semantic similarity matrix
}
\description{
Calculate Disease Ontology (DO) semantic similarity matrix
}
\usage{
DO_similarity(do_id, measure = "Rel", remove_orphan_terms = FALSE)
}
\arguments{

  \item{do_id}{A vector of DO IDs.}
  \item{measure}{Semantic measure for the DO similarity, pass to \code{\link[DOSE]{doSim}}.}
  \item{remove_orphan_terms}{Whether to remove terms that have zero similarity to all other terms?}

}
\details{
This function is basically a wrapper on \code{\link[DOSE]{doSim}}.
}
\value{
A symmetric matrix.
}
\examples{
\donttest{
require(DOSE)
do_id = random_DO(10)
DO_similarity(do_id)
}
}
