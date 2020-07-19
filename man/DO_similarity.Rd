\name{DO_similarity}
\alias{DO_similarity}
\title{
Calculate Disease Ontology (DO) similarity matrix
}
\description{
Calculate Disease Ontology (DO) similarity matrix
}
\usage{
DO_similarity(do_id, measure = "Rel")
}
\arguments{

  \item{do_id}{A vector of DO IDs.}
  \item{measure}{Measurement for the DO similarity, pass to \code{\link[DOSE]{doSim}}.}

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
