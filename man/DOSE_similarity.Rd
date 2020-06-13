\name{DOSE_similarity}
\alias{DOSE_similarity}
\title{
Calculate DOSE similarity matrix
}
\description{
Calculate DOSE similarity matrix
}
\usage{
DOSE_similarity(do_id, measure = "Rel")
}
\arguments{

  \item{do_id}{A vector of DOSE IDs.}
  \item{measure}{Measurement for the DOSE similarity, pass to \code{\link[GOSemSim]{termSim}}.}

}
\details{
This function is basically a wrapper on \code{\link[GOSemSim]{termSim}}.
}
\value{
A symmetric matrix.
}
\examples{
# There is no example
NULL

}
