\name{mesh_similarity}
\alias{mesh_similarity}
\title{
Calculate mesh similarity matrix
}
\description{
Calculate mesh similarity matrix
}
\usage{
mesh_similarity(mesh_id, measure = "Rel")
}
\arguments{

  \item{mesh_id}{A vector of mesh IDs.}
  \item{measure}{Measurement for the mesh similarity, pass to \code{\link[GOSemSim]{termSim}}.}

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
