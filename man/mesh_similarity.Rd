\name{mesh_similarity}
\alias{mesh_similarity}
\title{
Calculate MeSH similarity matrix
}
\description{
Calculate MeSH similarity matrix
}
\usage{
mesh_similarity(mesh_id, measure = "Rel")
}
\arguments{

  \item{mesh_id}{A vector of MeSH IDs.}
  \item{measure}{Measurement for the MeSH similarity, pass to \code{\link[meshes]{meshSim}}.}

}
\details{
This function is basically a wrapper on \code{\link[meshes]{meshSim}}.
}
\value{
A symmetric matrix.
}
\examples{
\donttest{
mesh_id = random_mesh(10)
mesh_similarity(mesh_id)
}
}
