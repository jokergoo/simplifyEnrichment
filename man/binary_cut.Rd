\name{binary_cut}
\alias{binary_cut}
\title{
Cluster GO terms by binary cutting the similarity matrix
}
\description{
Cluster GO terms by binary cutting the similarity matrix
}
\usage{
binary_cut(mat, value_fun = median, cutoff = 0.85, n_run = 1)
}
\arguments{

  \item{mat}{A GO similarity matrix.}
  \item{value_fun}{Value function to calculate the score for each node in the dendrogram.}
  \item{cutoff}{The cutoff for splitting the dendrogram.}
  \item{n_run}{If the value is larger than one, \code{\link{binary_cut}} is executed multiple times and generates a consensus clustering.}

}
\value{
A vector of cluster labels (in numeric).
}
\examples{
mat = readRDS(system.file("extdata", "similarity_mat.rds", package = "simplifyGO"))
binary_cut(mat)
}
