\name{binary_cut}
\alias{binary_cut}
\title{
Cluster functional terms by recursively binary cutting the similarity matrix
}
\description{
Cluster functional terms by recursively binary cutting the similarity matrix
}
\usage{
binary_cut(mat, value_fun = median, partition_fun = partition_by_pam,
    cutoff = 0.85)
}
\arguments{

  \item{mat}{A similarity matrix.}
  \item{value_fun}{Value function to calculate the score for each node in the dendrogram.}
  \item{partition_fun}{A function to split each node into two groups. Pre-defined functions in this package are \code{\link{partition_by_kmeans}}, \code{\link{partition_by_pam}} (the default) and \code{\link{partition_by_hclust}}.}
  \item{cutoff}{The cutoff for splitting the dendrogram.}

}
\value{
A vector of cluster labels (in numeric).
}
\examples{
mat = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds",
    package = "simplifyEnrichment"))
binary_cut(mat)
}
