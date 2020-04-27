\name{difference_score}
\alias{difference_score}
\title{
Difference score
}
\description{
Difference score
}
\usage{
difference_score(mat, cl)
}
\arguments{

  \item{mat}{The similarity matrix.}
  \item{cl}{Cluster labels.}

}
\details{
This function measures the different between the similarity values for the terms
that belong to the same clusters and in different clusters. The difference score
is the Kolmogorov-Smirnov statistic between the two distributions.
}
\examples{
mat = readRDS(system.file("extdata", "similarity_mat.rds", package = "simplifyEnrichment"))
cl = binary_cut(mat)
difference_score(mat, cl)
}
