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
This function measures the different between the similarity values for the GO terms
that belong to the same clusters and in different clusters. The difference score
is the Kolmogorov-Smirnov statistic between the two distributions.
}
\examples{
# There is no example
NULL

}
