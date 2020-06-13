\name{cluster_by_mclust}
\alias{cluster_by_mclust}
\title{
Cluster similarity matrix by mclust
}
\description{
Cluster similarity matrix by mclust
}
\usage{
cluster_by_mclust(mat, G = seq_len(max(2, min(round(nrow(mat)/5), 100))), ...)
}
\arguments{

  \item{mat}{The similarity matrix.}
  \item{G}{Passed to the \code{G} argument in \code{\link[mclust]{Mclust}}.}
  \item{...}{Other arguments passed to \code{\link[mclust]{Mclust}}.}

}
\value{
A vector of cluster labels (in numeric).
}
\examples{
# There is no example
NULL
}
