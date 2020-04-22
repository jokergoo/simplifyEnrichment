\name{cluster_by_mclust}
\alias{cluster_by_mclust}
\title{
Cluster GO similarity matrix by mclust
}
\description{
Cluster GO similarity matrix by mclust
}
\usage{
cluster_by_mclust(mat, ...)
}
\arguments{

  \item{mat}{The GO similarity matrix.}
  \item{...}{Other arguments passed to \code{\link[mclust]{Mclust}}.}

}
\details{
The value of \code{G} in \code{\link[mclust]{Mclust}} is set to \code{1:min(round(nrow(mat)/5), 100)}.
}
\value{
A vector of cluster labels (in numeric).
}
\examples{
# There is no example
NULL

}
