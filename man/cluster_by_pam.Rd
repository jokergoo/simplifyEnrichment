\name{cluster_by_pam}
\alias{cluster_by_pam}
\title{
Cluster similarity matrix by pam clustering
}
\description{
Cluster similarity matrix by pam clustering
}
\usage{
cluster_by_pam(mat, max_k = max(2, min(round(nrow(mat)/10), 100)), ...)
}
\arguments{

  \item{mat}{The similarity matrix.}
  \item{max_k}{maximal k for pam clustering.}
  \item{...}{Other arguments passed to \code{\link[fpc]{pamk}}.}

}
\details{
PAM is applied by \code{\link[fpc]{pamk}} which can automatically select the best k.
}
\value{
A vector of cluster labels (in numeric).
}
\examples{
# There is no example
NULL

}
