\name{cluster_by_igraph}
\alias{cluster_by_igraph}
\title{
Cluster GO similarity matrix by community detection methods
}
\description{
Cluster GO similarity matrix by community detection methods
}
\usage{
cluster_by_igraph(mat,
    method = c("cluster_fast_greedy",
    "cluster_leading_eigen",
    "cluster_louvain",
    "cluster_walktrap"),
    ...)
}
\arguments{

  \item{mat}{The Go similarity matrix.}
  \item{method}{method}
  \item{...}{Other arguments passed to.}

}
\examples{
# There is no example
NULL

}
