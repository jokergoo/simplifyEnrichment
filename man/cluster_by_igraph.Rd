\name{cluster_by_igraph}
\alias{cluster_by_igraph}
\title{
Cluster similarity matrix by graph community detection methods
}
\description{
Cluster similarity matrix by graph community detection methods
}
\usage{
cluster_by_igraph(mat,
    method = c("fast_greedy",
    "leading_eigen",
    "louvain",
    "walktrap"),
    ...)
}
\arguments{

  \item{mat}{The similarity matrix.}
  \item{method}{The community detection method.}
  \item{...}{Other arguments passed to the corresponding community detection function, see Details.}

}
\details{
The symmetric similarity matrix can be treated as an adjacency matrix and constructed as a graph/network with the similarity values as the weight of hte edges.
Thus, clustering the similarity matrix can be treated as detecting clusters/modules/communities from the graph.

Four methods implemented in igraph package can be used here:

\describe{
  \item{\code{fast_greedy}}{uses \code{\link[igraph]{cluster_fast_greedy}}.}
  \item{\code{leading_eigen}}{uses \code{\link[igraph]{cluster_leading_eigen}}.}
  \item{\code{louvain}}{uses \code{\link[igraph]{cluster_louvain}}.}
  \item{\code{walktrap}}{uses \code{\link[igraph]{cluster_walktrap}}.}
}
}
\value{
A vector of cluster labels (in numeric).
}
\examples{
# There is no example
NULL

}
