\name{edit_node}
\alias{edit_node}
\title{
Modify every node in a dendrogram
}
\description{
Modify every node in a dendrogram
}
\usage{
edit_node(dend, fun = function(d, index) d)
}
\arguments{

  \item{dend}{A dendrogram.}
  \item{fun}{A self-defined function.}

}
\details{
if \code{fun} only has one argument, it is basically the same as \code{\link[stats]{dendrapply}},
but it can have a second argument which is the index of the node in the dendrogram,
which makes it possible to get information of children nodes and parent nodes for
a specified node.
}
\examples{
# There is no example
NULL

}
