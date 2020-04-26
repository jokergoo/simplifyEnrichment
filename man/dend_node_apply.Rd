\name{dend_node_apply}
\alias{dend_node_apply}
\title{
Apply functions on every node in a dendrogram
}
\description{
Apply functions on every node in a dendrogram
}
\usage{
dend_node_apply(dend, fun)
}
\arguments{

  \item{dend}{A dendrogram.}
  \item{fun}{A self-defined function.}

}
\details{
The function return a vector or a list as the same length as the number of nodes in the dendrogram.

The self-defined function can have one single argument which is the sub-dendrogram at a certain node.
E.g. to get the number of members at every node:

  \preformatted{
    dend_node_apply(dend, function(d) attr(d, "members"))  }

The self-defined function can have a second argument, which is the index of current sub-dendrogram in 
the complete dendrogram. E.g. \code{dend[[1]]} is the first child node of the complete dendrogram and
\code{dend[[c(1, 2)]]} is the second child node of \code{dend[[1]]}, et al. This makes that at a certain node,
it is possible to get information of its child nodes and parent nodes.

  \preformatted{
    dend_node_apply(dend, function(d, index) \{
        d[[c(index, 1)]] # is the first child node of d
        d[[index[length(index)-1]]] # is the parent node of d
        ...
    \})  }

Node for the top node, the value of \code{index} is \code{NULL}.
}
\value{
A vector or a list, depends on whether \code{fun} returns a scalar or more complex values.
}
\examples{
mat = matrix(rnorm(100), 10)
dend = as.dendrogram(hclust(dist(mat)))
# number of members on every node
dend_node_apply(dend, function(d) attr(d, "members"))
# the depth on every node
dend_node_apply(dend, function(d, index) length(index))
}
