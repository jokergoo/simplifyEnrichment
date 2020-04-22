\name{dend_node_apply}
\alias{dend_node_apply}
\title{
Apply on every node in a dendrogram
}
\description{
Apply on every node in a dendrogram
}
\usage{
dend_node_apply(dend, fun)
}
\arguments{

  \item{dend}{A dendrogram.}
  \item{fun}{A self-defined function.}

}
\details{
The function return a vector or a list with the same length as the number of nodes in the denrogram.

The self-defined function can have one single argument which is the sub-dendrogram at a certain node,
e.g. to get the number of members at each node:

  \preformatted{
    dend_node_apply(dend, function(d) attr(d, "members"))  }

The self-defined function can have a second argument, which is the index of current sub-dendrogram in 
the complete dendrogram. E.g. \code{dend[[1]]} is the first child node of the complete dendrogram and
\code{dend[[c(1, 2)]]} is the second child node of \code{dend[[1]]}, et al. This makes that at a certain node,
it is possible to get informatino of its children nodes and parent nodes.

  \preformatted{
    dend_node_apply(dend, function(d, index) \{
        d[[c(index, 1)]] # is the first child node of d
        d[[index[-length(index)]]] # is the parent node of d
        ...
    \})  }
}
\value{
A vector or a list, depends on whether \code{fun} returns a scalar or more complex values.
}
\examples{
# There is no example
NULL

}
