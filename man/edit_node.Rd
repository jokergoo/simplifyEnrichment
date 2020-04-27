\name{edit_node}
\alias{edit_node}
\title{
Modify nodes in a dendrogram
}
\description{
Modify nodes in a dendrogram
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
which makes it possible to get information of child nodes and parent nodes for
a specified node.

As an example, we first assign random values to every node in the dendrogram:

  \preformatted{
    mat = matrix(rnorm(100), 10)
    dend = as.dendrogram(hclust(dist(mat)))
    dend = edit_node(dend, function(d) \{attr(d, 'score') = runif(1); d\})  }

Then for every node, we take the maximal absolute difference to all its child nodes
and parent node as the attribute \code{abs_diff}

  \preformatted{
    dend = edit_node(dend, function(d, index) \{
        n = length(index)
        s = attr(d, "score")
        if(is.null(index)) \{  # d is the top node
            s_children = sapply(d, function(x) attr(x, "score"))
            s_parent = NULL
        \} else if(is.leaf(d)) \{ # d is the leaf
            s_children = NULL
            s_parent = attr(dend[[index[-n]]], "score")
        \} else \{
            s_children = sapply(d, function(x) attr(x, "score"))
            s_parent = attr(dend[[index[-n]]], "score")
        \}
        abs_diff = max(abs(s - c(s_children, s_parent)))
        attr(d, "abs_diff") = abs_diff
        return(d)
    \})  }
}
\value{
A dendrogram object.
}
\examples{
# There is no example
NULL

}
