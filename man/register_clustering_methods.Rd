\name{register_clustering_methods}
\alias{register_clustering_methods}
\title{
Register new clustering methods
}
\description{
Register new clustering methods
}
\usage{
register_clustering_methods(...)
}
\arguments{

  \item{...}{A named list of clustering functions, see Details}

}
\details{
The user-defined functions should accept at least one argument which is the input matrix. 
The second optional argument should always be \code{...} so that parameters
for the clustering function can be passed by \code{control} argument from \code{\link{cluster_terms}}, \code{\link{simplifyGO}} or \code{\link{simplifyEnrichment}}.
If users forget to add \code{...}, it is added internally.

Please note, the user-defined function should automatically identify the optimized
number of clusters.

The function should return a vector of cluster labels. Internally it is converted to numeric labels.
}
\value{
No value is returned.
}
\examples{
register_clustering_methods(
	# assume there are 5 groups
	random = function(mat, ...) sample(5, nrow(mat), replace = TRUE)
)
all_clustering_methods()
remove_clustering_methods("random")
}
