\name{compare_methods}
\alias{compare_methods}
\title{
Compare clustering methods
}
\description{
Compare clustering methods
}
\usage{
compare_methods(...)
}
\arguments{

  \item{...}{All pass to \code{\link{compare_clustering_methods}}.}

}
\details{
This function will be removed soon.
}
\value{
No value is returned.

}
\examples{
\dontrun{
mat = readRDS(system.file("extdata", "similarity_mat.rds", package = "simplifyEnrichment"))
compare_methods(mat)
compare_methods(mat, plot_type = "heatmap")
}

}
