\name{export_to_shiny_app}
\alias{export_to_shiny_app}
\title{
Interactively visualize the similarity heatmap
}
\description{
Interactively visualize the similarity heatmap
}
\usage{
export_to_shiny_app(mat, cl = binary_cut(mat))
}
\arguments{

  \item{mat}{A similarity matrix.}
  \item{cl}{Cluster labels inferred from the similarity matrix, e.g. from \code{\link{cluster_terms}} or \code{\link{binary_cut}}.}

}
\examples{
mat = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds", package = "simplifyEnrichment"))
cl = binary_cut(mat)
if(interactive()) {
    export_to_shiny_app(mat, cl)
}
}
