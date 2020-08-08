\name{cmp_make_plot}
\alias{cmp_make_plot}
\title{
Make plots for comparing clustering methods
}
\description{
Make plots for comparing clustering methods
}
\usage{
cmp_make_plot(mat, clt, plot_type = c("mixed", "heatmap"), nrow = 3)
}
\arguments{

  \item{mat}{A similarity matrix.}
  \item{clt}{A list of clusterings from \code{\link{cmp_make_clusters}}.}
  \item{plot_type}{What type of plots to make. See Details.}
  \item{nrow}{Number of rows of the layout when \code{plot_type} is set to \code{heatmap}.}

}
\details{
If \code{plot_type} is the default value \code{mixed}, a figure with three panels generated:

\itemize{
  \item A heatmap of the similarity matrix with different classifications as row annotations.
  \item A heatmap of the pair-wise concordance of the classifications of every two clustering methods.
  \item Barplots of the difference scores for each method (calculated by \code{\link{difference_score}}), the number of clusters (total clusters and the clusters with size >= 5) and the mean similarity of the terms  that are in the same clusters.
}

If \code{plot_type} is \code{heatmap}. There are heatmaps for the similarity matrix under clusterings
from different methods. The last panel is a table with the number of clusters under different
clusterings.
}
\value{
No value is returned.
}
\examples{
\dontrun{
mat = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds",
    package = "simplifyEnrichment"))
clt = cmp_make_clusters(mat)
cmp_make_plot(mat, clt)
cmp_make_plot(mat, clt, plot_type = "heatmap")
}
}
