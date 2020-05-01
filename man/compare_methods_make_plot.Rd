\name{compare_methods_make_plot}
\alias{compare_methods_make_plot}
\title{
make plots for comparing clustering methods
}
\description{
make plots for comparing clustering methods
}
\usage{
compare_methods_make_plot(mat, clt, plot_type = c("mixed", "heatmap"))
}
\arguments{

  \item{mat}{A similarity matrix.}
  \item{clt}{A list of clusterings from \code{\link{compare_methods_make_clusters}}.}
  \item{plot_type}{What type of plots to make. See Details.}

}
\details{
If \code{plot_type} is the default value \code{mixed}, a figure with three panels generated:

\itemize{
  \item A heatmap of the similarity matrix with different classifications as row annotations.
  \item A heatmap of the pair-wise concordance of the classifications of every two clustering methods.
  \item Barplots of the difference scores for each method (calculated by \code{\link{difference_score}}), the number of clusters (all clusters and the clusters with size >= 5) and the mean similarity in the terms  that are in the same cluster.
}

If \code{plot_type} is \code{heatmap}. There are heatmaps for the similarity matrix under clusterings
from different methods. The last panel is a table with the number of clusters under different
clusterings.
}
\value{
No value is returned.
}
\examples{
# There is no example
NULL

}
