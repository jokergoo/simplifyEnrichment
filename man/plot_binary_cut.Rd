\name{plot_binary_cut}
\alias{plot_binary_cut}
\title{
Visualize the process of binary cut
}
\description{
Visualize the process of binary cut
}
\usage{
plot_binary_cut(mat, value_fun = median, cutoff = 0.85,
    partition_fun = partition_by_pam, dend = NULL, dend_width = unit(3, "cm"),
    depth = NULL, show_heatmap_legend = TRUE, ...)
}
\arguments{

  \item{mat}{The similarity matrix.}
  \item{value_fun}{Value function to calculate the score for each node in the dendrogram.}
  \item{cutoff}{The cutoff for splitting the dendrogram.}
  \item{partition_fun}{A function to split each node into two groups. Pre-defined functions in this package are \code{\link{partition_by_kmeans}}, \code{\link{partition_by_pam}} (the default) and \code{\link{partition_by_hclust}}.}
  \item{dend}{A dendrogram object, used internally.}
  \item{depth}{Depth of the recursive binary cut process.}
  \item{dend_width}{Width of the dendrogram.}
  \item{show_heatmap_legend}{Whether to show the heatmap legend.}
  \item{...}{Other arguments.}

}
\details{
After the functions which performs clustering are executed, such as \code{\link{simplifyGO}} or
\code{\link{binary_cut}}, the dendrogram is temporarily saved and \code{\link{plot_binary_cut}} directly
uses this dendrogram. So, if the partition function brings randomness, it makes sure
the clustering is the same as the one made by e.g. \code{\link{simplifyGO}}.
}
\examples{
\donttest{
mat = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds", 
    package = "simplifyEnrichment"))
plot_binary_cut(mat, depth = 1)
plot_binary_cut(mat, depth = 2)
plot_binary_cut(mat)
}
}
