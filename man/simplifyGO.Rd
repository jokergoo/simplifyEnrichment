\name{simplifyGO}
\alias{simplifyGO}
\title{
Simplify Gene Ontology (GO) enrichment results
}
\description{
Simplify Gene Ontology (GO) enrichment results
}
\usage{
simplifyGO(mat, method = "binary_cut", control = list(),
    plot = TRUE, term = NULL, verbose = TRUE,
    column_title = qq("@{nrow(mat)} GO terms clustered by '@{method}'"),
    ht_list = NULL, ...)
}
\arguments{

  \item{mat}{A GO similarity matrix.}
  \item{method}{Method for clustering the matrix. See \code{\link{cluster_terms}}.}
  \item{control}{A list of parameters for controlling the clustering method, passed to \code{\link{cluster_terms}}.}
  \item{plot}{Whether to make the heatmap.}
  \item{term}{The full name or the description of the corresponding GO IDs. }
  \item{column_title}{Column title for the heatmap.}
  \item{verbose}{Whether to print messages.}
  \item{ht_list}{A list of additional heatmaps added to the left of the similarity heatmap.}
  \item{...}{Arguments passed to \code{\link{ht_clusters}}.}

}
\details{
This is basically a wrapper function that it first runs \code{\link{cluster_terms}} to cluster
GO terms and then runs \code{\link{ht_clusters}} to visualize the clustering.

The arguments in \code{\link{simplifyGO}} passed to \code{\link{ht_clusters}} are:

\describe{
  \item{\code{draw_word_cloud}}{Whether to draw the word clouds.}
  \item{\code{min_term}}{Minimal number of GO terms in a cluster. All the clusters with size less than \code{min_term} are all merged into one single cluster in the heatmap.}
  \item{\code{order_by_size}}{Whether to reorder GO clusters by their sizes. The cluster that is merged from small clusters (size < \code{min_term}) is always put to the bottom of the heatmap.}
  \item{\code{exclude_words}}{Words that are excluded in the word cloud.}
  \item{\code{max_words}}{Maximal number of words visualized in the word cloud.}
  \item{\code{word_cloud_grob_param}}{A list of graphic parameters passed to \code{\link{word_cloud_grob}}.}
  \item{\code{fontsize_range}}{The range of the font size. The value should be a numeric vector with length two. The minimal font size is mapped to word frequency value of 1 and the maximal font size is mapped to the maximal word frequency. The font size interlopation is linear.}
}
}
\value{
A data frame with three columns: GO IDs, GO term names and cluster labels.
}
\examples{
\donttest{
set.seed(123)
go_id = random_GO(500)
mat = GO_similarity(go_id)
df = simplifyGO(mat, word_cloud_grob_param = list(max_width = 80))
head(df)
}
}
