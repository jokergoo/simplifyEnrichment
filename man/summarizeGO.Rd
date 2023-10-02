\name{summarizeGO}
\alias{summarizeGO}
\title{
A simplified way to visualize enrichment in GO clusters
}
\description{
A simplified way to visualize enrichment in GO clusters
}
\usage{
summarizeGO(go_id, value = NULL, aggregate = mean,
    method = "binary_cut", control = list(), verbose = TRUE,
    axis_label = "Value", title = "", legend_title = axis_label,
    
    min_term = round(nrow(mat)*0.01),
    stat = "pvalue",
    min_stat = ifelse(stat == "count", 5, 0.05),
    exclude_words = character(0),
    max_words = 6,
    word_cloud_grob_param = list(),
    fontsize_range = c(4, 16),
    bg_gp = gpar(fill = "#DDDDDD", col = "#AAAAAA")
    )
}
\arguments{

  \item{go_id}{A vector of GO IDs.}
  \item{value}{A list of numeric value associate with \code{go_id}. We suggest to use -log10(p.adjust) or -log2(fold enrichment) as the values.}
  \item{aggregate}{Function to aggregate values in each GO cluster.}
  \item{method}{Method for clustering the matrix. See \code{\link{cluster_terms}}.}
  \item{control}{A list of parameters for controlling the clustering method, passed to \code{\link{cluster_terms}}.}
  \item{verbose}{Whether to print messages.}
  \item{axis_label}{X-axis label.}
  \item{title}{Title for the whole plot.}
  \item{legend_title}{Title for the legend.}
  \item{min_term}{Minimal number of functional terms in a cluster. All the clusters with size less than \code{min_term} are all merged into one separated cluster in the heatmap.}
  \item{order_by_size}{Whether to reorder clusters by their sizes. The cluster that is merged from small clusters (size < \code{min_term}) is always put to the bottom of the heatmap.}
  \item{stat}{Type of value for mapping to the font size of keywords in the word clouds. There are two options: "count": simply number of keywords; "pvalue": enrichment on keywords is performed (by fisher's exact test) and -log10(pvalue) is used to map to font sizes.}
  \item{min_stat}{Minimal value for \code{stat} for selecting keywords.}
  \item{exclude_words}{Words that are excluded in the word cloud.}
  \item{max_words}{Maximal number of words visualized in the word cloud.}
  \item{word_cloud_grob_param}{A list of graphic parameters passed to \code{\link{word_cloud_grob}}.}
  \item{fontsize_range}{The range of the font size. The value should be a numeric vector with length two. The font size interpolation is linear.}
  \item{bg_gp}{Graphics parameters for controlling word cloud annotation background.}

}
\details{
There are several other ways to specify GO IDs and the associated values.

1. specify \code{value} as a named vector where GO IDs are the names.
2. specify \code{value} as a list of numeric named vectors. In this case, \code{value} contains multiple enrichment results.
}
\examples{
# There is no example
NULL
}
