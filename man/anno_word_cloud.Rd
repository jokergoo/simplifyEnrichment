\name{anno_word_cloud}
\alias{anno_word_cloud}
\title{
Word cloud annotations
}
\description{
Word cloud annotations
}
\usage{
anno_word_cloud(align_to, term, exclude_words = NULL, max_words = 10,
    word_cloud_grob_param = list(), fontsize_range = c(4, 16),
    bg_gp = gpar(fill = "#DDDDDD", col = "#AAAAAA"), side = c("right", "left"), ...)
}
\arguments{

  \item{align_to}{How to align the annotations to the heatmap. Similar as in \code{\link[ComplexHeatmap]{anno_link}}, the value of \code{align_to} can be a list of row indices or a categorical vector where each vector in the list corresponds to a word cloud.  If it is a categorical vector, rows with the same level correspond to a same word cloud.  If \code{align_to} is a categorical vector and \code{term} is a list, names of \code{term} should have overlap to the levels in \code{align_to}. When \code{align_to} is set as a categorical vector, normally the same value is set to \code{row_split} in the main heatmap so that each row slice can correspond to a word cloud.}
  \item{term}{The description text used for constructing the word clouds. The value should have the same format as \code{align_to}. If \code{align_to} is a list, \code{term} should also be a list. In this case, the length of vectors in \code{term} is not necessarily the same as in \code{align_to}. E.g. \code{length(term[[1]])} is not necessarily equal to \code{length(align_to[[1]]}. If \code{align_to} is a categorical vector, \code{term} should also be a character vector with the same length as \code{align_to}.}
  \item{exclude_words}{The words excluced for construcing word cloud.}
  \item{max_words}{Maximal number of words visualized in the word cloud.}
  \item{word_cloud_grob_param}{A list of graphics parameters passed to \code{\link{word_cloud_grob}}.}
  \item{fontsize_range}{The range of the font size. The value should be a numeric vector with length two. The minimal font size is mapped to word frequency value of 1 and the maximal font size is mapped to the maximal word frequency. The font size interlopation is linear.}
  \item{bg_gp}{Graphics parameters for controlling the background.}
  \item{side}{Side of the annotation relative to the heatmap.}
  \item{...}{Other parameters.}

}
\details{
The word cloud annotation is constructed by \code{\link[ComplexHeatmap]{anno_link}}.

If the annotation is failed to construct or no keyword is found, the function returns a \code{\link[ComplexHeatmap]{anno_empty}} with 1px width.
}
\examples{
gm = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds", package = "simplifyEnrichment"))
go_id = rownames(gm)
go_term = AnnotationDbi::select(GO.db::GO.db, keys = go_id, columns = "TERM")$TERM

split = sample(letters[1:4], 100, replace = TRUE)
align_to = split(1:100, split)
term = lapply(letters[1:4], function(x) sample(go_term, sample(100:400, 1)))
names(term) = letters[1:4]

require(ComplexHeatmap)
mat = matrix(rnorm(100*10), nrow = 100)
Heatmap(mat, cluster_rows = FALSE, row_split = split, 
	right_annotation = rowAnnotation(foo = anno_word_cloud(align_to, term)))
}
