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
    word_cloud_grob_param = list(), fontsize_range = c(4, 16), value_range = NULL,
    bg_gp = gpar(fill = "#DDDDDD", col = "#AAAAAA"), side = c("right", "left"),
    add_new_line = FALSE, count_words_param = list(), ..., return_gbl = FALSE)
}
\arguments{

  \item{align_to}{How to align the annotations to the heatmap. Similar as in \code{\link[ComplexHeatmap]{anno_link}}, the value of \code{align_to} can be a list of row indices or a categorical vector where each vector in the list corresponds to a word cloud.  If it is a categorical vector, rows with the same level correspond to a same word cloud.  If \code{align_to} is a categorical vector and \code{term} is a list, names of \code{term} should have overlap to the levels in \code{align_to}. When \code{align_to} is set as a categorical vector, normally the same value is set to \code{row_split} in the main heatmap so that each row slice can correspond to a word cloud.}
  \item{term}{The description text used for constructing the word clouds. The value should have the same format as \code{align_to}. If \code{align_to} is a list, \code{term} should also be a list. In this case, the length of vectors in \code{term} is not necessarily the same as in \code{align_to}. E.g. \code{length(term[[1]])} is not necessarily equal to \code{length(align_to[[1]]}. If \code{align_to} is a categorical vector, \code{term} should also be a character vector with the same length as \code{align_to}. To make it more genrall, when \code{align_to} is a list, \code{term} can also be a list of data frames where the first column contains  keywords and the second column contains numeric values that will be mapped to font sizes in the word clouds.}
  \item{exclude_words}{The words excluced for construcing word cloud.}
  \item{max_words}{Maximal number of words visualized in the word cloud.}
  \item{word_cloud_grob_param}{A list of graphics parameters passed to \code{\link{word_cloud_grob}}.}
  \item{fontsize_range}{The range of the font size. The value should be a numeric vector with length two. The font size interpolation is linear.}
  \item{value_range}{The range of values to map to font sizes.}
  \item{bg_gp}{Graphics parameters for controlling the background.}
  \item{side}{Side of the annotation relative to the heatmap.}
  \item{add_new_line}{Whether to add new line after every word? If \code{TRUE}, each word will be in a separated line.}
  \item{count_words_param}{A list of parameters passed to \code{\link{count_words}}.}
  \item{...}{Other parameters.}
  \item{return_gbl}{Internally used.}

}
\details{
The word cloud annotation is constructed by \code{\link[ComplexHeatmap]{anno_link}}.

If the annotation is failed to construct or no keyword is found, the function returns a \code{\link[ComplexHeatmap]{anno_empty}} with 1px width.

English stop words, punctuation and numbers are removed by default when counting words. As specific stop words might
coincide with gene or pathway names, and numbers in genes names might be meaningful it is recommended to adjust this
behaviour by passing appropriate arguments to the \code{\link{count_words}} function using \code{count_words_param}.
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
