\name{anno_word_cloud_from_GO}
\alias{anno_word_cloud_from_GO}
\title{
Word cloud annotations from GO
}
\description{
Word cloud annotations from GO
}
\usage{
anno_word_cloud_from_GO(align_to, go_id, term = NULL, exclude_words = NULL, ...)
}
\arguments{

  \item{align_to}{The same format as in \code{\link{anno_word_cloud}}.}
  \item{go_id}{The value should be in the same format as \code{align_to}. If \code{go_id} is a vector, it should have the same length as \code{align_to}, and if \code{go_id} is a list, note, e.g. \code{length(go_id[[1]])} is not necessarily equal to \code{length(align_to[[1]]}. If \code{align_to} is a categorical vector and \code{go_id} is a list, names of \code{go_id} should have overlap to the levels in \code{align_to}.}
  \item{term}{Alternatively the GO description can be set via the \code{term} argument. The same format as in \code{\link{anno_word_cloud}}.}
  \item{exclude_words}{The words excluced for construcing word cloud. Some words are internally exclucded: \code{c("via", "protein", "factor", "side", "type", "specific")}.}
  \item{...}{All other arguments passed to \code{\link{anno_word_cloud}}.}

}
\examples{
# There is no example
NULL

}
