\name{count_word}
\alias{count_word}
\title{
Calculate word frequency
}
\description{
Calculate word frequency
}
\usage{
count_word(go_id, term = NULL, exclude_words = NULL)
}
\arguments{

  \item{go_id}{A vector of GO IDs.}
  \item{term}{The corresponding names or description of terms if the input are not GO terms.}
  \item{exclude_words}{The words that should be excluded.}

}
\details{
The input can be simply set with a vector of GO id to \code{go_id} argument that the GO names
are automatically extracted. If the input are not GO terms, users need to provide a vector of long names/descriptions by \code{term} argument.

If the input is GO id, the following words are excluded: \code{c("via", "protein", "factor", "side", "type", "specific")}.
They are analyzed by \code{simplifyEnrichment:::all_GO_word_count()}.

The text preprocessing followings the instructions from \url{http://www.sthda.com/english/wiki/word-cloud-generator-in-r-one-killer-function-to-do-everything-you-need} .
}
\value{
A data frame with words and frequencies.
}
\examples{
\donttest{
go_id = random_GO(500)
head(count_word(go_id))
}
}
