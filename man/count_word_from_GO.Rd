\name{count_word_from_GO}
\alias{count_word_from_GO}
\title{
Calculate word frequency From GO
}
\description{
Calculate word frequency From GO
}
\usage{
count_word_from_GO(go_id, term = NULL, exclude_words = NULL)
}
\arguments{

  \item{go_id}{A vector of GO IDs.}
  \item{term}{The corresponding names or description of terms if the input are not GO terms.}
  \item{exclude_words}{The words that should be excluded.}

}
\details{
The input can be simply set with a vector of GO id to \code{go_id} argument so that the GO names
are automatically extracted. Users can also provide a vector of long names/descriptions by \code{term} argument.

If the input is GO id, the following words are excluded: \code{c("via", "protein", "factor", "side", "type", "specific")}.
They are analyzed by \code{simplifyEnrichment:::all_GO_word_count()}.
}
\value{
A data frame with words and frequencies.
}
\seealso{
\code{\link{count_word}}
}
\examples{
gm = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds", package = "simplifyEnrichment"))
go_id = rownames(gm)
head(count_word_from_GO(go_id))
}
