\name{count_word}
\alias{count_word}
\title{
Calculate word frequency
}
\description{
Calculate word frequency
}
\usage{
count_word(term, exclude_words = NULL)
}
\arguments{

  \item{term}{A vector of description texts.}
  \item{exclude_words}{The words that should be excluded.}

}
\details{
The text preprocessing followings the instructions from \url{http://www.sthda.com/english/wiki/word-cloud-generator-in-r-one-killer-function-to-do-everything-you-need} .
}
\value{
A data frame with words and frequencies.
}
\examples{
gm = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds", package = "simplifyEnrichment"))
go_id = rownames(gm)
go_term = AnnotationDbi::select(GO.db::GO.db, keys = go_id, columns = "TERM")$TERM
count_word(go_term)
}
