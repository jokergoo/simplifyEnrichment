\name{count_words}
\alias{count_words}
\title{
Calculate word frequency
}
\description{
Calculate word frequency
}
\usage{
count_words(term,
    exclude_words = NULL, stop_words = stopwords(),
    min_word_length = 1, tokenizer = 'words', transform_case = tolower,
    remove_numbers = TRUE, remove_punctuation = TRUE, custom_transformer = NULL,
    stemming = FALSE, dictionary = NULL)
}
\arguments{

  \item{term}{A vector of description texts.}
  \item{exclude_words}{The words that should be excluded.}
  \item{stop_words}{The stop words that should be be removed.}
  \item{min_word_length}{Minimum length of the word to be counted.}
  \item{tokenizer}{The tokenizer function, one of the values accepted by \code{tm::termFreq}.}
  \item{transform_case}{The function normalizing lettercase of the words.}
  \item{remove_numbers}{Whether to remove numbers.}
  \item{remove_punctuation}{Whether to remove punctuation.}
  \item{custom_transformer}{Custom function that transforms words.}
  \item{stemming}{Whether to only keep the roots of inflected words.}
  \item{dictionary}{A vector of words to be counted (if given all other words will be excluded).}

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
count_words(go_term)
}
