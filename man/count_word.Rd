\name{count_word}
\alias{count_word}
\title{
Calculate word cloud
}
\description{
Calculate word cloud
}
\usage{
count_word(go_id, term = NULL, exclude_words = NULL)
}
\arguments{

  \item{go_id}{A vector of GO IDs.}
  \item{term}{The corresponding GO terms. If provided, it will speed up this function.}
  \item{exclude_words}{The words that should be excluded.}

}
\details{
The text preprocessing followings the instructions from \url{http://www.sthda.com/english/wiki/word-cloud-generator-in-r-one-killer-function-to-do-everything-you-need} .
}
\value{
A data frame with words and frequencies.
}
\examples{
# There is no example
NULL

}
