\name{count_word}
\alias{count_word}
\title{
Calculate word frequency
}
\description{
Calculate word frequency
}
\usage{
count_word(id, term = NULL, exclude_words = NULL)
}
\arguments{

  \item{id}{A vector of term IDs.}
  \item{term}{The corresponding names or description of terms.}
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
