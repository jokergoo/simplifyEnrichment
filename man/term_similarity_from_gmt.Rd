\name{term_similarity_from_gmt}
\alias{term_similarity_from_gmt}
\title{
Similarity between terms from a gmt file
}
\description{
Similarity between terms from a gmt file
}
\usage{
term_similarity_from_gmt(term_id, gmt, extract_term_id = NULL, ...)
}
\arguments{

  \item{term_id}{A vector of terms.}
  \item{gmt}{The path of the gmt file.}
  \item{extract_term_id}{If the term ID in contained in the first column only as a substring, setting a function to extract this substring.}

}
\value{
A symmetric matrix
}
\examples{
# There is no example
NULL

}
