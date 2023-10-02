\name{term_similarity}
\alias{term_similarity}
\title{
Similarity between terms based on the overlap of genes
}
\description{
Similarity between terms based on the overlap of genes
}
\usage{
term_similarity(gl, method = c("kappa", "jaccard", "dice", "overlap"), all = NULL, remove_negative = TRUE)
}
\arguments{

  \item{gl}{A list of genes that are in the terms.}
  \item{method}{The similarity measurement.}
  \item{all}{The universe set.}
  \item{remove_negative}{If the value is negative, reset to zero}

}
\details{
The definition of the four similarity measurements can be found at \url{https://jokergoo.github.io/simplifyEnrichment_supplementary/supplS01_coefficient_definition/supplS01_coefficient_definition.html} .
}
\value{
A symmetric matrix.
}
\examples{
# There is no example
NULL

}
