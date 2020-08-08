\name{term_similarity}
\alias{term_similarity}
\title{
Similarity between terms based on the overlap of genes
}
\description{
Similarity between terms based on the overlap of genes
}
\usage{
term_similarity(gl, method = c("kappa", "jaccard", "dice", "overlap"))
}
\arguments{

  \item{gl}{A list of genes that are in the terms.}
  \item{method}{The similarity measurement.}

}
\details{
The definition of the four similarity measurements can be found at \url{https://simplifyenrichment.github.io/supplementary/suppl1_coefficient_definition/suppl1_coefficient_definition.html} .
}
\value{
A symmetric matrix.
}
\examples{
# There is no example
NULL

}
