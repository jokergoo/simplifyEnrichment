\name{term_similarity}
\alias{term_similarity}
\title{
Similarity between terms based on the overlap of genes
}
\description{
Similarity between terms based on the overlap of genes
}
\usage{
term_similarity(gl, method = c("kappa", "jaccard"))
}
\arguments{

  \item{gl}{A list of genes that are in the terms.}
  \item{method}{Kappa coefficient or Jaccard coefficient.}

}
\details{
For two sets of genes, the Jaccard coefficient is calculated as:

  \preformatted{
    length(intersect(set1, set2))/length(union(set1, set2))  }

For the calculation of Kappa coefficient, see \url{https://en.wikipedia.org/wiki/Cohen\%27s_kappa} .
}
\value{
A symmetric matrix.
}
\examples{
# There is no example
NULL

}
