\name{random_GO}
\alias{random_GO}
\title{
Generate random GO IDs
}
\description{
Generate random GO IDs
}
\usage{
random_GO(n, ont = "BP", db = 'org.Hs.eg.db')
}
\arguments{

  \item{n}{Number of GO IDs.}
  \item{ont}{GO ontology. Value should be one of "BP", "CC" or "MF".}
  \item{db}{Annotation database. It should be from \url{https://bioconductor.org/packages/3.10/BiocViews.html#___OrgDb}}

}
\value{
A vector of GO IDs.
}
\examples{
# There is no example
NULL

}
