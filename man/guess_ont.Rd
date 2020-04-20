\name{guess_ont}
\alias{guess_ont}
\title{
Guess the ontology of the input GO IDs
}
\description{
Guess the ontology of the input GO IDs
}
\usage{
guess_ont(go_id, db = 'org.Hs.eg.db')
}
\arguments{

  \item{go_id}{A vector of GO IDs.}
  \item{db}{Annotation database. It should be from \url{https://bioconductor.org/packages/3.10/BiocViews.html#___OrgDb}}

}
\details{
10 GO IDs are randomly sampled and tested.
}
\value{
A single character scalar of "BP", "CC" or "MF".

If there are more than one ontologies detected. It returns \code{NULL}.
}
\examples{
# There is no example
NULL

}
