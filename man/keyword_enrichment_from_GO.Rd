\name{keyword_enrichment_from_GO}
\alias{keyword_enrichment_from_GO}
\title{
Keyword enrichment for GO terms
}
\description{
Keyword enrichment for GO terms
}
\usage{
keyword_enrichment_from_GO(go_id, min_bg = 5, min_term = 2)
}
\arguments{

  \item{go_id}{A vector of GO IDs.}
  \item{min_bg}{Minimal number of GO terms (in the background, i.e. all GO temrs in the GO database) that contain a specific keyword.}
  \item{min_term}{Minimal number of GO terms (GO terms in \code{go_id}) that contain a specific keyword.}

}
\details{
The enrichment is applied by Fisher's exact test. For a keyword, there is the following 2x2 contigency table:

  \preformatted{
                      | contains the keyword | does not contain the keyword
    In the GO set     |          s11         |          s12
    Not in the GO set |          s21         |          s22  }

where s11, s12, s21 and s22 are number of GO terms in each category.
}
\value{
A data frame with keyword enrichment results.
}
\examples{
\dontrun{
go_id = random_GO(100)
keyword_enrichment_from_GO(go_id)
}
}
