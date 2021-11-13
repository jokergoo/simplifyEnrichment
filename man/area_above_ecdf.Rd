\name{area_above_ecdf}
\alias{area_above_ecdf}
\title{
Area above the eCDF curve
}
\description{
Area above the eCDF curve
}
\usage{
area_above_ecdf(x)
}
\arguments{

  \item{x}{A vector of similarity values.}

}
\details{
Denote F(x) as the eCDF (empirical Cumulative Distribution Function) of the similarity vector \code{x}, this function calculates
the area above the eCDF curve, which is 1 - \\int_0^1 F(x)dx.
}
\value{
A numeric value.
}
\examples{
# There is no example
NULL

}
