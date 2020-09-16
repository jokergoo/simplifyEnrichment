\name{partition_by_pam}
\alias{partition_by_pam}
\title{
Partition by PAM
}
\description{
Partition by PAM
}
\usage{
partition_by_pam(mat)
}
\arguments{

  \item{mat}{The similarity matrix.}

}
\details{
The clustering is performed by \code{\link[cluster]{pam}} with setting \code{pamonce} argument to 5.

This function is used to set to the \code{partition_fun} argument in \code{\link{binary_cut}}.
}
\examples{
# There is no example
NULL

}
