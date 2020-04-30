\name{scale_fontsize}
\alias{scale_fontsize}
\title{
Scale font size
}
\description{
Scale font size
}
\usage{
scale_fontsize(x, rg = c(1, 30), fs = c(4, 16))
}
\arguments{

  \item{x}{A numeric vector.}
  \item{rg}{The range.}
  \item{fs}{Range of th font size.}

}
\section{Detaisl}{
It is a linear interpolation.}
\value{
A numeric vector.
}
\examples{
scale_fontsize(runif(10, min = 1, max = 20))
}
