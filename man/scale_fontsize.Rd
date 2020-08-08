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
  \item{fs}{Range of the font size.}

}
\section{Detaisl}{
It is a linear interpolation.}
\value{
A numeric vector.
}
\examples{
x = runif(10, min = 1, max = 20)
# scale x to fontsize 4 to 16.
scale_fontsize(x)
}
