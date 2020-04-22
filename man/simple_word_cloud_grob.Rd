\name{simple_word_cloud_grob}
\alias{simple_word_cloud_grob}
\title{
A simple grob for word cloud
}
\description{
A simple grob for word cloud
}
\usage{
simple_word_cloud_grob(text, fontsize, max_width = 40)
}
\arguments{

  \item{text}{A vector of words.}
  \item{fontsize}{The corresponding font size.}
  \item{max_width}{The maximal width of the viewport to put the word cloud.  The value should be numeric. It is measured in mm.}

}
\value{
A \code{\link[grid:grid.grob]{grob}} object. The width and height of the grob can be get by \code{\link[grid]{grobWidth}} and \code{\link[grid:grobWidth]{grobHeight}}.
}
\examples{
if(!exists("strrep")) {
	strrep = function(x, i) paste(rep(x, i), collapse = "")
}
words = sapply(1:30, function(x) strrep(sample(letters, 1), sample(3:10, 1)))
require(grid)
gb = simple_word_cloud_grob(words, fontsize = runif(30, min = 5, max = 30), max_width = 80)
grid.newpage()
w = grobWidth(gb)
h = grobHeight(gb)
pushViewport(viewport(width = w, height = grobHeight(gb)))
grid.draw(gb)
grid.rect()
popViewport()
}
