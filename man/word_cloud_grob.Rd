\name{word_cloud_grob}
\alias{word_cloud_grob}
\title{
A simple grob for the word cloud
}
\description{
A simple grob for the word cloud
}
\usage{
word_cloud_grob(text, fontsize,
    line_space = unit(4, "pt"), word_space = unit(4, "pt"), max_width = unit(80, "mm"),
    col = function(fs) circlize::rand_color(length(fs), luminosity = "dark"),
    test = FALSE)
}
\arguments{

  \item{text}{A vector of words.}
  \item{fontsize}{The corresponding font size. With the frequency of the words known, \code{\link{scale_fontsize}} can be used to linearly interpolate frequencies to font sizes.}
  \item{line_space}{Space between lines. The value can be a \code{\link[grid]{unit}} object or a numeric scalar which is measured in mm.}
  \item{word_space}{Space between words. The value can be a \code{\link[grid]{unit}} object or a numeric scalar which is measured in mm.}
  \item{max_width}{The maximal width of the viewport to put the word cloud. The value can be a \code{\link[grid]{unit}} object or a numeric scalar which is measured in mm. Note this might be larger than the final width of the returned grob object.}
  \item{col}{Colors for the words. The value can be a vector, in numeric or character, which should have the same length as \code{text}. Or it is a self-defined function that takes the font size vector as  the only argument. The function should return a color vector. See Examples.}
  \item{test}{Internally used. It basically adds borders to the words and the viewport.}

}
\value{
A \code{\link[grid:grid.grob]{grob}} object. The width and height of the grob can be get by \code{\link[grid]{grobWidth}} and \code{\link[grid:grobWidth]{grobHeight}}.
}
\examples{
# very old R versions do not have strrep() function
if(!exists("strrep")) {
    strrep = function(x, i) paste(rep(x, i), collapse = "")
}
words = sapply(1:30, function(x) strrep(sample(letters, 1), sample(3:10, 1)))
require(grid)
gb = word_cloud_grob(words, fontsize = runif(30, min = 5, max = 30), 
    max_width = 100)
grid.newpage(); grid.draw(gb)

# color as a single scalar
gb = word_cloud_grob(words, fontsize = runif(30, min = 5, max = 30), 
    max_width = 100, col = 1)
grid.newpage(); grid.draw(gb)

# color as a vector
gb = word_cloud_grob(words, fontsize = runif(30, min = 5, max = 30), 
    max_width = 100, col = 1:30)
grid.newpage(); grid.draw(gb)

# color as a function
require(circlize)
col_fun = colorRamp2(c(5, 17, 30), c("blue", "black", "red"))
gb = word_cloud_grob(words, fontsize = runif(30, min = 5, max = 30), 
    max_width = 100, col = function(fs) col_fun(fs))
grid.newpage(); grid.draw(gb)
}
