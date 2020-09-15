\name{select_cutoff}
\alias{select_cutoff}
\title{
Select the cutoff for binary cut
}
\description{
Select the cutoff for binary cut
}
\usage{
select_cutoff(mat, cutoff = seq(0.6, 0.98, by = 0.01), verbose = TRUE, ...)
}
\arguments{

  \item{mat}{A similarity matrix.}
  \item{cutoff}{A list of cutoffs to test. Note the range of the cutoff values should be inside [0.5, 1].}
  \item{verbose}{Whether to print messages.}
  \item{...}{Pass to \code{\link{binary_cut}}.}

}
\details{
Binary cut is applied to each of the cutoff and the clustering results are evaluated by following metrics:

\itemize{
  \item difference score, calculated by \code{\link{difference_score}}.
  \item number of clusters.
  \item block mean, which is the mean similarity in the blocks in the diagonal of the heatmap.
}
}
\examples{
\donttest{
mat = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds",
    package = "simplifyEnrichment"))
select_cutoff(mat)
}
}
