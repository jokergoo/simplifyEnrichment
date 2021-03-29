\name{simplifyGOFromMultipleLists}
\alias{simplifyGOFromMultipleLists}
\title{
Perform simplifyGO analysis with multiple lists of GO IDs
}
\description{
Perform simplifyGO analysis with multiple lists of GO IDs
}
\usage{
simplifyGOFromMultipleLists(lt, go_id_column = NULL, padj_column = NULL, padj_cutoff = 1e-2,
    filter = function(x) any(x < padj_cutoff), default = 1,
    heatmap_param = list(NULL),
    method = "binary_cut", control = list(partial = TRUE),
    min_term = NULL, verbose = TRUE, column_title = NULL, ...)
}
\arguments{

  \item{lt}{A data frame, a list of numeric vectors (e.g. adjusted p-values) where each numeric vector has GO IDs as names, or a list of GO IDs.}
  \item{go_id_column}{Column index of GO ID if \code{lt} contains a list of data frames.}
  \item{padj_column}{Column index of adjusted p-values if \code{lt} contains a list of data frames.}
  \item{padj_cutoff}{Cut off for adjusted p-values}
  \item{filter}{A self-defined function for filtering GO IDs. By default it requires GO IDs should be significant in at least one list.}
  \item{default}{The default value for the adjusted p-values. See Details.}
  \item{heatmap_param}{Parameters for controlling the heatmap, see Details.}
  \item{method}{Pass to \code{\link{simplifyGO}}.}
  \item{control}{Pass to \code{\link{simplifyGO}}.}
  \item{min_term}{Pass to \code{\link{simplifyGO}}.}
  \item{verbose}{Pass to \code{\link{simplifyGO}}.}
  \item{column_title}{Pass to \code{\link{simplifyGO}}.}
  \item{...}{Pass to \code{\link{simplifyGO}}.}

}
\details{
The input data can have three types of formats:

\itemize{
  \item A list of numeric vectors of adjusted p-values where each vector has the GO IDs as names.
  \item A data frame. The column of the GO IDs can be specified with \code{go_id_column} argument and the column of the adjusted p-values can be specified with \code{padj_column} argument. If these columns are not specified, they are automatically identified. The GO ID column is found by checking whether a column contains all GO IDs. The adjusted p-value column is found by comparing the column names of the  data frame to see whether it might be a column for adjusted p-values. These two columns are used to construct a numeric vector with GO IDs as names.
  \item A list of character vectors of GO IDs. In this case, each character vector is changed to a numeric vector where all values take 1 and the original GO IDs are used as names of the vector.
}

Now let's assume there are \code{n} GO lists, we first construct a global matrix where columns correspond to the \code{n} GO lists and rows correspond
to the "union" of all GO IDs in the lists. The value for the ith GO ID and in the jth list are taken from the corresponding numeric vector
in \code{lt}. If the jth vector in \code{lt} does not contain the ith GO ID, the value defined by \code{default} argument is taken there (e.g. in most cases the numeric
values are adjusted p-values, \code{default} is set to 1). Let's call this matrix as \code{M0}.

Next step is to filter \code{M0} so that we only take a subset of GO IDs of interest. We define a proper function via argument \code{filter} to remove
GO IDs that are not important for the analysis. Functions for \code{filter} is applied to every row in \code{M0} and \code{filter} function needs
to return a logical value to decide whether to remove the current GO ID. For example, if the values in \code{lt} are adjusted p-values, the \code{filter} function
can be set as \code{function(x) any(x < padj_cutoff)} so that the GO ID is kept as long as it is signfiicant in at least one list. After the filter, let's call
the filtered matrix \code{M1}.

GO IDs in \code{M1} (row names of \code{M1}) are used for clustering. A heatmap of \code{M1} is attached to the left of the GO similarity heatmap so that
the group-specific (or list-specific) patterns can be easily observed and to corresponded to GO functions.

Argument \code{heatmap_param} controls several parameters for heatmap \code{M1}:

\itemize{
  \item \code{transform}: A self-defined function to transform the data for heatmap visualization. The most typical case is to transform adjusted p-values by \code{-log10(x)}.
  \item \code{breaks}: break values for color interpolation.
  \item \code{col}: The corresponding values for \code{breaks}.
  \item \code{labels}: The corresponding labels.
  \item \code{name}: Legend title.
}
}
\examples{
\donttest{
# perform functional enrichment on the signatures genes from cola anlaysis 
require(cola)
data(golub_cola) 
res = golub_cola["ATC:skmeans"]
require(hu6800.db)
x = hu6800ENTREZID
mapped_probes = mappedkeys(x)
id_mapping = unlist(as.list(x[mapped_probes]))
lt = functional_enrichment(res, k = 3, id_mapping = id_mapping) # you can check the value of `lt`

# a list of data frames
simplifyGOFromMultipleLists(lt, padj_cutoff = 0.001)

# a list of numeric values
lt2 = lapply(lt, function(x) structure(x$p.adjust, names = x$ID))
simplifyGOFromMultipleLists(lt2, padj_cutoff = 0.001)

# a list of GO IDS
lt3 = lapply(lt, function(x) x$ID[x$p.adjust < 0.001])
simplifyGOFromMultipleLists(lt3)
}
}
