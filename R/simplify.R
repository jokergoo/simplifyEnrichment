

# == title
# Simplify GO enrichment results
#
# == param
# -mat A GO similarity matrix.
# -method Method for clustering the matrix. See `cluster_GO`.
# -control A list of parameters passed to `cluster_GO`.
# -plot Whether to make the heatmap.
# -verbose Whether to print messages.
# -... Arguments passed to `ht_GO_clusters`.
#
# == details
# This is basically a wrapper function that it first runs `cluster_GO` to cluster
# GO terms and then runs `ht_GO_clusters` to visualize the clustering.
#
# The arguments passed to `ht_GO_clusters` are:
#
# -``draw_word_cloud`` Whether to draw the word clouds.
# -``min_term`` Minimal number of GO terms in a cluster. All the clusters
#     with size less than ``min_term`` are all merged into one single cluster in the heatmap.
# -``order_by_size`` Whether to reorder GO clusters by their sizes. The cluster
#      that is merged from small clusters (size < 5) is always put to the bottom of the heatmap.
# -``exclude_words`` Words that are excluded in the word cloud.
# -``max_words`` Maximal number of words visualized in the word cloud.
# -``word_cloud_grob_param`` A list of parameters passed to `word_cloud_grob`.
# -``fontsize_range`` The range of the font size. The value should be a numeric vector with length two.
#       The minimal font size is mapped to word frequency value of 1 and the maximal font size is mapped
#       to the maximal word frequency. The font size interlopation is linear.
#
# == value
# A data frame with three columns: GO IDs, GO term names and cluster labels.
#
# == example
# \dontrun{
# set.seed(123)
# go_id = simplifyGO:::random_GO(500)
# mat = GO_similarity(go_id)
# df = simplifyGO(mat)
# }
simplifyGO = function(mat, method = "binary_cut", control = list(), 
	plot = TRUE, verbose = TRUE, ...) {
	
	cl = do.call(cluster_GO, c(list(mat = mat, method = method, verbose = verbose), control))

	if(plot) ht_GO_clusters(mat, cl, column_title = qq("@{nrow(mat)} GO terms are clustered by '@{method}'"), ...)

	go_id = rownames(mat)
	term = select(GO.db::GO.db, keys = go_id, columns = "TERM")$TERM
	
	return(data.frame(id = go_id, name = term, cluster = cl, stringsAsFactors = FALSE))
}

