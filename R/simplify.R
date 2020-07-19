

# == title
# Simplify Gene Ontology (GO) enrichment results
#
# == param
# -mat A GO similarity matrix.
# -method Method for clustering the matrix. See `cluster_terms`.
# -control A list of parameters for controlling the clustering method, passed to `cluster_terms`.
# -plot Whether to make the heatmap.
# -term The full name or the description of the corresponding GO IDs. 
# -column_title Column title for the heatmap.
# -verbose Whether to print messages.
# -ht_list A list of additional heatmaps added to the left of the similarity heatmap.
# -... Arguments passed to `ht_clusters`.
#
# == details
# This is basically a wrapper function that it first runs `cluster_terms` to cluster
# GO terms and then runs `ht_clusters` to visualize the clustering.
#
# The arguments in `simplifyGO` passed to `ht_clusters` are:
#
# -``draw_word_cloud`` Whether to draw the word clouds.
# -``min_term`` Minimal number of GO terms in a cluster. All the clusters
#     with size less than ``min_term`` are all merged into one single cluster in the heatmap.
# -``order_by_size`` Whether to reorder GO clusters by their sizes. The cluster
#      that is merged from small clusters (size < ``min_term``) is always put to the bottom of the heatmap.
# -``exclude_words`` Words that are excluded in the word cloud.
# -``max_words`` Maximal number of words visualized in the word cloud.
# -``word_cloud_grob_param`` A list of graphic parameters passed to `word_cloud_grob`.
# -``fontsize_range`` The range of the font size. The value should be a numeric vector with length two.
#       The minimal font size is mapped to word frequency value of 1 and the maximal font size is mapped
#       to the maximal word frequency. The font size interlopation is linear.
#
# == value
# A data frame with three columns: GO IDs, GO term names and cluster labels.
#
# == example
# \donttest{
# set.seed(123)
# go_id = random_GO(500)
# mat = GO_similarity(go_id)
# df = simplifyGO(mat, word_cloud_grob_param = list(max_width = 80))
# head(df)
# }
simplifyGO = function(mat, method = "binary_cut", control = list(), 
	plot = TRUE, term = NULL, verbose = TRUE, 
	column_title = qq("@{nrow(mat)} GO terms clustered by '@{method}'"),
	ht_list = NULL, ...) {
	
	cl = do.call(cluster_terms, list(mat = mat, method = method, verbose = verbose, control = control))
	go_id = rownames(mat)

	if(!all(grepl("^GO:\\d+$", go_id))) {
		stop_wrap("Please ensure GO IDs are the row names of the similarity matrix and should be matched to '^GO:\\d+$'.")
	}

	if(is.null(term)) {
		suppressMessages(term <- select(GO.db::GO.db, keys = go_id, columns = "TERM")$TERM)
	}

	if(plot) ht_clusters(mat, cl, term = term, column_title = column_title, ht_list = ht_list, ...)

	return(invisible(data.frame(id = go_id, term = term, cluster = cl, stringsAsFactors = FALSE)))
}

# == title
# Simplify functional enrichment results
#
# == param
# -mat A similarity matrix.
# -method Method for clustering the matrix. See `cluster_terms`.
# -control A list of parameters for controlling the clustering method, passed to `cluster_terms`.
# -plot Whether to make the heatmap.
# -term The full name or the description of the corresponding terms. 
# -column_title Column title for the heatmap.
# -verbose Whether to print messages.
# -ht_list A list of additional heatmaps added to the left of the similarity heatmap.
# -... Arguments passed to `ht_clusters`.
#
# == details
# The usage is the same as `simplifyGO`, except you need to manually provide the term names by ``term`` argument.
#
simplifyEnrichment = function(mat, method = "binary_cut", control = list(), 
	plot = TRUE, term = NULL, verbose = TRUE, 
	column_title = qq("@{nrow(mat)} terms clustered by '@{method}'"),
	ht_list = NULL, ...) {
	
	cl = do.call(cluster_terms, list(mat = mat, method = method, verbose = verbose, control = control))
	term_id = rownames(mat)
	
	if(plot) ht_clusters(mat, cl, term = term, column_title = column_title, ht_list = ht_list, ...)

	if(is.null(term)) {
		return(data.frame(id = term_id, cluster = cl, stringsAsFactors = FALSE))
	} else {
		return(data.frame(id = term_id, term = term, cluster = cl, stringsAsFactors = FALSE))
	}
}
