

# == title
# Simplify GO enrichment results
#
# == param
# -mat A GO similarity matrix
# -cutoff cutoff
# -plot Whether to make the plot
#
# == example
# \dontrun{
# set.seed(123)
# go_id = simplifyGO:::random_GO(500)
# mat = GO_similarity(go_id)
# lt = simplify(mat)
# }
simplify = function(mat, cutoff = 0.8, plot = TRUE) {
	
	dend = cluster_mat(mat)
	cl = cut_dend(dend, cutoff)

	if(plot) plot_heatmap(mat, cl, min = 5)

	return(split(rownames(mat), cl))
}

