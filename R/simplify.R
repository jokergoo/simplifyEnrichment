

# == title
# Simplify GO enrichment results
#
# == param
# -mat A GO similarity matrix
# -vavlue_fun
# -cutoff cutoff
# -plot Whether to make the plot
# -min_term
# -order_by_sizez
# -exclude_words
#
# == example
# \dontrun{
# set.seed(123)
# go_id = simplifyGO:::random_GO(500)
# mat = GO_similarity(go_id)
# lt = simplify(mat)
# }
simplify = function(mat, value_fun = median, cutoff = 0.8, plot = TRUE, min_term = 5, 
	order_by_size = TRUE, exclude_words = character(0)) {
	
	dend = cluster_mat(mat, value_fun = value_fun)
	cl = cut_dend(dend, cutoff)

	if(plot) plot_heatmap(mat, cl, min_term = min_term, order_by_size = order_by_size, exclude_words = exclude_words)

	go_id = rownames(mat)
	term = select(GO.db::GO.db, keys = go_id, columns = "TERM")$TERM
	
	return(split(data.frame(ID = go_id, Name = term, stringsAsFactors = FALSE), cl))
}

