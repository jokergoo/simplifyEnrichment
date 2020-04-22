

compare_methods_make_clusters = function(mat, method = "all") {

	clt = list()
	if(method == "all") {
		method = c("kmeans", 
			       "dynamicTreeCut", 
			       # "mclust", 
			       "apcluster", 
			       "fast_greedy", 
			       "leading_eigen", 
			       "louvain", 
			       "walktrap", 
			       "binary_cut")
	}

	clt = lapply(method, function(me) as.character(cluster_GO(mat, me)))
	names(clt) = method
	clt = as.data.frame(clt)

	clt
}

compare_methods_make_plot = function(mat, clt) {

	ht = Heatmap(mat, col = colorRamp2(c(0, 1), c("white", "red")),
		name = "Similarity", column_title = "GO Similarity",
		show_row_names = FALSE, show_column_names = FALSE, 
		# cluster_rows = dend, cluster_columns = dend,
		show_row_dend = FALSE, show_column_dend = FALSE,
		right_annotation = rowAnnotation(df = clt, show_legend = FALSE))
	p1 = grid.grabExpr(draw(ht))

	x = sapply(clt, function(x) difference_score(mat, x))

	if(!requireNamespace("ggplot2", quietly = TRUE)) {
		stop_wrap("Package ggplot2 should be installed.")
	}
	p2 = ggplot2::ggplot(NULL, ggplot2::aes(x = names(x), y = x)) + ggplot2::geom_bar(stat = "identity") +
		ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

	cm = compare_methods_calc_concordance(clt)
	ht2 = Heatmap(cm, column_title = "concordance")
	p3 = grid.grabExpr(draw(ht2))

	if(!requireNamespace("cowplot", quietly = TRUE)) {
		stop_wrap("Package cowplot should be installed.")
	}
	print(cowplot::plot_grid(p1, cowplot::plot_grid(p2, p3, nrow = 1), nrow = 2))
}

# == title
# Difference score
#
# == param
# -mat The similarity matrix.
# -cl Cluster labels.
#
# == details
# This function measures the different between the similarity values for the GO terms
# that belong to the same clusters and in different clusters. The difference score
# is the Kolmogorov-Smirnov statistic between the two distributions.
#
difference_score = function(mat, cl) {
	n = nrow(mat)
	l_block = matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat))

	for(le in unique(cl)) {
		l = cl == le
		l_block[l, l] = TRUE
	}
	l_block2 = l_block
	l_block2[upper.tri(mat)] = FALSE
	x1 = mat[l_block2]

	l_block2 = !l_block
	l_block2[upper.tri(mat)] = FALSE
	x2 = mat[l_block2]

	ecdf1 = ecdf(x1)
	ecdf2 = ecdf(x2)

	p = seq(0, 1, length = 1000)
	max(abs(ecdf2(p) - ecdf1(p)))
}

block_mean = function(mat, cl) {
	n = nrow(mat)
	l_block = matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat))

	for(le in unique(cl)) {
		l = cl == le
		l_block[l, l] = TRUE
	}
	l_block2 = l_block
	l_block2[upper.tri(mat)] = FALSE
	x1 = mat[l_block2]
	mean(x1)
}

other_mean = function(mat, cl) {
	n = nrow(mat)
	l_block = matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat))

	for(le in unique(cl)) {
		l = cl == le
		l_block[l, l] = TRUE
	}

	l_block2 = !l_block
	l_block2[upper.tri(mat)] = FALSE
	x2 = mat[l_block2]
	mean(x2)
}

compare_methods_calc_concordance = function(clt) {

	concordance = function(cl1, cl2) {
		cl1 = as.vector(cl1)
		cl2 = as.vector(cl2)
		nle1 = length(unique(cl1))
		nle2 = length(unique(cl2))

		if(nle1 > nle2) {
			cl2 = relabel_class(cl2, cl1, return_map = FALSE)
		} else {
			cl1 = relabel_class(cl1, cl2, return_map = FALSE)
		}

		sum(cl1 == cl2)/length(cl1)
	}

	n = length(clt)
	mm = matrix(1, nrow = n, ncol = n)
	rownames(mm) = names(clt)
	colnames(mm) = names(clt)
	for(i in 1:(n-1)) {
		for(j in (i+1):n) {
			mm[i, j] = concordance(clt[[i]], clt[[j]])
			mm[j, i] = mm[i, j]
		}
	}

	mm
}

# == title
# Compare clustering methods
#
# == param
# -mat The GO similarity matrix.
#
# == details
# The function compares following clustering methods:
#
# -``binary_cut`` see `binary_cut`.
# -``kmeans`` see `cluster_by_kmeans`.
# -``dynamicTreeCut`` see `cluster_by_dynamicTreeCut`.
# -``mclust`` see `cluster_by_mclust`.
# -``apcluster`` see `cluster_by_apcluster`.
# -``fast_greedy`` see `cluster_by_igraph`.
# -``leading_eigen`` see `cluster_by_igraph`.
# -``louvain`` see `cluster_by_igraph`.
# -``walktrap`` see `cluster_by_igraph`.
#
# The function produces a plot with four panels.
#
compare_methods = function(mat) {
	clt = compare_methods_make_clusters(mat, "all")
	compare_methods_make_plot(mat, clt)
}

