
# == title
# Apply various clustering methods
#
# == param
# -mat The similarity matrix.
# -method Which methods to compare. All available methods are in `all_clustering_methods`.
#         A value of ``all`` takes all available methods. By default ``mclust`` is excluded because its long runtime.
# -verbose Whether to print messages.
#
# == details
# The function compares following default clustering methods:
#
# -``binary_cut`` see `binary_cut`.
# -``kmeans`` see `cluster_by_kmeans`.
# -``dynamicTreeCut`` see `cluster_by_dynamicTreeCut`.
# -``mclust`` see `cluster_by_mclust`. By default it is not included.
# -``apcluster`` see `cluster_by_apcluster`.
# -``fast_greedy`` see `cluster_by_igraph`.
# -``leading_eigen`` see `cluster_by_igraph`.
# -``louvain`` see `cluster_by_igraph`.
# -``walktrap`` see `cluster_by_igraph`.
#
# Also the user-defined methods in `all_clustering_methods` are also compared.
#
# == value
# A list of cluster label vectors for different clustering methods.
#
compare_methods_make_clusters = function(mat, method = setdiff(all_clustering_methods(), "mclust"),
	verbose = TRUE) {

	clt = list()
	if("all" %in% method) {
		method = all_clustering_methods()
	}

	clt = lapply(method, function(me) cluster_terms(mat, me, verbose = verbose))
	names(clt) = method
	clt = as.data.frame(clt)

	clt
}

# == title
# Make plots for comparing clustering methods
#
# == param
# -mat A similarity matrix.
# -clt A list of clusterings from `compare_methods_make_clusters`.
# -plot_type What type of plots to make. See Details.
# -nrow Number of rows of the layout when ``plot_type`` is set to ``heatmap``.
#
# == details
# If ``plot_type`` is the default value ``mixed``, a figure with three panels generated:
#
# - A heatmap of the similarity matrix with different classifications as row annotations.
# - A heatmap of the pair-wise concordance of the classifications of every two clustering methods.
# - Barplots of the difference scores for each method (calculated by `difference_score`), the number
#    of clusters (total clusters and the clusters with size >= 5) and the mean similarity of the terms 
#    that are in the same clusters.
#
# If ``plot_type`` is ``heatmap``. There are heatmaps for the similarity matrix under clusterings
# from different methods. The last panel is a table with the number of clusters under different
# clusterings.
#
# == value
# No value is returned.
compare_methods_make_plot = function(mat, clt, plot_type = c("mixed", "heatmap"), nrow = 2) {

	clt = lapply(clt, as.character)
	clt = as.data.frame(clt)
	methods = names(clt)
	
	plot_type = match.arg(plot_type)[1]

	if(!requireNamespace("cowplot", quietly = TRUE)) {
		stop_wrap("Package cowplot should be installed.")
	}
		
	if(plot_type == "mixed") {
		
		if("binary_cut" %in% names(clt)) {
			ref_class = clt[, "binary_cut"]
		} else {
			ref_class = clt[, which.min(sapply(clt, function(x) length(unique(x))))]
		}
		clt2 = lapply(clt, function(x) relabel_class(x, ref_class, return_map = FALSE))
		clt2 = as.data.frame(clt2)

		ht1 = Heatmap(mat, col = colorRamp2(c(0, 1), c("white", "red")),
			name = "Similarity",
			show_row_names = FALSE, show_column_names = FALSE, 
			# cluster_rows = dend, cluster_columns = dend,
			show_row_dend = FALSE, show_column_dend = FALSE) +
			Heatmap(as.matrix(clt2), show_heatmap_legend = FALSE, 
				width = unit(5, "mm")*ncol(clt2), column_names_rot = 45)
		p0 = grid.grabExpr(draw(ht1))

		stats = compare_methods_calc_stats(mat, clt)
		stats$method = factor(rownames(stats), levels = rownames(stats))

		if(!requireNamespace("ggplot2", quietly = TRUE)) {
			stop_wrap("Package ggplot2 should be installed.")
		}
		suppressWarnings(
			p1 <- ggplot2::ggplot(stats, ggplot2::aes(x = stats$method, y = stats$diff_s)) +
			ggplot2::geom_bar(stat = "identity") + ggplot2::ylab("Difference score") +
			ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank())
		)

		df1 = stats[, c("method", "n_all")]; colnames(df1) = c("method", "value")
		df2 = stats[, c("method", "n_large")]; colnames(df2) = c("method", "value")
		df1$type = "All sizes"
		df2$type = "size >= 5"
		df = rbind(df1, df2)
		suppressWarnings(
			p2 <- ggplot2::ggplot(df, ggplot2::aes(x = df$method, y = df$value, col = df$type, fill = df$type)) +
			ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) + ggplot2::ylab("Cluster number") + ggplot2::labs(col = "Type", fill = "Type") +
			ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank())
		)

		suppressWarnings(
			p3 <- ggplot2::ggplot(stats, ggplot2::aes(x = stats$method, y = stats$block_mean)) +
			ggplot2::geom_bar(stat = "identity") + ggplot2::ylab("Block mean") +
			ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
		)

		cm = compare_methods_calc_concordance(clt)
		p4 = grid.grabExpr(draw(Heatmap(cm, name = "Concordance", column_names_rot = 45)))

		suppressWarnings(cowplot::plot_grid(
			cowplot::plot_grid(p0, p4, ncol = 1), 
			cowplot::plot_grid(p1, p2, p3, nrow = 3, align = "v", axis = "lr", rel_heights = c(1, 1, 1.5)),
			nrow = 1
		))

} else {
		pl = list()
		for(i in seq_along(methods)) {
			pl[[i]] = grid.grabExpr(ht_clusters(mat, clt[[i]], draw_word_cloud = FALSE, column_title = qq("@{nrow(mat)} terms clustered by '@{methods[i]}'")))
		}

		stats = compare_methods_calc_stats(mat, clt)

		tb = data.frame(method = methods, "#clusters" = stats["n_all"], "#cluster(size >= 5)" = stats["n_large"], check.names = FALSE)

		if(!requireNamespace("gridExtra", quietly = TRUE)) {
			stop_wrap("Package gridExtra should be installed.")
		}
		pl[[length(pl) + 1]] = gridExtra::tableGrob(tb, rows = NULL)

		print(cowplot::plot_grid(plotlist = pl, nrow = nrow))
	}
}

# == title
# Difference score
#
# == param
# -mat The similarity matrix.
# -cl Cluster labels.
#
# == details
# This function measures the different between the similarity values for the terms
# that belong to the same clusters and in different clusters. The difference score
# is the Kolmogorov-Smirnov statistic between the two distributions.
#
# == value
# A numeric scalar.
#
# == examples
# mat = readRDS(system.file("extdata", "similarity_mat.rds", package = "simplifyEnrichment"))
# cl = binary_cut(mat)
# difference_score(mat, cl)
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

	clt = lapply(clt, as.character)
	clt = as.data.frame(clt)

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
	for(i in seq(1, (n-1))) {
		for(j in seq(i+1, n)) {
			mm[i, j] = concordance(clt[[i]], clt[[j]])
			mm[j, i] = mm[i, j]
		}
	}

	mm
}

compare_methods_calc_stats = function(mat, clt) {
	x = data.frame("diff_s" = sapply(clt, function(x) difference_score(mat, x)),
		 "n_all" = sapply(clt, function(x) length(table(x))),
	     "n_large" = sapply(clt, function(x) {tb = table(x); sum(tb >= 5)}),
		 "block_mean" = sapply(clt, function(x) block_mean(mat, x)))
	return(x)
	
}

# == title
# Compare clustering methods
#
# == param
# -mat The similarity matrix.
# -method Which methods to compare. All available methods are in `all_clustering_methods`.
#         A value of ``all`` takes all available methods. By default ``mclust`` is excluded because its long runtime.
# -plot_type See explanation in `compare_methods_make_plot`.
# -verbose Whether to print messages.
#
# == details
# The function compares following clustering methods:
#
# -``binary_cut`` see `binary_cut`.
# -``kmeans`` see `cluster_by_kmeans`.
# -``dynamicTreeCut`` see `cluster_by_dynamicTreeCut`.
# -``mclust`` see `cluster_by_mclust`. By default it is not included.
# -``apcluster`` see `cluster_by_apcluster`.
# -``fast_greedy`` see `cluster_by_igraph`.
# -``leading_eigen`` see `cluster_by_igraph`.
# -``louvain`` see `cluster_by_igraph`.
# -``walktrap`` see `cluster_by_igraph`.
#
# This functon is basically a wrapper function. It calls following two functions:
#
# - `compare_methods_make_clusters`: applies clustering by different methods.
# - `compare_methods_make_plot`: makes the plots.
#
# == value
# No value is returned.
#
# == example
# \dontrun{
# mat = readRDS(system.file("extdata", "similarity_mat.rds", package = "simplifyEnrichment"))
# compare_methods(mat)
# compare_methods(mat, plot_type = "heatmap")
# }
compare_methods = function(mat, method = setdiff(all_clustering_methods(), "mclust"),
	plot_type = c("mixed", "heatmap"), verbose = TRUE) {

	clt = compare_methods_make_clusters(mat, method, verbose = verbose)

	plot_type = match.arg(plot_type)[1]
	compare_methods_make_plot(mat, clt, plot_type = plot_type)
}

