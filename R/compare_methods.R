
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
# The function compares following default clustering methods by default:
#
# -``kmeans`` see `cluster_by_kmeans`.
# -``dynamicTreeCut`` see `cluster_by_dynamicTreeCut`.
# -``mclust`` see `cluster_by_mclust`. By default it is not included.
# -``apcluster`` see `cluster_by_apcluster`.
# -``hdbscan`` see `cluster_by_hdbscan`.
# -``fast_greedy`` see `cluster_by_igraph`.
# -``leading_eigen`` see `cluster_by_igraph`.
# -``louvain`` see `cluster_by_igraph`.
# -``walktrap`` see `cluster_by_igraph`.
# -``MCL`` see `cluster_by_MCL`.
# -``binary_cut`` see `binary_cut`.
#
# Also the user-defined methods in `all_clustering_methods` are also compared.
#
# == value
# A list of cluster label vectors for different clustering methods.
#
# == examples
# \dontrun{
# mat = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds",
#     package = "simplifyEnrichment"))
# clt = cmp_make_clusters(mat)
# }
cmp_make_clusters = function(mat, method = setdiff(all_clustering_methods(), "mclust"),
	verbose = TRUE) {

	clt = list()
	if("all" %in% method) {
		method = all_clustering_methods()
	}

	clt = lapply(method, function(me) {
		oe = try(cl <- cluster_terms(mat, me, verbose = verbose))
		if(inherits(oe, "try-error")) {
			rep(NA, nrow(mat))
		} else {
			cl
		}
	})
	names(clt) = method
	clt = as.data.frame(clt)

	clt
}

# == title
# Make plots for comparing clustering methods
#
# == param
# -mat A similarity matrix.
# -clt A list of clusterings from `cmp_make_clusters`.
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
#
# == examples
# \dontrun{
# mat = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds",
#     package = "simplifyEnrichment"))
# clt = cmp_make_clusters(mat)
# cmp_make_plot(mat, clt)
# cmp_make_plot(mat, clt, plot_type = "heatmap")
# }
cmp_make_plot = function(mat, clt, plot_type = c("mixed", "heatmap"), nrow = 3) {

	clt = as.data.frame(clt)
	methods = names(clt)
	
	plot_type = match.arg(plot_type)[1]

	check_pkg("cowplot", bioc = FALSE)
		
	if(plot_type == "mixed") {

		l = sapply(clt, function(x) any(is.na(x)))
		clt2 = clt[!l]
		methods = methods[!l]
		
		if("binary_cut" %in% names(clt)) {
			ref_class = clt2[, "binary_cut"]
		} else {
			ref_class = clt2[, which.min(vapply(clt2, function(x) length(unique(x)), 0))]
		}
		clt2 = lapply(clt2, function(x) as.character(relabel_class(x, ref_class, return_map = FALSE)))
		clt2 = as.data.frame(clt2)

		ht1 = Heatmap(mat, col = colorRamp2(c(0, quantile(mat, 0.975)), c("white", "red")),
			name = "Similarity",
			show_row_names = FALSE, show_column_names = FALSE, 
			# cluster_rows = dend, cluster_columns = dend,
			show_row_dend = FALSE, show_column_dend = FALSE) +
			Heatmap(as.matrix(clt2), show_heatmap_legend = FALSE, 
				width = unit(5, "mm")*ncol(clt2), column_names_rot = 45)
		p0 = grid.grabExpr(draw(ht1))

		stats = cmp_calc_stats(mat, clt2)
		stats$method = factor(rownames(stats), levels = rownames(stats))

		check_pkg("ggplot2", bioc = FALSE)

		suppressWarnings(
			p1 <- ggplot2::ggplot(stats, ggplot2::aes(x = stats$method, y = stats$diff_s)) +
			ggplot2::geom_bar(stat = "identity") + ggplot2::ylab("Difference score") +
			ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank())
		)

		df1 = stats[, c("method", "n_all")]; colnames(df1) = c("method", "value")
		df2 = stats[, c("method", "n_large")]; colnames(df2) = c("method", "value")
		df1$type = "All sizes"
		df2$type = "Size \u2265 5"
		df = rbind(df1, df2)
		suppressWarnings(
			p2 <- ggplot2::ggplot(df, ggplot2::aes(x = df$method, y = df$value, col = df$type, fill = df$type)) +
			ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) + ggplot2::ylab("Number of clusters") + ggplot2::labs(col = "Type", fill = "Type") +
			ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank())
		)

		suppressWarnings(
			p3 <- ggplot2::ggplot(stats, ggplot2::aes(x = stats$method, y = stats$block_mean)) +
			ggplot2::geom_bar(stat = "identity") + ggplot2::ylab("Block mean") +
			ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
		)

		cm = cmp_calc_concordance(clt2)
		p4 = grid.grabExpr(draw(Heatmap(cm, col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), name = "Concordance", column_names_rot = 45)))

		suppressWarnings(print(cowplot::plot_grid(
			cowplot::plot_grid(p0, p4, ncol = 1), 
			cowplot::plot_grid(p1, p2, p3, nrow = 3, align = "v", axis = "lr", rel_heights = c(1, 1, 1.5)),
			nrow = 1
		)))

	} else if(tolower(plot_type) == "heatmap") {
		pl = list()
		lgd = NULL

		for(i in seq_along(methods)) {
			if(any(is.na(clt[[i]]))) {
				pl[[i]] = textGrob(qq("@{methods[i]}\nan error occured."))
			} else {
				pl[[i]] = grid.grabExpr(ht <- ht_clusters(mat, clt[[i]], draw_word_cloud = FALSE, 
					column_title = qq("@{nrow(mat)} terms clustered by '@{methods[i]}'"),
					show_heatmap_legend = FALSE))
				lgd1 = color_mapping_legend(ht@ht_list[[1]]@matrix_color_mapping, plot = FALSE,
					legend_direction = "horizontal", title_position = "lefttop")
				lgd2 = Legend(labels = "Small clusters (size < 5)", legend_gp = gpar(fill = "darkgreen"))
				lgd = packLegend(lgd1, lgd2)
			}
		}

		n_all = sapply(clt, function(x) {
			if(any(is.na(x))) {
				NA
			} else {
				length(unique(x))
			}
		})
		n_big = sapply(clt, function(x) {
			if(any(is.na(x))) {
				NA
			} else {
				tb = table(x)
				sum(tb >= 5)
			}
		})

		tb = data.frame(Method = methods, "All clusters" = n_all, "Large clusters (size \u2265 5)" = n_big, check.names = FALSE)

		check_pkg("gridExtra", bioc = FALSE)

		pl[[length(pl) + 1]] = gridExtra::tableGrob(tb, rows = NULL)

		np = length(pl)

		ncol = ceiling(np/nrow)
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow = nrow, ncol = ncol)))
		for(i in 1:np) {
			ir = ceiling(i/ncol)
			ic = i %% ncol; if(ic == 0) ic = ncol
			pushViewport(viewport(layout.pos.row = ir, layout.pos.col = ic))
			if(i < np) {
				grid.draw(pl[[i]])
			} else {
				pushViewport(viewport(x = unit(0, "npc") + unit(2, "mm"), y = unit(1, "npc") - unit(1, "cm"),
					width = sum(pl[[i]]$widths), height = sum(pl[[i]]$heights),
					just = c("left", "top")))
				grid.draw(pl[[i]])
				grid.text("Number of clusters", y = unit(1, "npc") + unit(2.5, "mm"), just = "bottom", gp = gpar(fontsize = 14))
				popViewport()

				if(!is.null(lgd)) {
					draw(lgd, x = unit(0, "npc") + unit(2, "mm"), y = unit(1, "npc") - sum(pl[[i]]$heights) - unit(1.5, "cm"), just = c("left", "top"))
				}
			}
			popViewport()
		}
	}
}

cmp_calc_concordance = function(clt) {

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

cmp_calc_stats = function(mat, clt) {
	x = data.frame("diff_s" = vapply(clt, function(x) difference_score(mat, x), 0),
		 "n_all" = vapply(clt, function(x) length(table(x)), 0),
	     "n_large" = vapply(clt, function(x) {tb = table(x); sum(tb >= 5)}, 0),
		 "block_mean" = vapply(clt, function(x) block_mean(mat, x), 0))
	return(x)
}

# == title
# Compare clustering methods
#
# == param
# -mat The similarity matrix.
# -method Which methods to compare. All available methods are in `all_clustering_methods`.
#         A value of ``all`` takes all available methods. By default ``mclust`` is excluded because its long runtime.
# -plot_type See explanation in `cmp_make_plot`.
# -nrow Number of rows of the layout when ``plot_type`` is set to ``heatmap``.
# -verbose Whether to print messages.
#
# == details
# The function compares following clustering methods by default:
#
# -``kmeans`` see `cluster_by_kmeans`.
# -``dynamicTreeCut`` see `cluster_by_dynamicTreeCut`.
# -``mclust`` see `cluster_by_mclust`. By default it is not included.
# -``apcluster`` see `cluster_by_apcluster`.
# -``hdbscan`` see `cluster_by_hdbscan`.
# -``fast_greedy`` see `cluster_by_igraph`.
# -``leading_eigen`` see `cluster_by_igraph`.
# -``louvain`` see `cluster_by_igraph`.
# -``walktrap`` see `cluster_by_igraph`.
# -``MCL`` see `cluster_by_MCL`.
# -``binary_cut`` see `binary_cut`.
#
# This functon is basically a wrapper function. It calls the following two functions:
#
# - `cmp_make_clusters`: applies clustering with different methods.
# - `cmp_make_plot`: makes the plots.
#
# == value
# No value is returned.
#
# == example
# \dontrun{
# mat = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds",
#     package = "simplifyEnrichment"))
# compare_clustering_methods(mat)
# compare_clustering_methods(mat, plot_type = "heatmap")
# }
compare_clustering_methods = function(mat, method = setdiff(all_clustering_methods(), "mclust"),
	plot_type = c("mixed", "heatmap"), nrow = 3, verbose = TRUE) {

	clt = cmp_make_clusters(mat, method, verbose = verbose)

	plot_type = match.arg(plot_type)[1]
	cmp_make_plot(mat, clt, plot_type = plot_type, nrow = nrow)
}

