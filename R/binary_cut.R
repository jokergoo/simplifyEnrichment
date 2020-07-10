
dend_max_depth = function(dend) {
	max(unlist(dend_node_apply(dend, function(d, index) {
		length(index) + 1
	})))
}

# cluster the similarity matrix and assign scores to nodes
cluster_mat = function(mat, value_fun = median) {

	env = new.env()
	env$value_fun = value_fun
	.cluster_mat(mat, dist_mat = dist(mat), .env = env)

	dend = env$dend

	max_d = dend_max_depth(dend)
	dend = dendrapply(dend, function(d) {
		attr(d, "height") = max_d - attr(d, "height")
		d
	})

	# to correctly assign midpoint for each node
	dend2 = as.dendrogram(as.hclust(dend))

	dend = edit_node(dend, function(d, index) {
		if(is.null(index)) {  # top node
			if(!is.leaf(d)) {
				attr(d, "midpoint") = attr(dend2, "midpoint")
			}
		} else {
			if(!is.leaf(d)) {
				attr(d, "midpoint") = attr(dend2[[index]], "midpoint")
			}
		}
		return(d)
	})

	## s2 = max(s, s_child_none_leaf)
	dend = edit_node(dend, function(d, index) {
		s = attr(d, "score")
		if(is.leaf(d)) {
			attr(d, "score2") = s
		} else {
			l = vapply(d, is.leaf, TRUE)
			if(all(l)) {
				attr(d, "score2") = s
			} else {
				attr(d, "score2") = max(s, vapply(d[!l], function(x) attr(x, "score"), 0))
			}
		}

		return(d)
	})

	return(dend)
}

consensus_kmeans = function(mat, centers, km_repeats = 10) {
    partition_list = lapply(seq_len(km_repeats), function(i) {
        as.cl_hard_partition(kmeans(mat, centers))
    })
    partition_list = cl_ensemble(list = partition_list)
    partition_consensus = cl_consensus(partition_list)
    cl = as.vector(cl_class_ids(partition_consensus))
    if(length(unique(cl)) == 1) {
    	cl = partition_list[[1]]$.Data$cluster
    }
    cl
}

.cluster_mat = function(mat, dist_mat = dist(mat), .env, index = seq_len(nrow(mat)), 
	depth = 0, dend_index = NULL) {

	nr = nrow(mat)

	if(nr == 1) {
		.env$dend[[dend_index]] = index
		attributes(.env$dend[[dend_index]]) = list(
			members = 1,
			label = "",
			leaf = TRUE,
			height = depth,
			score = 0.5,
			index = dend_index,
			class = "dendrogram"
		)
		return(NULL)
	}
	if(nrow(mat) == 2) {
		cl = c(1, 2)
	} else {
		oe = try(suppressWarnings(cl <- consensus_kmeans(mat, centers = 2)), silent = TRUE)
		if(inherits(oe, "try-error")) {
			cl = rep(2, nr)
			cl[seq_len(ceiling(nr/2))] = 1
		}
	}
	l1 = cl == 1
	l2 = cl == 2
	if(mean(mat[l1, l1]) > mean(mat[l2, l2])) {
		l3 = l1
		l1 = l2
		l2 = l3
	}
	m11 = mat[l1, l1, drop = FALSE]
	m11 = m11[ lower.tri(m11) | upper.tri(m11)]
	x11 = .env$value_fun(m11)
	x12 = .env$value_fun(mat[l1, l2])
	x21 = .env$value_fun(mat[l2, l1])
	m22 = mat[l2, l2, drop = FALSE]
	m22 = m22[ lower.tri(m22) | upper.tri(m22)]
	x22 = .env$value_fun(m22)

	if(is.na(x11)) x11 = 0
	if(is.na(x22)) x22 = 0

	s = (x11 + x22)/(x11 + x12 + x21 + x22)
	if(is.na(s)) s = 1

	# sr = numeric(1)
	# for(i in seq_len(10)) {
	# 	clr = sample(cl, length(cl))
	# 	l1r = clr == 1
	# 	l2r = clr == 2
	# 	x11r = sample(as.vector(mat[l1r, l1r]), 1)
	# 	x12r = sample(as.vector(mat[l1r, l2r]), 1)
	# 	x21r = sample(as.vector(mat[l2r, l1r]), 1)
	# 	x22r = sample(as.vector(mat[l2r, l2r]), 1)

	# 	sr[i] = (x11r + x22r)/(x11r + x12r + x21r + x22r)
	# 	if(is.na(sr[i])) sr[i] = 1
	# }
	# sr = mean(sr)

	if(is.null(dend_index)) {
		.env$dend = list()
		attributes(.env$dend) = list(
			members = nrow(mat),
			height = depth,
			score = s,
			# random_score = sr,
			index = NULL,
			class = "dendrogram"
		)
	} else {
		.env$dend[[dend_index]] = list()
		attributes(.env$dend[[dend_index]]) = list(
			members = nrow(mat),
			height = depth,
			score = s,
			# random_score = sr,
			index = dend_index,
			class = "dendrogram"
		)
	}

	.cluster_mat(mat[l1, l1, drop = FALSE], as.dist(as.matrix(dist_mat)[l1, l1, drop = FALSE]), 
		.env, index = index[l1], depth = depth + 1, dend_index = c(dend_index, 1))
	.cluster_mat(mat[l2, l2, drop = FALSE], as.dist(as.matrix(dist_mat)[l2, l2, drop = FALSE]), 
		.env, index = index[l2], depth = depth + 1, dend_index = c(dend_index, 2))
}

cut_dend = function(dend, cutoff = 0.85, field = "score2", return = "cluster") {

	children_score = function(dend, field) {
		if(is.leaf(dend)) {
			-Inf
		} else {
			d1 = dend[[1]]
			d2 = dend[[2]]
			c(attr(d1, field), attr(d2, field))
		}
	}

	dont_split = function(dend, field, cutoff) {
		s = attr(dend, field)

		if(s >= cutoff) {
			return(FALSE)
		} else {
			s_children = children_score(dend, field)
			all(s_children < cutoff)
		}
	}

	# if the top node
	if(dont_split(dend, field, cutoff)) {
		dend2 = dendrapply(dend, function(d) {
			attr(d, "height") = 0
			d
		})
		if(plot) {
			plot(dend)
			box()
		}
		return(dend2)
	}

	dend2 = edit_node(dend, function(d, index) {
		if(dont_split(d, field, cutoff)) {
			attr(d, "height") = 0
		}
		d
	})

	## make sure all sub-nodes having height 0 if the node is 0 height
	is_parent_zero_height = function(index) {
		h = vapply(seq_along(index), function(i) {
			attr(dend2[[ index[seq_len(i)] ]], "height")
		}, 0)
		any(h == 0)
	}
	dend2 = edit_node(dend2, function(d, index) {
		if(is_parent_zero_height(index)) {
			attr(d, "height") = 0
			attr(d, "nodePar") = NULL
		}
		d
	})

	if(return == "dend") {
		dend2
	} else {
		cl = cutree(as.hclust(dend2), h = 0.1)

		od = order.dendrogram(dend2)
		cl = factor(cl, levels = unique(cl[od]))
		tb = table(cl)
		tb = sort(tb, decreasing = TRUE)

		map = structure(1:length(tb), names = names(tb))
		cl = map[as.character(cl)]
		return(cl)
	}
}

render_dend = function(dend, field = "score2", cutoff = 0.85, align_leaf = FALSE, depth = NULL) {

	if(!is.null(depth)) {
		dend = edit_node(dend, function(d, index) {
			if(length(index) + 1 > depth) {
				d = dendrapply(d, function(d) {
					attr(d, "height") = 0
					d
				})
			}
			return(d)
		})
		return(dend)
	}

	dend2 = cut_dend(dend, field = field, cutoff = cutoff, return = "dend")
	col_fun = colorRamp2(c(0.5, cutoff, 1), c("blue", "yellow", "red"))
	dend = edit_node(dend, function(d, index) {
		if(is.null(index)) {
			if(!is.leaf(d)) {
				s = attr(d, field)
				attr(d, "edgePar") = list(col = col_fun(s))
				
				if(attr(dend2, "height") > 0.5) {
					attr(d, "nodePar") = list(pch = 4, cex = 0.5)
				}
			} else {
				attr(d, "edgePar") = list(col = "#DDDDDD")
			}
		} else {
			if(!is.leaf(d)) {
				s = attr(d, field)
				attr(d, "edgePar") = list(col = col_fun(s))	
				if(attr(dend2[[index]], "height") > 0.5) {
					if(length(index) > 1) {
						if(attr(dend2[[index[-length(index)]]], "height") > 0.5) {
							attr(d, "nodePar") = list(pch = 4, cex = 0.5)
						}
					} else {
						if(attr(dend2, "height") > 0.5) {
							attr(d, "nodePar") = list(pch = 4, cex = 0.5)
						}
					}
				}
			} else {
				attr(d, "edgePar") = list(col = "#DDDDDD")
			}
		}

		if(align_leaf) {
			if(is.leaf(d)) {
				attr(d, "height") = 0
			}
		}
		return(d)
	})
	
	attr(dend, "col_fun") = col_fun
	dend
}

dend_env = new.env()
plot_binary_cut = function(mat, value_fun = median, cutoff = 0.85, dend = NULL, 
	border = "#404040", depth = NULL, show_heatmap_legend = TRUE, ...) {

	if(!requireNamespace("gridGraphics", quietly = TRUE)) {
		stop_wrap("Package gridGraphics should be installed.")
	}

	hash = digest::digest(list(mat, value_fun))
	if(is.null(dend)) {
		if(identical(hash, dend_env$hash)) {
			dend = dend_env$dend
		} else {
			dend_env$hash = NULL
		}
	} else {
		dend_env$dend = NULL
		dend_env$hash = NULL
	}
	if(is.null(dend)) {
		dend = cluster_mat(mat, value_fun = value_fun)
		dend_env$dend = dend
		dend_env$hash = hash
	}

	dend2 = render_dend(dend, cutoff = cutoff, depth = depth, ...)
	score_col_fun = attr(dend2, "col_fun")
	
	if(is.null(depth)) {
		cl = cut_dend(dend, cutoff = cutoff)
	} else {
		cl = cutree(as.hclust(dend2), h = 0.1)
	}
	
	dend2 = rev(dend2)
	f = function() {
		op =  par(c("mar","xpd"))
		par(mar = c(0, 0, 0, 0), xpd = NA)
		plot(rev(dend2), horiz = TRUE, axes = FALSE, ann = FALSE, ylim = c(0.5, nobs(dend2)+0.5), xaxs = "i", yaxs = "i")
		if(is.null(score_col_fun)) {
			text(par("usr")[2], par("usr")[4], qq("depth = @{depth}"), adj = c(1.2, 1))
		}
		par(op)
	}
	od = order.dendrogram(dend2)
	p2 = grid.grabExpr({
		cl = factor(cl, levels = unique(cl[od]))
		col_fun = colorRamp2(c(0, 1), c("white", "red"))
		ht = Heatmap(mat, name = "Similarity", col = col_fun,
			show_row_names = FALSE, show_column_names = FALSE,
			row_order = od, column_order = od,
			row_split = cl, column_split = cl,
			row_title = NULL, column_title = NULL,
			row_gap = unit(0, "mm"), column_gap = unit(0, "mm"),
			border = border, show_heatmap_legend = show_heatmap_legend) + NULL

		if(is.null(score_col_fun)) {
			draw(ht)
		} else {
			draw(ht, heatmap_legend_list = list(Legend(title = "Score", col_fun = score_col_fun)))
		}
	})

	grid.newpage()
	pushViewport(viewport(x = 0, width = 0.3, just = "left"))
	pushViewport(viewport(height = unit(1, "npc") - unit(4, "mm"), x = unit(2, "mm"), width = unit(1, "npc") - unit(2, "mm"), just = "left"))
	gridGraphics::grid.echo(f, newpage = FALSE)
	popViewport(2)
	pushViewport(viewport(x = 0.3, width = 0.7, just = "left"))
	grid.draw(p2)
	popViewport()

}

# == title
# Cluster functional terms by binary cutting the similarity matrix
#
# == param
# -mat A similarity matrix.
# -value_fun Value function to calculate the score for each node in the dendrogram.
# -cutoff The cutoff for splitting the dendrogram.
# -n_run If the value is larger than one, `binary_cut` is executed multiple times
#      and generates a consensus clustering.
#
# == value
# A vector of cluster labels (in numeric). 
#
# == example
# mat = readRDS(system.file("extdata", "similarity_mat.rds", package = "simplifyEnrichment"))
# binary_cut(mat)
binary_cut = function(mat, value_fun = median, cutoff = 0.85, n_run = 1) {
	if(n_run == 1) {
		dend = cluster_mat(mat, value_fun = value_fun)
		cl = cut_dend(dend, cutoff)
		return(unname(cl))
	} else {
		partition_list = lapply(seq_len(n_run), function(i) {
			dend = cluster_mat(mat, value_fun = value_fun)
			cl = cut_dend(dend, cutoff)
            as.cl_hard_partition(cl)
        })
        partition_list = cl_ensemble(list = partition_list)
        partition_consensus = cl_consensus(partition_list)
        as.vector(cl_class_ids(partition_consensus)) 
	}
}


