
dend_max_depth = function(dend) {
	max(unlist(dend_node_apply(dend, function(d, index) {
		length(index) + 1
	})))
}

# == title
# Area above the eCDF curve
#
# == param
# -x A vector of similarity values.
#
# == details
# Denote F(x) as the eCDF (empirical Cumulative Distribution Function) of the similarity vector ``x``, this function calculates
# the area above the eCDF curve, which is 1 - \\int_0^1 F(x)dx.
#
# == value
# A numeric value.
area_above_ecdf = function(x) { 
	f = ecdf(x)
	i = 1:100/100
	1 - sum(1/100*f(i))
}

# cluster the similarity matrix and assign scores to nodes
# two scores attached to each ndoe
# - score: the original score
# - score2: the score which have checked the children nodes' scores
cluster_mat = function(mat, value_fun = area_above_ecdf, partition_fun = partition_by_pam, cutoff = 0.85, return_dend = FALSE) {
	dend = NULL

	dend_ind_list = list(NULL)
	mat_ind_list = list(seq_len(nrow(mat)))

	# breadth-first
	while(1) {

		if(length(mat_ind_list) == 0) break

		mat_ind_list2 = list()
		dend_ind_list2 = list()

		n = length(mat_ind_list)
		for(i in seq_len(n)) {
			mat_ind = mat_ind_list[[i]]
			dend_ind = dend_ind_list[[i]]

			lt = .cluster_mat(mat[mat_ind, mat_ind, drop = FALSE], value_fun, partition_fun)
			lt$attr$height = length(dend_ind)
			lt$attr$members = length(mat_ind)

			if(length(mat_ind) == 1) {
				if(is.null(dend_ind)) {
					dend = mat_ind
				} else {
					dend[[dend_ind]] = mat_ind
				}
			} else {
				if(is.null(dend_ind)) {
					dend = list()
				} else {
					dend[[dend_ind]] = list()
				}
			}
			if(is.null(dend_ind)) {
				attributes(dend) = lt$attr
			} else {
				attributes(dend[[dend_ind]]) = lt$attr
			}

			if(length(mat_ind) > 1) {
				mat_ind_list2 = c(mat_ind_list2, list(mat_ind[lt$ind1], mat_ind[lt$ind2]))
				dend_ind_list2 = c(dend_ind_list2, list(c(dend_ind, 1), c(dend_ind, 2)))
			}
		}

		mat_ind_list = mat_ind_list2
		dend_ind_list = dend_ind_list2
	}

	dend = rev(dend)

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

	dend = assign_score2(dend, cutoff)

	hash = digest::digest(list(mat, value_fun, partition_fun, cutoff))
	dend_env$dend = dend
	dend_env$hash = hash

	cl = cut_dend(dend, cutoff)
	cl = as.numeric(as.vector(unname(cl)))
	
	if(return_dend)
	  return(list(cl = cl, dend = dend))
	return(cl)
}

cluster_mat2 = function(mat, value_fun = median, partition_fun = partition_by_pam, cutoff = 0.85, return_dend = FALSE) {

	dend = NULL
	dend_env$dend = NULL

	dend_ind_list = list(NULL)
	mat_ind_list = list(seq_len(nrow(mat)))

	# breadth-first
	while(1) {

		if(length(mat_ind_list) == 0) break

		mat_ind_list2 = list()
		dend_ind_list2 = list()

		merged_dend_ind_list = list()

		n = length(mat_ind_list)
		for(i in seq_len(n)) {
			mat_ind = mat_ind_list[[i]]
			dend_ind = dend_ind_list[[i]]

			do_partitioning = TRUE
			if(length(mat_ind) == 1) {
				lt = .cluster_mat(mat[mat_ind, mat_ind, drop = FALSE], value_fun, partition_fun)
			} else {
				if(length(dend_ind) >= 3) {
					# parent node
					pd = dend[[dend_ind[seq(1, length(dend_ind)-2)]]]
					
					s2 = max(attr(pd[[1]], "score"), attr(pd[[2]], "score"))
					s = attr(pd, "score")
					if( min(nobs(pd[[1]]), nobs(pd[[2]]))/nobs(pd) < 0.1 ) {
						
					} else if(s > cutoff*0.8 && s < s2) {
						s = s2
					}

					if(s < cutoff) {
						do_partitioning = FALSE
						merged_dend_ind_list[[ length(merged_dend_ind_list) + 1 ]] = dend_ind[seq(1, length(dend_ind)-2)]
					} else {
						lt = .cluster_mat(mat[mat_ind, mat_ind, drop = FALSE], value_fun, partition_fun)
					}
				} else {
					lt = .cluster_mat(mat[mat_ind, mat_ind, drop = FALSE], value_fun, partition_fun)
				}
			}

			if(do_partitioning) {
				lt$attr$height = length(dend_ind)
				lt$attr$members = length(mat_ind)

				if(length(mat_ind) == 1) {
					if(is.null(dend_ind)) {
						dend = mat_ind
					} else {
						dend[[dend_ind]] = mat_ind
					}
				} else {
					if(is.null(dend_ind)) {
						dend = list()
					} else {
						dend[[dend_ind]] = list()
					}
				}
				if(is.null(dend_ind)) {
					attributes(dend) = lt$attr
				} else {
					attributes(dend[[dend_ind]]) = lt$attr
				}

				if(length(mat_ind) > 1) {
					mat_ind_list2 = c(mat_ind_list2, list(mat_ind[lt$ind1], mat_ind[lt$ind2]))
					dend_ind_list2 = c(dend_ind_list2, list(c(dend_ind, 1), c(dend_ind, 2)))
				}
			} else {
				# construct dendrogram for mat[mat_ind, mat_ind, drop = FALSE]
				dend[[dend_ind]] = mat_ind
				attributes(dend[[dend_ind]]) = list(
					label = "",
					leaf = TRUE,
					members = 1,
					height = 0,
					score = 0.5,
					class = "dendrogram"
				)
			}
		}

		mat_ind_list = mat_ind_list2
		dend_ind_list = dend_ind_list2

		merged_dend_ind_list = unique(merged_dend_ind_list)
		for(i in seq_along(merged_dend_ind_list)) {
			h = attr(dend[[ merged_dend_ind_list[[i]] ]], "height")
			dend[[ merged_dend_ind_list[[i]] ]] = unlist(dend[[ merged_dend_ind_list[[i]] ]])
			attributes(dend[[ merged_dend_ind_list[[i]] ]]) = list(
				label = "",
				leaf = TRUE,
				members = 1,
				height = h,
				score = 0.5,
				class = "dendrogram"
			)
		}
	}

	dend = rev(dend)

	lt = dend_node_apply(dend, function(d, ind) {
		if(is.leaf(d)) {
			as.vector(d) 
		} else {
			NULL
		}
	})
	lt = lt[!sapply(lt, is.null)]
	cl = numeric(nrow(mat))
	for(i in seq_along(lt)) {
		cl[lt[[i]]] = i
	}
	
	cl = as.numeric(as.vector(unname(cl)))
	if(return_dend)
	  return(list(cl = cl, dend = dend))
	return(cl)
}

## bottom-up
assign_score2 = function(d, cutoff) {
	if(is.leaf(d)) {
		attr(d, "score2") = 0.5
		return(d)
	}
	if(length(d) == 0) {
		attr(d, "score2") = 0.5
		return(d)
	}
	# check the two child node
	if(is.null(attr(d[[1]], "score2"))) {
		d[[1]] = assign_score2(d[[1]], cutoff)
	}
	if(is.null(attr(d[[2]], "score2"))) {
		d[[2]] = assign_score2(d[[2]], cutoff)
	}

	s2 = max(attr(d[[1]], "score2"), attr(d[[2]], "score2"))
	# the node is assigned with its child nodes' score only if
	# 1. the size of the two children should not be too extremely different
	# 2. the score should not be too smaller than its child nodes' 
	s = attr(d, "score")
	if( min(nobs(d[[1]]), nobs(d[[2]]))/nobs(d) < 0.1 ) {
		attr(d, "score2")  = s
	} else if(s > cutoff*0.8 && s < s2) {
		attr(d, "score2") = s2
	} else {
		attr(d, "score2")  = s
	}
	return(d)
}

# assign current node and do partitioning to generate indices of two subsets
.cluster_mat = function(mat, value_fun, partition_fun) {
	nr = nrow(mat)

	if(nr == 1) {
		attr = list(
			label = "",
			leaf = TRUE,
			score = 0.5,
			class = "dendrogram"
		)
		return(list(attr = attr))
	}

	if(nrow(mat) == 2) {
		cl = c(1, 2)
	} else {
		oe = try(suppressWarnings(cl <- partition_fun(mat)), silent = TRUE)
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
	if(nrow(m11) == 1) {
		x11 = 1
		x12 = value_fun(mat[l1, l2])
		x21 = value_fun(mat[l2, l1])
	} else {
		m11 = m11[ lower.tri(m11) | upper.tri(m11)]
		x11 = value_fun(m11)
		x12 = value_fun(mat[l1, l2])
		x21 = value_fun(mat[l2, l1])
	}
	m22 = mat[l2, l2, drop = FALSE]
	if(nrow(m22) == 1) {
		x22 = 1
	} else {
		m22 = m22[ lower.tri(m22) | upper.tri(m22)]
		x22 = value_fun(m22)
	}
	
	s = (x11 + x22)/(x11 + x12 + x21 + x22)
	if(is.na(s)) s = 1

	attr = list(
		members = NA,
		height = NA,
		score = s,
		class = "dendrogram"
	)
	return(list(attr = attr, ind1 = which(l1), ind2 = which(l2)))

}

cut_dend = function(dend, cutoff = 0.85, field = "score2", return = "cluster") {

	## add a split attribute on each node
	assign_child_node = function(d, split = NULL) {
		if(is.null(split)) {
			if(attr(d, field) >= cutoff) {
				attr(d, "split") = TRUE
				if(!is.leaf(d)) {
					d[[1]] = assign_child_node(d[[1]], NULL)
					d[[2]] = assign_child_node(d[[2]], NULL)
				}
			} else {
				attr(d, "split") = FALSE
				if(!is.leaf(d)) {
					d[[1]] = assign_child_node(d[[1]], FALSE)
					d[[2]] = assign_child_node(d[[2]], FALSE)
				}
			}
			
		} else {
			attr(d, "split") = split
			if(!is.leaf(d)) {
				d[[1]] = assign_child_node(d[[1]], split)
				d[[2]] = assign_child_node(d[[2]], split)
			}
		}
		d
	}

	dend = assign_child_node(dend)

	# if the top node
	if(!attr(dend, "split")) {
		dend2 = dendrapply(dend, function(d) {
			attr(d, "height") = 0
			d
		})

		if(return == "dend") {
			return(dend2)
		} else {
			return(rep(1, nobs(dend)))
		}
	}

	dend2 = edit_node(dend, function(d, index) {
		if(!attr(d, "split")) {
			attr(d, "height") = 0
		}
		d
	})

	if(return == "dend") {
		dend2
	} else {
		cl = cutree(as.hclust(dend2), h = 0.1)
		cl = factor(cl, levels = unique(cl[order.dendrogram(dend)]))
		return(cl)
	}
}

render_dend = function(dend, field = "score2", cutoff = 0.85, align_leaf = FALSE, 
	depth = NULL, add_label = FALSE) {

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

	# note: in dend2 branches that are not split all have height 0
	dend2 = cut_dend(dend, field = field, cutoff = cutoff, return = "dend")
	col_fun = colorRamp2(c(0.5, cutoff, 1), c("blue", "yellow", "red"))
	dend = edit_node(dend, function(d, index) {
		if(is.null(index)) {
			if(!is.leaf(d)) {
				s = attr(d, field)
				attr(d, "edgePar") = list(col = col_fun(s))
				
				if(attr(dend2, "height") > 0.5) {
					attr(d, "nodePar") = list(pch = 4, cex = 0.5)
					if(add_label) attr(d, "edgetext") = sprintf("%.2f", attr(d, "score"))
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
							if(add_label) attr(d, "edgetext") = sprintf("%.2f", attr(d, "score"))
						}
					} else {
						if(attr(dend2, "height") > 0.5) {
							attr(d, "nodePar") = list(pch = 4, cex = 0.5)
							if(add_label) attr(d, "edgetext") = sprintf("%.2f", attr(d, "score"))
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

# == title
# Visualize the process of binary cut
#
# == param
# -mat The similarity matrix.
# -value_fun A function that calculates the scores for the four submatrices on a node.
# -cutoff The cutoff for splitting the dendrogram.
# -partition_fun A function to split each node into two groups. Pre-defined functions
#                in this package are `partition_by_kmeanspp`, `partition_by_pam` and `partition_by_hclust`.
# -dend A dendrogram object, used internally.
# -depth Depth of the recursive binary cut process.
# -dend_width Width of the dendrogram on the plot.
# -show_heatmap_legend Whether to show the heatmap legend.
# -... Other arguments. 
#
# == details
# After the functions which perform clustering are executed, such as `simplifyGO` or
# `binary_cut`, the dendrogram is temporarily saved and `plot_binary_cut` directly
# uses this dendrogram.
#
# == example
# \donttest{
# mat = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds", 
#     package = "simplifyEnrichment"))
# plot_binary_cut(mat, depth = 1)
# plot_binary_cut(mat, depth = 2)
# plot_binary_cut(mat)
# }
plot_binary_cut = function(mat, value_fun = area_above_ecdf, cutoff = 0.85, 
	partition_fun = partition_by_pam, dend = NULL, dend_width = unit(3, "cm"),
	depth = NULL, show_heatmap_legend = TRUE, ...) {

	check_pkg("gridGraphics", bioc = FALSE)

	hash = digest::digest(list(mat, value_fun, partition_fun, cutoff))

	# if it was already generated
	if(is.null(dend)) {
		if(identical(hash, dend_env$hash)) {
			dend = dend_env$dend
			if(se_opt$verbose) {
				message("use the cached dendrogram.")
			}
		} else {
			dend_env$hash = NULL
		}
	} else {
		dend_env$dend = NULL
		dend_env$hash = NULL
	}
	if(is.null(dend)) {
		if(se_opt$verbose) {
			message("create a new dendrogram.")
		}
		dend = cluster_mat(mat, value_fun = value_fun, partition_fun = partition_fun, 
		                   cutoff = cutoff, return_dend = TRUE)$dend
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
	
	od = order.dendrogram(dend2)
	
	col_fun = colorRamp2(c(0, quantile(mat, 0.95)), c("white", "red"))
	ht = Heatmap(mat, name = "Similarity", col = col_fun,
		show_row_names = FALSE, show_column_names = FALSE,
		cluster_rows = dend2, row_dend_width = dend_width,
		column_order = od,
		row_title = NULL, column_title = NULL,
		show_heatmap_legend = show_heatmap_legend,
		use_raster = TRUE) + NULL

	if(is.null(score_col_fun)) {
		draw(ht)
	} else {
		at = c(0.5, cutoff, 1)
		labels = c(0.5, paste0(cutoff, " (cutoff)"), 1)
		draw(ht, heatmap_legend_list = list(Legend(title = "Score", col_fun = score_col_fun, at = at, labels = labels)))
	}
	decorate_heatmap_body("Similarity", {
		grid.rect(gp = gpar(fill = NA, col = "#404040"))
		cl = factor(cl, levels = unique(cl[od]))
		tbcl = table(cl)
		ncl = length(cl)
		x = cumsum(c(0, tbcl))/ncl
		grid.segments(x, 0, x, 1, default.units = "npc", gp = gpar(col = "#404040"))
		grid.segments(0, 1 - x, 1, 1 - x, default.units = "npc", gp = gpar(col = "#404040"))
	})
}

# == title
# Cluster functional terms by recursively binary cutting the similarity matrix
#
# == param
# -mat A similarity matrix.
# -value_fun A function that calculates the scores for the four submatrices on a node.
# -partition_fun A function to split each node into two groups. Pre-defined functions
#                in this package are `partition_by_kmeanspp`, `partition_by_pam`  and `partition_by_hclust`.
# -cutoff The cutoff for splitting the dendrogram.
# -try_all_partition_fun Different ``partition_fun`` gives different clusterings. If the vaule
#      of ``try_all_partition_fun`` is set to ``TRUE``, the similarity matrix is clustered by three
#      partitioning method: `partition_by_pam`, `partition_by_kmeanspp` and `partition_by_hclust`.
#      The clustering with the highest difference score is finally selected as the final clustering.
# -partial Whether to generate the complete clustering or the clustering stops when sub-matrices
#     cannot be split anymore.
# -return_dend Whether to return the generated dendogram together with the cluster labels or not. 
#
# == value
# If ``return_dend`` is set to ``FALSE`` return a vector of cluster labels (in numeric). Otherwise 
# returns a named list with a vector of cluster labels (in numeric) and the generated dendrogram
#
# == example
# mat = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds",
#     package = "simplifyEnrichment"))
# binary_cut(mat)
binary_cut = function(mat, value_fun = area_above_ecdf, partition_fun = partition_by_pam,
	cutoff = 0.85, try_all_partition_fun = FALSE, partial = FALSE, return_dend = FALSE) {

	if(try_all_partition_fun) {
	  clt = list(
	    by_pam = binary_cut(mat, value_fun = value_fun, partition_fun = partition_by_pam,
	                        cutoff = cutoff, try_all_partition_fun = FALSE, partial = partial),
	    by_kmeanspp = binary_cut(mat, value_fun = value_fun, partition_fun = partition_by_kmeanspp,
	                             cutoff = cutoff, try_all_partition_fun = FALSE, partial = partial),
	    by_hclust = binary_cut(mat, value_fun = value_fun, partition_fun = partition_by_hclust,
	                           cutoff = cutoff, try_all_partition_fun = FALSE, partial = partial)
		)
		i = which.max(sapply(clt, function(cl) difference_score(mat, cl)))

		# qqcat("@{names(clt)[i]} gives the highest difference score\n")
		if(length(i) == 0) i = 1
		return(clt[[i]])
	}

  cluster_method <- ifelse(partial, cluster_mat2, cluster_mat)
  clustering <- do.call(cluster_method, 
                        list(mat = mat,value_fun = value_fun, 
                             partition_fun = partition_fun, cutoff = cutoff,
                             return_dend = return_dend))
  return(clustering)
}

# == title
# Select the cutoff for binary cut
#
# == param
# -mat A similarity matrix.
# -cutoff A list of cutoffs to test. Note the range of the cutoff values should be inside [0.5, 1].
# -verbose Whether to print messages.
# -... Pass to `binary_cut`. Of node, ``return_dend`` is set to ``FALSE`` and cannot be passed
#
# == details
# Binary cut is applied to each of the cutoff and the clustering results are evaluated by following metrics:
#
# - difference score, calculated by `difference_score`.
# - number of clusters.
# - block mean, which is the mean similarity in the blocks in the diagonal of the heatmap.
#
# == example
# \donttest{
# mat = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds",
#     package = "simplifyEnrichment"))
# select_cutoff(mat)
# }
select_cutoff = function(mat, cutoff = seq(0.6, 0.98, by = 0.01), verbose = TRUE, ...) {

	cutoff = cutoff[cutoff >= 0.5 & cutoff <= 1]

	s1 = s2 = s3 = s4 = NULL
	for(i in seq_along(cutoff)) {
		if(verbose) message(qq("@{i}/@{length(cutoff)}, cutoff = @{cutoff[i]}..."))
		cl = binary_cut(mat, cutoff = cutoff[i], return_dend = FALSE, ...)
		s1[i] = difference_score(mat, cl)
		tb = table(cl)
		s2[i] = length(tb)
		s3[i] = sum(tb >= 5)
		s4[i] = block_mean(mat, cl)
	}

	check_pkg("cowplot", bioc = FALSE)
	check_pkg("ggplot2", bioc = FALSE)

	suppressWarnings(
		p1 <- ggplot2::ggplot(data = NULL, ggplot2::aes(x = cutoff, y = s1)) +
		ggplot2::geom_point() + ggplot2::ylab("Difference score") +
		ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank())
	)

	df1 = data.frame(cutoff, s2); colnames(df1) = c("method", "value")
	df2 = data.frame(cutoff, s3); colnames(df2) = c("method", "value")
	df1$type = "All sizes"
	df2$type = "Size >= 5"
	df = rbind(df1, df2)
	suppressWarnings(
		p2 <- ggplot2::ggplot(df, ggplot2::aes(x = df$method, y = df$value, col = df$type, fill = df$type)) +
		ggplot2::geom_point() + ggplot2::ylab("Number of clusters") + ggplot2::labs(col = "Type", fill = "Type") +
		ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank())
	)

	suppressWarnings(
		p3 <- ggplot2::ggplot(data = NULL, ggplot2::aes(x = cutoff, y = s4)) +
		ggplot2::geom_point() + ggplot2::ylab("Block mean") + ggplot2::xlab("Cutoff") +
		ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
	)

	suppressWarnings(print(cowplot::plot_grid(p1, p2, p3, nrow = 3, align = "v", axis = "lr", rel_heights = c(1, 1, 1.3))))
}

