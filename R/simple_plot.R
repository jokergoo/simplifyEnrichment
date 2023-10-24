
# == title
# A simplified way to visualize enrichment in GO clusters
#
# == param
# -go_id A vector of GO IDs.
# -value A list of numeric value associate with ``go_id``. We suggest to use -log10(p.adjust) or -log2(fold enrichment) as the values.
# -aggregate Function to aggregate values in each GO cluster.
# -method Method for clustering the matrix. See `cluster_terms`.
# -control A list of parameters for controlling the clustering method, passed to `cluster_terms`.
# -verbose Whether to print messages.
# -axis_label X-axis label.
# -title Title for the whole plot.
# -legend_title Title for the legend.
# -min_term Minimal number of functional terms in a cluster. All the clusters
#     with size less than ``min_term`` are all merged into one separated cluster in the heatmap.
# -stat Type of value for mapping to the font size of keywords in the word clouds. There are two options:
#       "count": simply number of keywords; "pvalue": enrichment on keywords is performed (by fisher's exact test) and -log10(pvalue) is used to map to font sizes.
# -min_stat Minimal value for ``stat`` for selecting keywords.
# -exclude_words Words that are excluded in the word cloud.
# -max_words Maximal number of words visualized in the word cloud.
# -word_cloud_grob_param A list of graphic parameters passed to `word_cloud_grob`.
# -fontsize_range The range of the font size. The value should be a numeric vector with length two.
#       The font size interpolation is linear.
# -bg_gp Graphics parameters for controlling word cloud annotation background.
#
# == details
# There are several other ways to specify GO IDs and the associated values.
#
# 1. specify ``value`` as a named vector where GO IDs are the names.
# 2. specify ``value`` as a list of numeric named vectors. In this case, ``value`` contains multiple enrichment results.
summarizeGO = function(go_id, value = NULL, aggregate = mean, 
	method = "binary_cut", control = list(), verbose = TRUE, 
	axis_label = "Value", title = "", legend_title = axis_label,

	min_term = round(nrow(mat)*0.01), 
	stat = "pvalue", 
	min_stat = ifelse(stat == "count", 5, 0.05),
	exclude_words = character(0), 
	max_words = 6,
	word_cloud_grob_param = list(), 
	fontsize_range = c(4, 16), 
	bg_gp = gpar(fill = "#DDDDDD", col = "#AAAAAA")

	) {

	if(missing(go_id)) {
		if(is.atomic(value)) {
			go_id = names(value)
		} else if(is.list(value)) {
			go_id = unique(unlist(lapply(value, names)))
			if(length(go_id) == 0) {
				stop("If `value` is set as list, each element vector should be a named vector.")
			}
			vv = matrix(0, nrow = length(go_id), ncol = length(value))
			rownames(vv) = go_id
			colnames(vv) = names(value)
			for(i in seq_along(value)) {
				vv[names(value[[i]]), i] = value[[i]]
			}
			value = vv
		}
	}

	if(is.null(value)) {
		value = rep(1, length(go_id))
		aggregate = sum
		axis_label = "Number of terms"
	}

	if(is.vector(value)) {
		value = cbind(value)
	}
	rownames(value) = go_id
	if(is.null(colnames(value))) {
		colnames(value) = paste0("C", seq_len(ncol(value)))
	}
	mat = GO_similarity(go_id)

	cl = do.call(cluster_terms, list(mat = mat, method = method, verbose = verbose, control = control))
	
	value = value[rownames(mat), , drop = FALSE]
	go_id = rownames(mat)

	cl = as.vector(cl)
	cl_tb = table(cl)
	cl[as.character(cl) %in% names(cl_tb[cl_tb < min_term])] = 0
	cl = factor(cl, levels = c(setdiff(sort(cl), 0), 0))

	l = cl != 0
	cl = cl[l]
	cl = as.vector(cl)
	go_id = go_id[l]
	value = value[l, , drop = FALSE]

	align_to = split(seq_along(cl), cl)
	go_id = split(go_id, cl)
	n = length(align_to)

	v2 = tapply(seq_along(cl), cl, function(ind) {
		apply(value[ind, , drop = FALSE], 2, aggregate)
	}, simplify = FALSE)
	v2 = do.call(rbind, v2)

	gbl = anno_word_cloud_from_GO(align_to, go_id, return_gbl = TRUE,
		stat = stat, min_stat = min_stat,
		exclude_words = exclude_words, max_words = max_words, word_cloud_grob_param = word_cloud_grob_param, 
		fontsize_range = fontsize_range, bg_gp = bg_gp)


	gbl_h = lapply(gbl, function(x) convertHeight(grobHeight(x), "cm") + unit(10, "pt"))
	gbl_h = do.call(unit.c, gbl_h)

	gbl_w = lapply(gbl, function(x) convertWidth(grobWidth(x), "cm"))
	gbl_w = do.call(unit.c, gbl_w)
	gbl_w = max(gbl_w) + unit(10, "pt")

	gap = rep( (unit(1, "npc") - sum(gbl_h))/(n-1), n - 1)
	gap = unit.c(unit(0, "mm"), gap)

	if(ncol(v2) > 1) {
		size_fun = generate_size_fun(range(v2), c(2, 20))
		size_breaks = grid.pretty(range(v2), 3)
		lgd = Legend(title = legend_title, at = size_breaks, type = "points", 
			size = unit(size_fun(size_breaks), "pt"), pch = 16, legend_gp = gpar(col = "#888888"),
			row_gap = unit(6, "pt"), background = "white")
	}
	
	grid.newpage()
	pushViewport(viewport(x = unit(5, "mm"), y = unit(1.6, "cm"), width = unit(1, "npc") - unit(1, "cm"), height = unit(1, "npc") - unit(3, "cm"), just = c("left", "bottom")))
	for(i in seq_along(gbl)) {
		y = sum(gbl_h[seq_len(i)]) + sum(gap[seq_len(i)]) - gbl_h[i]*0.5
		pushViewport(viewport(x = 0, y = y, 
			width = gbl_w, height = gbl_h[i], just = c("left")))
		grid.rect(gp = gpar(fill = "#EEEEEE"))
		gb = gbl[[i]]
	    gb$vp$x = gb$vp$width*0.5 + unit(5, "pt")
	    grid.draw(gb)
		popViewport()

		if(ncol(v2) == 1) {
			if(all(v2 > 0)) {
				pushViewport(viewport(x = gbl_w + unit(5, "pt"), y = y, width = unit(1, "npc") - gbl_w, height = gbl_h[i],
					xscale = c(0, max(v2)), just = c("left")))
			} else {
				pushViewport(viewport(x = gbl_w + unit(5, "pt"), y = y, width = unit(1, "npc") - gbl_w, height = gbl_h[i],
					xscale = range(v2), just = c("left")))
			}

			grid.rect(0, 0.5, width = unit(v2[i], "native"), height = unit(6, "mm"), just = c("left"),
				gp = gpar(fill = "#CCCCCC"))

			if(i == 1) {
				gb = xaxisGrob(gp = gpar(fontsize = 8))
				grid.draw(gb)
				grid.text(axis_label, 0.5, -unit(25, "pt"), just = "top")
			}
			popViewport()
		} else {
			pushViewport(viewport(x = gbl_w + unit(5, "pt"), y = y, width = unit(1, "npc") - gbl_w - grobWidth(lgd@grob) - unit(10, "pt"), height = gbl_h[i],
					xscale = c(0.5, ncol(v2) + 0.5), just = c("left")))
			grid.points(x = unit(1:ncol(v2), "native"), y = rep(0.5, ncol(v2)), pch = 16, gp = gpar(col = "#888888"), size = unit(size_fun(v2[i, ]), "pt"))
			if(i == 1) {
				grid.text(colnames(value), x = unit(1:ncol(v2), "native"), y = rep(unit(-4, "pt"), ncol(v2)), just = "top", gp = gpar(fontsize = 8))
				grid.text(axis_label, 0.5, -unit(20, "pt"), just = "top")
			}
			popViewport()
		}

	}
	grid.text(title, y = unit(1, "npc") + unit(10, "pt"), just = "bottom", gp = gpar(fontsize = 14))

	if(ncol(v2) > 1) {
		draw(lgd, x = unit(1, "npc"), y = unit(0.5, "npc"), just = c("right"))
	}
	popViewport()

}

generate_size_fun = function(rg, size) {
	function(x) {
		x[x < rg[1]] = rg[1]
		x[x > rg[2]] = rg[2]

		(x - rg[1]) * (size[2] - size[1])/(rg[2] - rg[1]) + size[1]
	}
}