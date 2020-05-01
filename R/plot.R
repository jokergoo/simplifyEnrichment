
# == title
# Visualize the similarity matrix and the clustering
#
# == param
# -mat A similarity matrix.
# -cl Cluster labels inferred from the similarity matrix, e.g. from `cluster_terms` or `binary_cut`.
# -dend Used internally.
# -draw_word_cloud Whether to draw the word clouds.
# -term The full name or the description of the corresponding GO IDs. 
# -min_term Minimal number of functional terms in a cluster. All the clusters
#     with size less than ``min_term`` are all merged into one separated cluster in the heatmap.
# -order_by_size Whether to reorder clusters by their sizes. The cluster
#      that is merged from small clusters (size < ``min_term``) is always put to the bottom of the heatmap.
# -exclude_words Words that are excluded in the word cloud.
# -max_words Maximal number of words visualized in the word cloud.
# -word_cloud_grob_param A list of graphic parameters passed to `word_cloud_grob`.
# -fontsize_range The range of the font size. The value should be a numeric vector with length two.
#       The minimal font size is mapped to word frequency value of 1 and the maximal font size is mapped
#       to the maximal word frequency. The font size interlopation is linear.
# -... Other arguments passed to `ComplexHeatmap::draw,HeatmapList-method`.
#
# == value
# A `ComplexHeatmap::HeatmapList-class` object.
#
# == example
# mat = readRDS(system.file("extdata", "similarity_mat.rds", package = "simplifyEnrichment"))
# cl = binary_cut(mat)
# ht_clusters(mat, cl, word_cloud_grob_param = list(max_width = 80))
ht_clusters = function(mat, cl, dend = NULL, 
	draw_word_cloud = TRUE, term = NULL, min_term = 5, order_by_size = FALSE,
	exclude_words = character(0), max_words = 10,
	word_cloud_grob_param = list(), fontsize_range = c(4, 16), ...) {

	if(inherits(cl, "try-error")) {
		grid.newpage()
		pushViewport(viewport())
		grid.text("Clustering has an error.")
		popViewport()
		return(invisible(NULL))
	}

	cl = as.character(cl)
	cl_tb = table(cl)
	cl[cl %in% names(cl_tb[cl_tb < min_term])] = "0"
	cl = factor(cl, levels = c(setdiff(names(sort(table(cl), decreasing = TRUE)), "0"), "0"))

	if(!is.null(dend)) {
		ht = Heatmap(mat, col = colorRamp2(c(0, 1), c("white", "red")),
			name = "Similarity", column_title = NULL,
			show_row_names = FALSE, show_column_names = FALSE,
			cluster_rows = dend, cluster_columns = dend, 
			show_row_dend = TRUE, show_column_dend = FALSE,
			row_dend_width = unit(4, "cm"),
			border = "#404040", row_title = NULL)
	} else {
		ht = Heatmap(mat, col = colorRamp2(c(0, 1), c("white", "red")),
			name = "Similarity", column_title = NULL,
			show_row_names = FALSE, show_column_names = FALSE,
			show_row_dend = FALSE, show_column_dend = FALSE,
			cluster_row_slices = !order_by_size, 
			cluster_column_slices = !order_by_size,
			row_split = cl, column_split = cl, 
			border = "#404040", row_title = NULL,
			row_gap = unit(0, "mm"), column_gap = unit(0, "mm"))

		if(draw_word_cloud) {
			go_id = rownames(mat)
			if(!is.null(term)) {
				if(length(term) != length(go_id)) {
					stop_wrap("Length of `term` should be the same as the nrow of `mat`.")
				}
				id2term = structure(term, names = go_id)
			}
			keywords = tapply(go_id, cl, function(term_id) {
				if(is.null(term)) {
					suppressMessages(suppressWarnings(df <- count_word(term_id, exclude_words = exclude_words)))	
				} else {
					suppressMessages(suppressWarnings(df <- count_word(term = id2term[term_id], exclude_words = exclude_words)))
				}
				df = df[df$freq > 1, , drop = FALSE]
				if(nrow(df) > max_words) {
					df = df[order(df$freq, decreasing = TRUE)[1:max_words], ]
				}
				df
			})
			keywords = keywords[sapply(keywords, nrow) > 0]

			align_to = split(1:nrow(mat), cl)
			align_to = align_to[names(align_to) != "0"]
			align_to = align_to[names(align_to) %in% names(keywords)]

			word_cloud_grob_param = word_cloud_grob_param[setdiff(names(word_cloud_grob_param), c("text", "fontsize"))]
			gbl = lapply(names(align_to), function(nm) {
				kw = rev(keywords[[nm]][, 1])
				freq = rev(keywords[[nm]][, 2])
				fontsize = scale_fontsize(freq, rg = c(1, max(10, freq)), fs = fontsize_range)

				lt = c(list(text = kw, fontsize = fontsize), word_cloud_grob_param)
				do.call(word_cloud_grob, lt)
			})
			names(gbl) = names(align_to)

			margin = unit(8, "pt")
			gbl_h = lapply(gbl, function(x) convertHeight(grobHeight(x), "cm") + margin)
			gbl_h = do.call(unit.c, gbl_h)

			gbl_w = lapply(gbl, function(x) convertWidth(grobWidth(x), "cm"))
			gbl_w = do.call(unit.c, gbl_w)
			gbl_w = max(gbl_w) + margin

			panel_fun = function(index, nm) {
				pushViewport(viewport())
				grid.rect(gp = gpar(fill = "#DDDDDD", col = NA))
				grid.lines(c(0, 1, 1, 0), c(0, 0, 1, 1), gp = gpar(col = "#AAAAAA"), default.units = "npc")
			    pushViewport(viewport(width = unit(1, "npc") - margin, height = unit(1, "npc") - margin))
			    gb = gbl[[nm]]
			    gb$vp$x = gb$vp$width*0.5
			    gb$vp$y = gb$vp$height*0.5
			    grid.draw(gb)
			    popViewport()
			    popViewport()
			}

			ht = ht + rowAnnotation(keywords = anno_link(align_to = align_to, which = "row", panel_fun = panel_fun, 
		    	size = gbl_h, gap = unit(2, "mm"), width = gbl_w + unit(5, "mm"),
		    	link_gp = gpar(fill = "#DDDDDD", col = "#AAAAAA"), internal_line = FALSE))
		} else {
			if(any(cl == "0")) {
				ht = ht + Heatmap(ifelse(cl == "0", "< 5", ">= 5"), col = c("< 5" = "darkgreen", ">= 5" = "white"), width = unit(1, "mm"),
					heatmap_legend_param = list(title = "", at = "< 5", labels = "Small clusters"),
					show_column_names = FALSE)
			}
		}
	}
	draw(ht, gap = unit(2, "pt"), ...)
}

# == title
# Scale font size
#
# == param
# -x A numeric vector.
# -rg The range.
# -fs Range of th font size.
#
# == detaisl
# It is a linear interpolation.
#
# == value
# A numeric vector.
#
# == example
# scale_fontsize(runif(10, min = 1, max = 20))
scale_fontsize = function(x, rg = c(1, 30), fs = c(4, 16)) {
	k = (fs[2] - fs[1])/(rg[2] - rg[1]) 
	b = fs[2] - k*rg[2]
	y = k*x + b
	y[y < fs[1]] = fs[1]
	y[y > fs[2]] = fs[2]
	round(y)
}

