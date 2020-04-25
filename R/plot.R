
# == title
# Simplify GO enrichment results
#
# == param
# -mat A GO similarity matrix.
# -cl Cluster labels.
# -dend Used internally.
# -draw_word_cloud Whether draw the word clouds.
# -min_term Minimal number of GO terms in a cluster.
# -order_by_size Whether to reorder GO clusters by their sizes.
# -exclude_words Words that are excluded in the word cloud.
# -max_words Maximal number of words in the word cloud.
# -word_cloud_width The maximal width of the viewport to put the word cloud. 
#            The value should be numeric. It is measured in mm.
# -... Other arguments
#
ht_GO_clusters = function(mat, cl, dend = NULL, 
	draw_word_cloud = TRUE, min_term = 5, order_by_size = FALSE,
	exclude_words = character(0), max_words = 10, word_cloud_width = 60, ...) {

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
			keywords = tapply(rownames(mat), cl, function(go_id) {
				suppressMessages(suppressWarnings(df <- count_word(go_id, exclude_words = exclude_words)))
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

			gbl = lapply(names(align_to), function(nm) {
				kw = rev(keywords[[nm]][, 1])
				freq = rev(keywords[[nm]][, 2])
				fontsize = scale_fontsize(freq, rg = c(1, max(10, freq)))
				simple_word_cloud_grob(kw, fontsize, max_width = word_cloud_width)
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
			    grid.draw(gbl[[nm]])
			    popViewport()
			    popViewport()
			}

			ht = ht + rowAnnotation(keywords = anno_link(align_to = align_to, which = "row", panel_fun = panel_fun, 
		    	size = gbl_h, gap = unit(2, "mm"), width = gbl_w + unit(5, "mm"),
		    	link_gp = gpar(fill = "#DDDDDD", col = "#AAAAAA"), internal_line = FALSE))
		} else {
			if(any(cl == "0")) {
				ht = ht + Heatmap(ifelse(cl == "0", "< 5", ">= 5"), col = c("< 5" = "darkgreen", ">= 5" = "white"), width = unit(2, "mm"),
					heatmap_legend_param = list(title = "", at = "< 5", lable = "Clusters with\nsize < 5"),
					show_column_names = FALSE)
			}
		}
	}
	draw(ht, gap = unit(2, "pt"), ...)
}


scale_fontsize = function(x, rg = c(1, 30), fs = c(4, 16)) {
	k = (fs[2] - fs[1])/(rg[2] - rg[1]) 
	b = fs[2] - k*rg[2]
	y = k*x + b
	y[y < fs[1]] = fs[1]
	y[y > fs[2]] = fs[2]
	y
}

