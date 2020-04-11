
plot_heatmap = function(mat, cl, min_term = 5, order_by_size = TRUE,
	exclude_words = character(0)) {

	cl = as.character(cl)
	cl_tb = table(cl)
	cl[cl %in% names(cl_tb[cl_tb < min_term])] = "0"
	cl = factor(cl, levels = c(setdiff(names(sort(table(cl), decreasing = TRUE)), "0"), "0"))

	keywords = tapply(rownames(mat), cl, function(go_id) {
		suppressMessages(suppressWarnings(df <- count_word(go_id, exclude_words = exclude_words)))
		df = df[df$freq > 1, , drop = FALSE]
		if(nrow(df) > 10) {
			df = df[order(df$freq, decreasing = TRUE)[1:10], ]
		}
		df
	})

	ht = Heatmap(mat, col = colorRamp2(c(0, 1), c("white", "red")),
		name = "Similarity", column_title = "GO Similarity",
		show_row_names = FALSE, show_column_names = FALSE, 
		show_row_dend = FALSE, show_column_dend = FALSE,
		cluster_row_slices = !order_by_size, cluster_column_slices = !order_by_size,
		row_split = cl, column_split = cl, 
		border = "#404040", row_title = NULL,
		row_gap = unit(0, "mm"), column_gap = unit(0, "mm"))

	align_to = split(1:nrow(mat), cl)
	align_to = align_to[names(align_to) != "0"]

	gbl = lapply(names(align_to), function(nm) {
		kw = rev(keywords[[nm]][, 1])
		freq = rev(keywords[[nm]][, 2])
		fontsize = scale_fontsize(freq, rg = c(1, max(10, freq)))
		simple_word_cloud_grob(kw, fontsize, max_width = 60)
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

	ht = ht + rowAnnotation(keywords = anno_zoom(align_to = align_to, which = "row", panel_fun = panel_fun, 
    	size = gbl_h, gap = unit(2, "mm"), width = gbl_w + unit(5, "mm"),
    	link_gp = gpar(fill = "#DDDDDD", col = "#AAAAAA"), internal_line = FALSE))
	draw(ht, gap = unit(2, "pt"))
}


scale_fontsize = function(x, rg = c(1, 30), fs = c(4, 16)) {
	k = (fs[2] - fs[1])/(rg[2] - rg[1]) 
	b = fs[2] - k*rg[2]
	y = k*x + b
	y[y < fs[1]] = fs[1]
	y[y > fs[2]] = fs[2]
	y
}

