
# make_rule = function(dend, cutoff = 0.8) {
# 	size = dendextend::get_nodes_attr(dend, "member")
	
# 	score = dendextend::get_nodes_attr(dend, "score2")

# 	tapply(seq_along(size), size, function(index)  {
# 		s = score[index]
# 		f = local({
# 			cutoff = cutoff
# 			function(x) {
# 				ifelse(x > cutoff, 2, 1)
# 			}
			
# 		})
# 		return(f)
# 	})
# }

# split_dend = function(dend, cutoff = 0.8) {
	
# 	rule = make_rule(dend, cutoff = cutoff)

# 	dend2 = dendrapply(dend, function(d) {
# 		size = attr(d, "member")
# 		score = attr(d, "score")
# 		if(rule[[ as.character(size) ]](score) == 1) {
# 			attr(d, "height") = 0
# 		}
# 		d
# 	})

# 	as.character(cutree(dend2, h = 0.1))
# }

cut_dend = function(dend, cutoff = 0.8, field = "score2", plot = FALSE) {

	if(attr(dend, field) < cutoff) {
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
		s = attr(d, field)
		if(s < cutoff) {
			attr(d, "height") = 0
		}
		d
	})

	## make sure all sub-nodes having height 0 if the node is 0 height
	is_parent_zero_height = function(index) {
		h = sapply(seq_along(index), function(i) {
			attr(dend2[[ index[1:i] ]], "height")
		})
		any(h == 0)
	}
	dend2 = edit_node(dend2, function(d, index) {
		if(is_parent_zero_height(index)) {
			attr(d, "height") = 0
			attr(d, "nodePar") = NULL
		}
		d
	})

	if(plot) {
		col_fun = colorRamp2(c(0.5, 0.75, 1), c("blue", "yellow", "red"))
		dend = edit_node(dend, function(d, index) {
			if(is.null(index)) {
				if(!is.leaf(d)) {
					if(attr(dend2, "height") > 0) {
						s = attr(d, field)
						attr(d, "nodePar") = list(pch = ifelse(s > cutoff, 16, 4), cex = 0.5, col = col_fun(s))
					}
				}
			} else {
				if(!is.leaf(d)) {
					if(attr(dend2[[index]], "height") > 0) {
						s = attr(d, field)
						attr(d, "nodePar") = list(pch = ifelse(s > cutoff, 16, 4), cex = 0.5, col = col_fun(s))
					}
				}
			}
			return(d)
		})
		
		plot(dend)
		box()
	}

	as.character(cutree(dend2, h = 0.1))
}

plot_heatmap = function(mat, cl, min = 5) {

	cl = as.character(cl)
	cl_tb = table(cl)
	cl[cl %in% names(cl_tb[cl_tb < min])] = "0"
	cl = factor(cl, levels = c(setdiff(names(sort(table(cl), decreasing = TRUE)), "0"), "0"))

	keywords = tapply(rownames(mat), cl, function(go_id) {
		df = count_word(go_id)
		df = df[df$freq > 1, , drop = FALSE]
		if(nrow(df) > 10) {
			df = df[order(df$freq, decreasing = TRUE)[1:10], ]
		}
		df
	})

	ht = Heatmap(mat, col = colorRamp2(c(0, 1), c("white", "red")),
		name = "Similarity",
		show_row_names = FALSE, show_column_names = FALSE, 
		show_row_dend = FALSE, show_column_dend = FALSE,
		# cluster_row_slices = FALSE, cluster_column_slices = FALSE,
		row_split = cl, column_split = cl, 
		border = "black", row_title = NULL, column_title = NULL,
		row_gap = unit(0, "mm"), column_gap = unit(0, "mm"))

	align_to = split(1:nrow(mat), cl)
	align_to = align_to[names(align_to) != "0"]

	gbl = lapply(names(align_to), function(nm) {
		kw = rev(keywords[[nm]][, 1])
		freq = rev(keywords[[nm]][, 2])
		fontsize = scale_fontsize(freq, rg = c(1, max(10, freq)))
		txt = qq("<span style='font-size:@{fontsize}px;line-height:10px;'>@{kw}</span> ")
		textbox_grob(txt, width = unit(8, "cm"))
	})
	names(gbl) = names(align_to)

	gbl_h = lapply(gbl, function(x) convertHeight(grobHeight(x), "cm") + unit(2, "mm"))
	gbl_h = do.call(unit.c, gbl_h)

	gbl_w = lapply(gbl, function(x) convertWidth(grobWidth(x), "cm"))
	gbl_w = do.call(unit.c, gbl_w)
	gbl_w = max(gbl_w)

	panel_fun = function(index, nm) {
		pushViewport(viewport())
	    grid.rect()
	    grid.draw(gbl[[nm]])
	    popViewport()
	}

	ht = ht + rowAnnotation(keywords = anno_zoom(align_to = align_to, which = "row", panel_fun = panel_fun, 
    	size = gbl_h, gap = unit(2, "mm"), width = gbl_w + unit(1, "cm")))
	draw(ht)
}


scale_fontsize = function(x, rg = c(1, 30), fs = c(4, 20)) {
	k = (fs[2] - fs[1])/(rg[2] - rg[1]) 
	b = fs[2] - k*rg[2]
	y = k*x + b
	y[y < fs[1]] = fs[1]
	y[y > fs[2]] = fs[2]
	round(y)
}
