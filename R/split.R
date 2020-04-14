
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
