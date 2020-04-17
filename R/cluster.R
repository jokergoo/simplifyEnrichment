
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
			l = sapply(d, is.leaf)
			if(all(l)) {
				attr(d, "score2") = s
			} else {
				attr(d, "score2") = max(s, sapply(d[!l], function(x) attr(x, "score")))
			}
		}

		return(d)
	})

	return(dend)
}

.cluster_mat = function(mat, dist_mat = dist(mat), .env, index = seq_len(nrow(mat)), 
	depth = 0, dend_index = NULL) {

	if(nrow(mat) == 1) {
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

	# cl = cutree(hclust(dist_mat), k = 2)
	if(nrow(mat) == 2) {
		cl = c(1, 2)
	} else {
		cl = kmeans(mat, centers = 2)$cluster
	}
	l1 = cl == 1
	l2 = cl == 2
	x11 = .env$value_fun(mat[l1, l1])
	x12 = .env$value_fun(mat[l1, l2])
	x21 = .env$value_fun(mat[l2, l1])
	x22 = .env$value_fun(mat[l2, l2])

	s = (x11 + x22)/(x11 + x12 + x21 + x22)
	if(is.na(s)) s = 1

	sr = numeric(1)
	for(i in 1:10) {
		clr = sample(cl, length(cl))
		l1r = clr == 1
		l2r = clr == 2
		x11r = sample(as.vector(mat[l1r, l1r]), 1)
		x12r = sample(as.vector(mat[l1r, l2r]), 1)
		x21r = sample(as.vector(mat[l2r, l1r]), 1)
		x22r = sample(as.vector(mat[l2r, l2r]), 1)

		sr[i] = (x11r + x22r)/(x11r + x12r + x21r + x22r)
		if(is.na(sr[i])) sr[i] = 1
	}
	sr = mean(sr)

	if(is.null(dend_index)) {
		.env$dend = list()
		attributes(.env$dend) = list(
			members = nrow(mat),
			height = depth,
			score = s,
			random_score = sr,
			index = NULL,
			class = "dendrogram"
		)
	} else {
		.env$dend[[dend_index]] = list()
		attributes(.env$dend[[dend_index]]) = list(
			members = nrow(mat),
			height = depth,
			score = s,
			random_score = sr,
			index = dend_index,
			class = "dendrogram"
		)
	}

	.cluster_mat(mat[l1, l1, drop = FALSE], as.dist(as.matrix(dist_mat)[l1, l1, drop = FALSE]), 
		.env, index = index[l1], depth = depth + 1, dend_index = c(dend_index, 1))
	.cluster_mat(mat[l2, l2, drop = FALSE], as.dist(as.matrix(dist_mat)[l2, l2, drop = FALSE]), 
		.env, index = index[l2], depth = depth + 1, dend_index = c(dend_index, 2))
}

render_dend = function(dend, field = "score2", cutoff = 0.8) {
	col_fun = colorRamp2(c(0.5, 0.75, 1), c("blue", "yellow", "red"))
	dend = dendrapply(dend, function(d) {
		if(!is.leaf(d)) {
			s = attr(d, field)
			attr(d, "nodePar") = list(pch = ifelse(s > cutoff, 16, 4), cex = 0.5, col = col_fun(s))
		}
		return(d)
	})
	
	plot(dend, main = field)
	box()
}

