
library(GOSemSim)
	
simplify_GO = function(go_id) {
	hsGO = godata('org.Hs.eg.db', ont = "BP")
	go_sim = mgoSim(all_go_id, all_go_id, measure = "Rel", semData = hsGO, combine = NULL)
	go_sim[is.na(go_sim)] = 0
	diag(go_sim) = 1
}

.cluster_mat = function(mat, dist_mat = dist(mat), .env, index = seq_len(nrow(mat)), 
	depth = 0, dend_index = NULL) {

	if(nrow(mat) == 1) {
		.env$dend[[dend_index]] = index
		attributes(.env$dend[[dend_index]]) = list(
			members = 1,
			label = index,
			leaf = TRUE,
			height = depth,
			score = 1,
			class = "dendrogram"
		)
		return(NULL)
	}

	cl = cutree(hclust(dist_mat), k = 2)
	l1 = cl == 1
	l2 = cl == 2
	x11 = median(mat[l1, l1])
	x12 = median(mat[l1, l2])
	x21 = median(mat[l2, l1])
	x22 = median(mat[l2, l2])

	s = (x11 + x22)/(x11 + x12 + x21 + x22)
	if(is.na(s)) s = 1

	if(is.null(dend_index)) {
		.env$dend = list()
		attributes(.env$dend) = list(
			members = nrow(mat),
			height = depth,
			score = s,
			class = "dendrogram"
		)
	} else {
		.env$dend[[dend_index]] = list()
		attributes(.env$dend[[dend_index]]) = list(
			members = nrow(mat),
			height = depth,
			score = s,
			class = "dendrogram"
		)
	}

	.cluster_mat(mat[l1, l1, drop = FALSE], as.dist(as.matrix(dist_mat)[l1, l1, drop = FALSE]), 
		.env, index = index[l1], depth = depth + 1, dend_index = c(dend_index, 1))
	.cluster_mat(mat[l2, l2, drop = FALSE], as.dist(as.matrix(dist_mat)[l2, l2, drop = FALSE]), 
		.env, index = index[l2], depth = depth + 1, dend_index = c(dend_index, 2))
}

cluster_mat = function(mat) {

	env = new.env()
	.cluster_mat(mat, dist_mat = dist(mat), .env = env)

	dend = env$dend
	max_d = dendextend::max_depth(dend)
	dend = dendrapply(dend, function(d) {
		attr(d, "height") = max_d - attr(d, "height")
		d
	})

	dend = ComplexHeatmap::adjust_dend_by_x(dend)
}


make_rule = function(dend, plot = FALSE) {
	size = get_nodes_attr(dend, "member")
	tb = table(size)

	score = get_nodes_attr(dend, "score")
	plot(size, score, log = "x")

	find_midpoint = function(x) {
		x = sort(x)
		i = which.max(diff(x))
		(x[i + 1] + x[i])/2
	}

	cluster = rep(NA, length(score))
	i = 1
	n = length(tb)
	while(i <= length(tb)) {
		if(tb[i] >= 50) {
			l = size == as.numeric(names(tb)[i])
			midpoint = find_midpoint(score[l])
			cl = rep(1, sum(l))
			cl[score[l] > midpoint] = 2
			cluster[l] = cl
		} else {
			sum_size = tb[i]
			i_size = i

			if(i < n) {
				if(sum(tb[seq(i+1, n)]) < 50) {
					i_size = c(i_size, seq(i+1, n))
					i = n+1
				} else {
					while(i < n) {
						if(sum_size + tb[i + 1] < 100) {
							sum_size = sum_size + tb[i+1]
							i_size = c(i_size, i+1)
							i = i + 1
						} else if(sum(tb[seq(i+1, n)]) < 50) {
							i_size = c(i_size, seq(i+1, n))
							i = n+1
						} else {
							break
						}
					}
				}
			}

			l = size %in% as.numeric(names(tb)[i_size])
			midpoint = find_midpoint(score[l])
			cl = rep(1, sum(l))
			cl[score[l] > midpoint] = 2
			cluster[l] = cl
		}
		i = i + 1
	}

	cluster[size == 1] = 2

	if(plot) plot(size, score, col = cluster, log = "x")

	tapply(seq_along(size), size, function(index)  {
		cl = cluster[index]
		s = score[index]
		f = local({
			if(all(cl == 1)) {
				function(x) 1
			} else if(all(cl == 2)) {
				function(x) 2
			} else {
				x1 = max(s[cl == 1])
				x2 = min(s[cl == 2])
				mid = (x1+x2)/2
				function(x) {
					ifelse(x > mid, 2, 1)
				}
			}
		})
		return(f)
	})
}


split_dend = function(dend) {
	
	rule = make_rule(dend)

	dend2 = dendrapply(dend, function(d) {
		size = attr(d, "member")
		score = attr(d, "score")
		if(rule[[ as.character(size) ]](score) == 1) {
			attr(d, "height") = 0
		}
		d
	})

	as.character(cutree(dend2, h = 0.1))
}

plot_dend = function(dend, mat) {

	cl = as.character(split_dend(dend))
	cl_tb = table(cl)
	cl[cl %in% names(cl_tb[cl_tb <= 5])] = "0"
	cl = factor(cl, levels = c(setdiff(names(sort(table(cl), decreasing = TRUE)), "0"), "0"))

	ht = Heatmap(mat, col = colorRamp2(c(0, 1), c("white", "red")),
		show_row_names = FALSE, show_column_names = FALSE, 
		show_row_dend = FALSE, show_column_dend = FALSE,
		cluster_row_slices = FALSE, cluster_column_slices = FALSE,
		row_split = cl, column_split = cl, 
		border = "black", row_title = NULL, column_title = NULL,
		row_gap = unit(0, "mm"), column_gap = unit(0, "mm"))
	draw(ht)
}


dend = cluster_mat(go_sim)


