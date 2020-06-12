
term_similarity = function(gl) {
	all = unique(unlist(gl))
	gl = lapply(gl, function(x) as.numeric(factor(x, levels = all)))
	n = length(gl)

	mg = matrix(0, ncol = length(all), nrow = n)
	for(i in seq_len(n)) {
		mg[i, gl[[i]]] = 1
	}

	mat = 1 - dist(mg, method = "binary")
	mat = as.matrix(mat)
	rownames(mat) = colnames(mat) = names(gl)
	mat
}
