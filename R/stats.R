
# == title
# Difference score
#
# == param
# -mat The similarity matrix.
# -cl Cluster labels.
#
# == details
# This function measures the different between the similarity values for the terms
# that belong to the same clusters and in different clusters. The difference score
# is the Kolmogorov-Smirnov statistic between the two distributions.
#
# == value
# A numeric scalar.
#
# == examples
# mat = readRDS(system.file("extdata", "similarity_mat.rds", package = "simplifyEnrichment"))
# cl = binary_cut(mat)
# difference_score(mat, cl)
difference_score = function(mat, cl) {
	n = nrow(mat)
	l_block = matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat))

	for(le in unique(cl)) {
		l = cl == le
		l_block[l, l] = TRUE
	}
	l_block2 = l_block
	l_block2[upper.tri(mat)] = FALSE
	x1 = mat[l_block2]

	l_block2 = !l_block
	l_block2[upper.tri(mat)] = FALSE
	x2 = mat[l_block2]

	ecdf1 = ecdf(x1)
	ecdf2 = ecdf(x2)

	p = seq(0, 1, length = 1000)
	max(abs(ecdf2(p) - ecdf1(p)))
}

block_mean = function(mat, cl) {
	n = nrow(mat)
	l_block = matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat))

	for(le in unique(cl)) {
		l = cl == le
		l_block[l, l] = TRUE
	}
	l_block2 = l_block
	l_block2[upper.tri(mat)] = FALSE
	x1 = mat[l_block2]
	mean(x1)
}

other_mean = function(mat, cl) {
	n = nrow(mat)
	l_block = matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat))

	for(le in unique(cl)) {
		l = cl == le
		l_block[l, l] = TRUE
	}

	l_block2 = !l_block
	l_block2[upper.tri(mat)] = FALSE
	x2 = mat[l_block2]
	mean(x2)
}


 # sapply(clt,function(cl)modularity(g, membership = cl, weights = get.edge.attribute(g, "weight")))
