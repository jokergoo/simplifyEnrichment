

##### old implementation
cluster_mat_old = function(mat, value_fun = area_above_ecdf, partition_fun = partition_by_pam,
	cutoff = 0.85) {

	env = new.env()
	env$value_fun = value_fun
	env$partition_fun = partition_fun

	.cluster_mat_old(mat, dist_mat = dist(mat), .env = env)
	dend = env$dend
	dend = rev(dend)

	max_d = simplifyEnrichment:::dend_max_depth(dend)
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

	## bottom-up
	assign_score2 = function(d) {
		if(is.leaf(d)) {
			attr(d, "score2") = 0.5
			return(d)
		}
		# check the two child node
		if(is.null(attr(d[[1]], "score2"))) {
			d[[1]] = assign_score2(d[[1]])
		}
		if(is.null(attr(d[[2]], "score2"))) {
			d[[2]] = assign_score2(d[[2]])
		}

		s2 = max(attr(d[[1]], "score2"), attr(d[[2]], "score2"))
		# the node is assigned with its child nodes' score only if
		# 1. the size of the two children should not be too extremely different
		# 2. the score should not be too smaller than its child nodes' 
		s = attr(d, "score")
		if( min(nobs(d[[1]]), nobs(d[[2]]))/nobs(d) < 0.1 ) {
			attr(d, "score2")  = s
		} else if(s > cutoff*0.8 && s < s2) {
			attr(d, "score2") = s2
		} else {
			attr(d, "score2")  = s
		}
		return(d)
	}
	dend = assign_score2(dend)

	hash = digest::digest(list(mat, value_fun, partition_fun, cutoff))
	# dend_env$dend = dend
	# dend_env$hash = hash

	return(dend)
}

.cluster_mat_old = function(mat, dist_mat = dist(mat), .env, index = seq_len(nrow(mat)), 
	depth = 0, dend_index = NULL) {
	nr = nrow(mat)

	if(nr == 1) {
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

	if(nrow(mat) == 2) {
		cl = c(1, 2)
	} else {
		oe = try(suppressWarnings(cl <- .env$partition_fun(mat)), silent = TRUE)
		if(inherits(oe, "try-error")) {
			cl = rep(2, nr)
			cl[seq_len(ceiling(nr/2))] = 1
		}
	}

	l1 = cl == 1
	l2 = cl == 2
	if(mean(mat[l1, l1]) > mean(mat[l2, l2])) {
		l3 = l1
		l1 = l2
		l2 = l3
	}
	m11 = mat[l1, l1, drop = FALSE]
	if(nrow(m11) == 1) {
		x11 = 1
		x12 = .env$value_fun(mat[l1, l2])
		x21 = .env$value_fun(mat[l2, l1])
	} else {
		m11 = m11[ lower.tri(m11) | upper.tri(m11)]
		x11 = .env$value_fun(m11)
		x12 = .env$value_fun(mat[l1, l2])
		x21 = .env$value_fun(mat[l2, l1])
	}
	m22 = mat[l2, l2, drop = FALSE]
	if(nrow(m22) == 1) {
		x22 = 1
	} else {
		m22 = m22[ lower.tri(m22) | upper.tri(m22)]
		x22 = .env$value_fun(m22)
	}
	
	s = (x11 + x22)/(x11 + x12 + x21 + x22)
	if(is.na(s)) s = 1

	if(is.null(dend_index)) {
		.env$dend = list()
		attributes(.env$dend) = list(
			members = nrow(mat),
			height = depth,
			score = s,
			index = NULL,
			class = "dendrogram"
		)
	} else {
		.env$dend[[dend_index]] = list()
		attributes(.env$dend[[dend_index]]) = list(
			members = nrow(mat),
			height = depth,
			score = s,
			index = dend_index,
			class = "dendrogram"
		)
	}

	.cluster_mat_old(mat[l1, l1, drop = FALSE], as.dist(as.matrix(dist_mat)[l1, l1, drop = FALSE]), 
		.env, index = index[l1], depth = depth + 1, dend_index = c(dend_index, 1))
	.cluster_mat_old(mat[l2, l2, drop = FALSE], as.dist(as.matrix(dist_mat)[l2, l2, drop = FALSE]), 
		.env, index = index[l2], depth = depth + 1, dend_index = c(dend_index, 2))
}

### test

mat = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds", package = "simplifyEnrichment"))

dend1 = cluster_mat_old(mat)

dend1 = dendrapply(dend1, function(d) {
	attr(d, "index") = NULL
	d
})
dend2 = simplifyEnrichment:::cluster_mat(mat)


context("test cluster_mat")

test_that("test cluster_mat", {
	expect_equal(attributes(dend1), attributes(dend2))
	expect_equal(attributes(dend1[[1]]), attributes(dend2[[1]]))
	expect_equal(attributes(dend1[[c(1, 1)]]), attributes(dend2[[c(1, 1)]]))

	expect_equal(dend1, dend2)
})
