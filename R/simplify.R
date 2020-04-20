

# == title
# Simplify GO enrichment results
#
# == param
# -mat A GO similarity matrix.
# -method method
# -... other arguments
# -plot Whether to make the heatmap.
# -min_term Minimal number of GO terms in a cluster.
# -order_by_size Whether to reorder GO clusters by their sizes.
# -exclude_words Words that are excluded in the word cloud.
# -max_words Maximal number of words in the word cloud.
#
# == example
# \dontrun{
# set.seed(123)
# go_id = simplifyGO:::random_GO(500)
# mat = GO_similarity(go_id)
# lt = simplifyGO(mat)
# }
simplifyGO = function(mat, method = "binary_cut", ..., 

	# parameters for plotting the heatmap
	plot = TRUE, min_term = 5, order_by_size = TRUE, 
	exclude_words = character(0), max_words = 10) {
	
	cl = cluster_GO(mat, method, ...)

	if(plot) ht_GO_clusters(mat, cl, min_term = min_term, order_by_size = order_by_size, 
		exclude_words = exclude_words, max_words = max_words, column_title = qq("GO terms are clustered by '@{method}'"))

	go_id = rownames(mat)
	term = select(GO.db::GO.db, keys = go_id, columns = "TERM")$TERM
	
	return(data.frame(id = go_id, name = term, cluster = cl, stringsAsFactors = FALSE))
}

cluster_GO = function(mat, method = "binary_cut", ...) {
	
	qqcat("cluster GO terms by @{method}...")

	if(any(method %in% c("cluster_by_kmeans", "kmeans"))) {
		cl = cluster_by_kmeans(mat, ...)
	} else if(any(method %in% c("cluster_by_dynamicTreeCut", "dynamicTreeCut"))) {
		cl = cluster_by_dynamicTreeCut(mat, ...)
	} else if(any(method %in% c("cluster_by_mclust", "mclust"))) {
		cl = cluster_by_mclust(mat, ...)
	} else if(any(method %in% c("cluster_by_apclust", "apclust"))) {
		cl = cluster_by_apclust(mat, ...)
	} else if(any(method %in% c("cluster_fast_greedy", "fast_greedy"))) {
		cl = cluster_by_igraph(mat, method = "cluster_fast_greedy", ...)
	} else if(any(method %in% c("cluster_leading_eigen", "leading_eigen"))) {
		cl = cluster_by_igraph(mat, method = "cluster_leading_eigen", ...)
	} else if(any(method %in% c("cluster_louvain", "louvain"))) {
		cl = cluster_by_igraph(mat, method = "cluster_louvain", ...)
	} else if(any(method %in% c("cluster_walktrap", "walktrap"))) {
		cl = cluster_by_igraph(mat, method = "cluster_walktrap", ...)
	} else if(any(method %in% c("cluster_by_binarycut", "binary_cut"))) {
		cl = binary_cut(mat, ...)
	} else {
		stop_wrap(qq("method '@{method}' is not supported."))
	}

	qqcat(" @{length(unique(cl))} clusters.\n")

	return(cl)
}

# == title
# Cluster GO similarity matrix by k-means clustering
#
# == param
# -mat The Go similarity matrix.
# -... Other arguments passed to `stats::kmeans`.
#
cluster_by_kmeans = function(mat, ...) {
	if(!requireNamespace("cluster")) {
		stop_wrap("Package cluster should be installed.")
	}
	cl = list()
	for(k in 2:min(round(nrow(mat)/5), 100)) {
		cl[[k-1]] = kmeans(mat, k, ...)$cluster
	}
	i = which.max(sapply(cl, function(x) {
		n = length(x)
		if(n <= 500) {
			mean(cluster::silhouette(x, dist(mat))[, "sil_width"])
		} else {
			ind = sort(sample(n, 500))
			mean(cluster::silhouette(x[ind], dist(mat[ind, ind]))[, "sil_width"])
		}
	}))
	cl[[i]]
}

# == title
# Cluster GO similarity matrix by dynamicTreeCut
#
# == param
# -mat The Go similarity matrix.
# -... Other arguments passed to `dynamicTreeCut::cutreeDynamic`.
#
cluster_by_dynamicTreeCut = function(mat, ...) {
	if(!requireNamespace("dynamicTreeCut")) {
		stop_wrap("Package dynamicTreeCut should be installed.")
	}
	cl = dynamicTreeCut::cutreeDynamic(hclust(dist(mat)), distM = 1 - mat, minClusterSize = 5, verbose = 0, ...)
	as.character(unname(cl))
}

# == title
# Cluster GO similarity matrix by community detection methods
#
# == param
# -mat The Go similarity matrix.
# -method method
# -... Other arguments passed to.
#
cluster_by_igraph = function(mat, 
    method = c("cluster_fast_greedy",
               "cluster_leading_eigen",
               "cluster_louvain",
               "cluster_walktrap"),
    ...) {

	if(!requireNamespace("igraph")) {
		stop_wrap("Package igraph should be installed.")
	}

	if(is.character(method)) {
		method = match.arg(method)[1]
		method = getFromNamespace(method, ns = "igraph")
	}
	g = igraph::graph_from_adjacency_matrix(mat, mode = "upper", weight = TRUE)

	imc = method(g, ...)
	as.vector(igraph::membership(imc))
}

# == title
# Cluster GO similarity matrix by mclust
#
# == param
# -mat The Go similarity matrix.
# -... Other arguments passed to `mclust::Mclust`.
#
cluster_by_mclust = function(mat, ...) {
	if(!requireNamespace("mclust")) {
		stop_wrap("Package mclust should be installed.")
	}
	fit = mclust::Mclust(as.matrix(mat), G = 1:min(round(nrow(mat)/5), 100), verbose = TRUE, ...)
	unname(fit$classification)
}

# == title
# Cluster GO similarity matrix by apcluster
#
# == param
# -mat The Go similarity matrix.
# -... Other arguments passed to `apcluster::apcluster`.
#
cluster_by_apclust = function(mat, ...) {
	if(!requireNamespace("apcluster")) {
		stop_wrap("Package apcluster should be installed.")
	}
	x = apcluster::apcluster(apcluster::negDistMat(r = 2), mat, ...)
	cl = numeric(nrow(mat))
	for(i in seq_along(x@clusters)) {
		cl[x@clusters[[i]]] = i
	}
	cl
}


