
# == title
# Cluster functional terms
#
# == param
# -mat A similarity matrix.
# -method Method for clustering the matrix.
# -catch_error Internally used.
# -verbose Whether to print messages.
# -... Other arguments passed to the clustering function.
#
# == details
# The following methods are supported:
#
# -``binary_cut`` see `binary_cut`.
# -``kmeans`` see `cluster_by_kmeans`.
# -``dynamicTreeCut`` see `cluster_by_dynamicTreeCut`.
# -``mclust`` see `cluster_by_mclust`.
# -``apcluster`` see `cluster_by_apcluster`.
# -``fast_greedy`` see `cluster_by_igraph`.
# -``leading_eigen`` see `cluster_by_igraph`.
# -``louvain`` see `cluster_by_igraph`.
# -``walktrap`` see `cluster_by_igraph`.
#
# Note the parametes for each clustering method are passes by ``...`` from `cluster_terms`.
#
# == value
# A numeric vector of cluster labels (in numeric).
#
# If ``catch_error`` is set to ``TRUE`` and if the clustering produces an error,
# the function returns a ``try-error`` object.
cluster_terms = function(mat, method = "binary_cut", catch_error = FALSE, verbose = TRUE, ...) {
	
	if(verbose) qqcat("cluster @{nrow(mat)} terms by @{method}...")

	if(any(method %in% c("cluster_by_kmeans", "kmeans"))) {
		oe = try(cl <- cluster_by_kmeans(mat, ...), silent = TRUE)
	} else if(any(method %in% c("cluster_by_dynamicTreeCut", "dynamicTreeCut"))) {
		oe  = try(cl <- cluster_by_dynamicTreeCut(mat, ...), silent = TRUE)
	} else if(any(method %in% c("cluster_by_mclust", "mclust"))) {
		oe = try(cl <- cluster_by_mclust(mat, ...), silent = TRUE)
	} else if(any(method %in% c("cluster_by_apcluster", "apcluster"))) {
		oe = try(cl <- cluster_by_apcluster(mat, ...), silent = TRUE)
	} else if(any(method %in% c("cluster_fast_greedy", "fast_greedy"))) {
		oe = try(cl <- cluster_by_igraph(mat, method = "cluster_fast_greedy", ...), silent = TRUE)
	} else if(any(method %in% c("cluster_leading_eigen", "leading_eigen"))) {
		oe = try(cl <- cluster_by_igraph(mat, method = "cluster_leading_eigen", ...), silent = TRUE)
	} else if(any(method %in% c("cluster_louvain", "louvain"))) {
		oe = try(cl <- cluster_by_igraph(mat, method = "cluster_louvain", ...), silent = TRUE)
	} else if(any(method %in% c("cluster_walktrap", "walktrap"))) {
		oe = try(cl <- cluster_by_igraph(mat, method = "cluster_walktrap", ...), silent = TRUE)
	} else if(any(method %in% c("cluster_by_binarycut", "binary_cut"))) {
		oe = try(cl <- binary_cut(mat, ...), silent = TRUE)
	} else {
		stop_wrap(qq("method '@{method}' is not supported."))
	}

	if(inherits(oe, "try-error")) {
		if(catch_error) {
			return(oe)
		} else {
			cat("\n")
			stop(oe)
		}
	}

	if(verbose) qqcat(" @{length(unique(cl))} clusters.\n")

	return(cl)
}

# == title
# Cluster similarity matrix by k-means clustering
#
# == param
# -mat The similarity matrix.
# -... Other arguments passed to `stats::kmeans`.
#
# == details
# The number of clusters are tried from 2 to ``min(round(nrow(mat)/5), 100)``. The best number
# of k for k-means clustering is identified according to the "elbow" or "knee" method on
# the distribution of within-cluster sum of squares at each k.
#
# == value
# A vector of cluster labels (in numeric).
#
cluster_by_kmeans = function(mat, ...) {
	
	cl = list()
	wss = NULL
	max_km = min(round(nrow(mat)/5), 100)

	for (i in 2:max_km) {
		km = kmeans(mat, centers = i, iter.max = 50)
		cl[[i - 1]] = km$cluster
		wss[i - 1] = sum(km$withinss)
	}
	best_km = min(elbow_finder(2:max_km, wss)[1], knee_finder(2:max_km, wss)[1])

	cl[[best_km - 1]]
}

# https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
elbow_finder = function(x_values, y_values) {
  # Max values to create line
  max_x_x = max(x_values)
  max_x_y = y_values[which.max(x_values)]
  max_y_y = max(y_values)
  max_y_x = x_values[which.max(y_values)]
  max_df = data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))

  # Creating straight line between the max values
  fit = lm(max_df$y ~ max_df$x)

  # Distance from point to line
  distances = c()
  for(i in seq_along(x_values)) {
    distances = c(distances, abs(coef(fit)[2]*x_values[i] - y_values[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2))
  }

  # Max distance point
  x_max_dist = x_values[which.max(distances)]
  y_max_dist = y_values[which.max(distances)]

  return(c(x_max_dist, y_max_dist))
}

# https://raghavan.usc.edu//papers/kneedle-simplex11.pdf
knee_finder = function(x, y) {
	n = length(x)
	a = (y[n] - y[1])/(x[n] - x[1])
	b = y[1] - a*x[1]
	d = a*x - y
	x[which.max(d)]
}

# == title
# Cluster similarity matrix by dynamicTreeCut
#
# == param
# -mat The similarity matrix.
# -minClusterSize Minimal number of objects in a cluster. Pass to `dynamicTreeCut::cutreeDynamic`.
# -... Other arguments passed to `dynamicTreeCut::cutreeDynamic`.
#
# == value
# A vector of cluster labels (in numeric).
#
cluster_by_dynamicTreeCut = function(mat, minClusterSize = 5, ...) {
	if(!requireNamespace("dynamicTreeCut", quietly = TRUE)) {
		stop_wrap("Package dynamicTreeCut should be installed.")
	}
	cl = dynamicTreeCut::cutreeDynamic(hclust(dist(mat)), distM = 1 - mat, minClusterSize = minClusterSize, verbose = 0, ...)
	as.character(unname(cl))
}

# == title
# Cluster similarity matrix by graph community detection methods
#
# == param
# -mat The similarity matrix.
# -method The community detection method.
# -... Other arguments passed to the corresponding community detection function, see Details.
#
# == details
# The symmetric similarity matrix can be treated as an adjacency matrix and constructed as a graph/network with the similarity values as the weight of hte edges.
# Thus, clustering the similarity matrix can be treated as detecting clusters/modules/communities from the graph.
#
# Four methods implemented in igraph package can be used here:
#
# -``cluster_fast_greedy`` uses `igraph::cluster_fast_greedy`.
# -``cluster_leading_eigen`` uses `igraph::cluster_leading_eigen`.
# -``cluster_louvain`` uses `igraph::cluster_louvain`.
# -``cluster_walktrap`` uses `igraph::cluster_walktrap`.
#
# == value
# A vector of cluster labels (in numeric).
#
cluster_by_igraph = function(mat, 
    method = c("cluster_fast_greedy",
               "cluster_leading_eigen",
               "cluster_louvain",
               "cluster_walktrap"),
    ...) {

	if(!requireNamespace("igraph", quietly = TRUE)) {
		stop_wrap("Package igraph should be installed.")
	}

	if(is.character(method)) {
		method = match.arg(method)[1]
		method = getFromNamespace(method, ns = "igraph")
	}
	g = igraph::graph_from_adjacency_matrix(mat, mode = "upper", weighted = TRUE)

	imc = method(g, ...)
	as.vector(igraph::membership(imc))
}

# == title
# Cluster similarity matrix by mclust
#
# == param
# -mat The similarity matrix.
# -... Other arguments passed to `mclust::Mclust`.
#
# == details
# The value of ``G`` in `mclust::Mclust` is set to ``1:min(round(nrow(mat)/5), 100)``.
#
# == value
# A vector of cluster labels (in numeric).
#
cluster_by_mclust = function(mat, ...) {
	if(!requireNamespace("mclust", quietly = TRUE)) {
		stop_wrap("Package mclust should be installed.")
	}
	mclustBIC = mclust::mclustBIC
	fit = mclust::Mclust(as.matrix(mat), G = 1:min(round(nrow(mat)/5), 100), verbose = FALSE, ...)
	unname(fit$classification)
}

# == title
# Cluster similarity matrix by apcluster
#
# == param
# -mat The similarity matrix.
# -... Other arguments passed to `apcluster::apcluster`.
#
# == value
# A vector of cluster labels (in numeric).
#
cluster_by_apcluster = function(mat, ...) {
	if(!requireNamespace("apcluster", quietly = TRUE)) {
		stop_wrap("Package apcluster should be installed.")
	}
	x = apcluster::apcluster(apcluster::negDistMat(r = 2), mat, ...)
	cl = numeric(nrow(mat))
	for(i in seq_along(x@clusters)) {
		cl[x@clusters[[i]]] = i
	}
	cl
}

