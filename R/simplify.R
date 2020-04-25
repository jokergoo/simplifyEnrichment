

# == title
# Simplify GO enrichment results
#
# == param
# -mat A GO similarity matrix.
# -method Method for clustering the matrix.
# -... Other arguments passed to the clustering function, see details.
# -plot Whether to make the heatmap.
# -draw_word_cloud Whether draw the word clouds.
# -min_term Minimal number of GO terms in a cluster.
# -order_by_size Whether to reorder GO clusters by their sizes.
# -exclude_words Words that are excluded in the word cloud.
# -max_words Maximal number of words in the word cloud.
# -word_cloud_width The maximal width of the viewport to put the word cloud. 
#            The value should be numeric. It is measured in mm.
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
# Note the parametes for each clustering method are passes by ``...`` from `simplifyGO`.
#
# == value
# A data frame with three columns:
#
# - GO ID
# - GO term name
# - Cluster label
#
# == example
# \dontrun{
# set.seed(123)
# go_id = simplifyGO:::random_GO(500)
# mat = GO_similarity(go_id)
# df = simplifyGO(mat)
# }
simplifyGO = function(mat, method = "binary_cut", ..., 

	# parameters for plotting the heatmap
	plot = TRUE, draw_word_cloud = TRUE, min_term = 5, order_by_size = FALSE, 
	exclude_words = character(0), max_words = 10,
	word_cloud_width = 60) {
	
	cl = cluster_GO(mat, method, ...)

	if(plot) ht_GO_clusters(mat, cl, draw_word_cloud = draw_word_cloud, min_term = min_term, order_by_size = order_by_size, 
		exclude_words = exclude_words, max_words = max_words, word_cloud_width = word_cloud_width,
		column_title = qq("@{nrow(mat)} GO terms are clustered by '@{method}'"))

	go_id = rownames(mat)
	term = select(GO.db::GO.db, keys = go_id, columns = "TERM")$TERM
	
	return(data.frame(id = go_id, name = term, cluster = cl, stringsAsFactors = FALSE))
}

# == title
# Cluster GO terms
#
# == param
# -mat A GO similarity matrix.
# -method Method for clustering the matrix.
# -catch_error Internally used.
# -verbose Whether print messages.
# -... Other arguments passed to the clustering function.
#
# == value
# A numeric vector of cluster labels.
#
# If ``catch_error`` is set to ``TRUE`` and if the clustering produces an error,
# the function returns a ``try-error`` object.
cluster_GO = function(mat, method = "binary_cut", catch_error = FALSE, verbose = TRUE, ...) {
	
	if(verbose) qqcat("cluster @{nrow(mat)} GO terms by @{method}...")

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
# Cluster GO similarity matrix by k-means clustering
#
# == param
# -mat The GO similarity matrix.
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
elbow_finder <- function(x_values, y_values) {
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
  for(i in 1:length(x_values)) {
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
# Cluster GO similarity matrix by dynamicTreeCut
#
# == param
# -mat The GO similarity matrix.
# -minClusterSize Minimal number of elements in a cluster. Pass to `dynamicTreeCut::cutreeDynamic`.
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
# Cluster GO similarity matrix by graph community detection methods
#
# == param
# -mat The GO similarity matrix.
# -method The community detection method.
# -... Other arguments passed to the corresponding community detection function, see details.
#
# == details
# The symmetric GO similarity matrix can be treated as an adjacency matrix for constructing a graph/network.
# Thus, clustering the GO similarity matrix is identical to detecting clusters/modules/communities from the graph.
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
	g = igraph::graph_from_adjacency_matrix(mat, mode = "upper", weight = TRUE)

	imc = method(g, ...)
	as.vector(igraph::membership(imc))
}

# == title
# Cluster GO similarity matrix by mclust
#
# == param
# -mat The GO similarity matrix.
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
# Cluster GO similarity matrix by apcluster
#
# == param
# -mat The GO similarity matrix.
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


