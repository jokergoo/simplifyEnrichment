
.ENV = new.env()
.ENV$ALL_CLUSTERING_FUN = list()
.ENV$ALL_CLUSTERING_METHODS = NULL

get_clustering_method = function(method, control = list()) {
	if(!method %in% .ENV$ALL_CLUSTERING_METHODS) {
		stop_wrap(qq("clustering method '@{method}' has not been defined yet."))
	}
	fun = .ENV$ALL_CLUSTERING_FUN[[method]]

	fun2 = function(mat) {
		cl = do.call(fun, c(list(mat), control))
		if(is.atomic(cl)) {
			if(length(cl) != nrow(mat)) {
				stop_wrap("Length of clusterings should be the same as number of matrix rows.")
			}
		} else {
			stop_wrap(qq("clustering method '@{method}' should return an atomic vector."))
		}
		cl = as.numeric(as.factor(cl))

		return(cl)
	}
	
	return(fun2)
}

# == title
# Register new clustering methods
#
# == param
# -... A named list of clustering functions, see Details
#
# == details
# The user-defined functions should accept at least one argument which is the input matrix. 
# The second optional argument should always be ``...`` so that parameters
# for the clustering function can be passed by ``control`` argument from `cluster_terms` or `simplifyGO`.
# If users forget to add ``...``, it is added internally.
#
# Please note, the user-defined function should automatically identify the optimized
# number of clusters.
#
# The function should return a vector of cluster labels. Internally it is converted to numeric labels.
#
# == value
# No value is returned.
#
register_clustering_methods = function(...) {
	
	lt = list(...)
	lt = lapply(lt, function(fun) {
		# just in case people forgot to add the ...
		if(length(formals(fun)) == 1) {
			function(mat, ...) {
				fun(mat)
			}
		} else {
			fun
		}
	})
	lt1 = lt[intersect(names(lt), .ENV$ALL_CLUSTERING_METHODS)]
	lt2 = lt[setdiff(names(lt), .ENV$ALL_CLUSTERING_METHODS)]
	if(length(lt1)) .ENV$ALL_CLUSTERING_FUN[names(lt1)] = lt1
	if(length(lt2)) .ENV$ALL_CLUSTERING_FUN = c(.ENV$ALL_CLUSTERING_FUN, lt2)
	.ENV$ALL_CLUSTERING_METHODS = names(.ENV$ALL_CLUSTERING_FUN)
}

# == title
# All clustering methods
#
# == details
# The default clustering methods are:
#
# -``binary_cut`` see `binary_cut`.
# -``kmeans`` see `cluster_by_kmeans`.
# -``dynamicTreeCut`` see `cluster_by_dynamicTreeCut`.
# -``mclust`` see `cluster_by_mclust`. By default it is not included.
# -``apcluster`` see `cluster_by_apcluster`.
# -``fast_greedy`` see `cluster_by_igraph`.
# -``leading_eigen`` see `cluster_by_igraph`.
# -``louvain`` see `cluster_by_igraph`.
# -``walktrap`` see `cluster_by_igraph`.
#
# == value
# A vector of method names
#
# == example
# all_clustering_methods()
all_clustering_methods = function() {
	x = .ENV$ALL_CLUSTERING_METHODS
	return(x)
}

# == title
# Remove clustering methods
#
# == param
# -method A vector of method names.
#
# == value
# No value is returned.
remove_clustering_methods = function(method) {
	nm_keep = setdiff(.ENV$ALL_CLUSTERING_METHODS, method)
	.ENV$ALL_CLUSTERING_FUN = .ENV$ALL_CLUSTERING_FUN[nm_keep]
	.ENV$ALL_CLUSTERING_METHODS = nm_keep
}

register_clustering_methods(
	kmeans = function(mat, ...) cluster_by_kmeans(mat, ...),
	dynamicTreeCut = function(mat, ...) cluster_by_dynamicTreeCut(mat, ...),
	mclust = function(mat, ...) cluster_by_mclust(mat, ...),
	apcluster = function(mat, ...) cluster_by_apcluster(mat, ...),
	fast_greedy = function(mat, ...) cluster_by_igraph(mat, method = "cluster_fast_greedy", ...),
	leading_eigen = function(mat, ...) cluster_by_igraph(mat, method = "cluster_leading_eigen", ...),
	louvain = function(mat, ...) cluster_by_igraph(mat, method = "cluster_louvain", ...),
	walktrap = function(mat, ...) cluster_by_igraph(mat, method = "cluster_walktrap", ...),
	binary_cut = function(mat, ...) binary_cut(mat, ...)
)

# == title
# Reset to default clustering methods
#
# == details
# The default methods are:
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
# == value
# No value is returned.
#
# == example
# all_clustering_methods()
# remove_clustering_methods(c("kmeans", "mclust"))
# all_clustering_methods()
# reset_clustering_methods()
# all_clustering_methods()
reset_clustering_methods = function() {
	remove_clustering_methods(all_clustering_methods())
	register_clustering_methods(
		kmeans = function(mat, ...) cluster_by_kmeans(mat, ...),
		dynamicTreeCut = function(mat, ...) cluster_by_dynamicTreeCut(mat, ...),
		mclust = function(mat, ...) cluster_by_mclust(mat, ...),
		apcluster = function(mat, ...) cluster_by_apcluster(mat, ...),
		fast_greedy = function(mat, ...) cluster_by_igraph(mat, method = "cluster_fast_greedy", ...),
		leading_eigen = function(mat, ...) cluster_by_igraph(mat, method = "cluster_leading_eigen", ...),
		louvain = function(mat, ...) cluster_by_igraph(mat, method = "cluster_louvain", ...),
		walktrap = function(mat, ...) cluster_by_igraph(mat, method = "cluster_walktrap", ...),
		binary_cut = function(mat, ...) binary_cut(mat, ...)
	)
}
