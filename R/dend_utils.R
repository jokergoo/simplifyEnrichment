
next_k = local({
	k = 0
	function(reset = FALSE) {
		if(reset) {
			k <<- 0
		} else {
			k <<- k + 1
		}
		k
	}
})

# == title
# Apply functions on every node in a dendrogram
#
# == param
# -dend A dendrogram.
# -fun A self-defined function.
#
# == details
# The function returns a vector or a list as the same length as the number of nodes in the dendrogram.
#
# The self-defined function can have one single argument which is the sub-dendrogram at a certain node.
# E.g. to get the number of members at every node:
#
#     dend_node_apply(dend, function(d) attr(d, "members"))
#
# The self-defined function can have a second argument, which is the index of current sub-dendrogram in 
# the complete dendrogram. E.g. ``dend[[1]]`` is the first child node of the complete dendrogram and
# ``dend[[c(1, 2)]]`` is the second child node of ``dend[[1]]``, et al. This makes that at a certain node,
# it is possible to get information of its child nodes and parent nodes. 
#
#     dend_node_apply(dend, function(d, index) {
#         dend[[c(index, 1)]] # is the first child node of d, or simply d[[1]]
#         dend[[index[-length(index)]]] # is the parent node of d
#         ...
#     })
#
# Note for the top node, the value of ``index`` is ``NULL``.
#
# == value
# A vector or a list, depends on whether ``fun`` returns a scalar or more complex values.
#
# == example
# mat = matrix(rnorm(100), 10)
# dend = as.dendrogram(hclust(dist(mat)))
# # number of members on every node
# dend_node_apply(dend, function(d) attr(d, "members"))
# # the depth on every node
# dend_node_apply(dend, function(d, index) length(index))
dend_node_apply = function(dend, fun) {

	next_k(reset = TRUE)

	assign_to = function(env, k, v) {
		n = length(env$var)
		if(n == 0) {
			env$var = list()
		}
		env$var[[k]] = v
	}


	if(length(as.list(formals(fun))) == 1) {
		fun2 = fun
		fun = function(d, index) fun2(d)
	}

	env = new.env()
	.do = function(dend, fun, index) {

		if(is.null(index)) {
			if(is.leaf(dend)) {
				assign_to(env, next_k(), fun(dend, index))
				return(NULL)
			} else {
				assign_to(env, next_k(), fun(dend, index))
			}
		} else {
			assign_to(env, next_k(), fun(dend[[index]], index))

			if(is.leaf(dend[[index]])) {
				return(NULL)
			}
		}

		if(is.null(index)) {
			n = length(dend)
		
		} else {
			n = length(dend[[index]])
		}
		for(i in seq_len(n)) {
			.do(dend, fun, c(index, i))
		}
	}

	.do(dend, fun, NULL)

	var = env$var
	if(all(vapply(var, is.atomic, TRUE))) {
		if(all(vapply(var, length, 0) == 1)) {
			var = unlist(var)
		}
	}

	return(var)
}


# == title
# Modify nodes in a dendrogram
#
# == param
# -dend A dendrogram.
# -fun A self-defined function.
#
# == details
# if ``fun`` only has one argument, it is basically the same as `stats::dendrapply`,
# but it can have a second argument which is the index of the node in the dendrogram,
# which makes it possible to get information of child nodes and parent nodes for
# a specific node.
#
# As an example, we first assign random values to every node in the dendrogram:
#
#     mat = matrix(rnorm(100), 10)
#     dend = as.dendrogram(hclust(dist(mat)))
#     dend = edit_node(dend, function(d) {attr(d, 'score') = runif(1); d})
#
# Then for every node, we take the maximal absolute difference to all its child nodes
# and parent node as the attribute ``abs_diff``
#
#     dend = edit_node(dend, function(d, index) {
#         n = length(index)
#         s = attr(d, "score")
#         if(is.null(index)) {  # d is the top node
#             s_children = sapply(d, function(x) attr(x, "score"))
#             s_parent = NULL
#         } else if(is.leaf(d)) { # d is the leaf
#             s_children = NULL
#             s_parent = attr(dend[[index[-n]]], "score")
#         } else {
#             s_children = sapply(d, function(x) attr(x, "score"))
#             s_parent = attr(dend[[index[-n]]], "score")
#         }
#         abs_diff = max(abs(s - c(s_children, s_parent)))
#         attr(d, "abs_diff") = abs_diff
#         return(d)
#     })
#
# == value
# A dendrogram object.
edit_node = function(dend, fun = function(d, index) d) {

	# breadth first

	env = new.env()
	env$dend = dend

	if(length(as.list(formals(fun))) == 1) {
		fun2 = fun
		fun = function(d, index) fun2(d)
	}

	traversal_dend = function(env, index = NULL) {
		# index is null means it is the top node
		if(is.null(index)) {
			d = env$dend
			if(is.leaf(d)) {
				env$dend = fun(d, NULL)
				return(NULL)
			} else {
				env$dend = fun(d, NULL)

				for(i in seq_along(d)) {
					traversal_dend(env, i)
				}
			}
		} else {
			d = env$dend[[index]]
			if(is.leaf(d)) {
				env$dend[[index]] = fun(d, index)
				return(NULL)
			} else {
				env$dend[[index]] = fun(d, index)
				for(i in seq_along(d)) {
					traversal_dend(env, c(index, i))
				}
			}
		}
	}

	traversal_dend(env)
	return(env$dend)
}
