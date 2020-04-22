
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
# Apply on every node in a dendrogram
#
# == param
# -dend A dendrogram.
# -fun A self-defined function.
#
# == details
# The function return a vector or a list with the same length as the number of nodes in the denrogram.
#
# The self-defined function can have one single argument which is the sub-dendrogram at a certain node,
# e.g. to get the number of members at each node:
#
#     dend_node_apply(dend, function(d) attr(d, "members"))
#
# The self-defined function can have a second argument, which is the index of current sub-dendrogram in 
# the complete dendrogram. E.g. ``dend[[1]]`` is the first child node of the complete dendrogram and
# ``dend[[c(1, 2)]]`` is the second child node of ``dend[[1]]``, et al. This makes that at a certain node,
# it is possible to get informatino of its children nodes and parent nodes. 
#
#     dend_node_apply(dend, function(d, index) {
#         d[[c(index, 1)]] # is the first child node of d
#         d[[index[-length(index)]]] # is the parent node of d
#         ...
#     })
#
# == value
# A vector or a list, depends on whether ``fun`` returns a scalar or more complex values.
#
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

		n = length(dend)
		for(i in seq_len(n)) {
			.do(dend, fun, c(index, i))
		}
	}

	.do(dend, fun, NULL)

	var = env$var
	if(all(sapply(var, is.atomic))) {
		if(all(sapply(var, length) == 1)) {
			var = unlist(var)
		}
	}

	return(var)
}


# == title
# Modify every node in a dendrogram
#
# == param
# -dend A dendrogram.
# -fun A self-defined function.
#
# == details
# if ``fun`` only has one argument, it is basically the same as `stats::dendrapply`,
# but it can have a second argument which is the index of the node in the dendrogram,
# which makes it possible to get information of children nodes and parent nodes for
# a specified node.
#
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
