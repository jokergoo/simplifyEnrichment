
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

dend_node_apply = function(dend, fun) {

	next_k(reset = TRUE)

	assign_to = function(env, k, v) {
		n = length(env$var)
		if(n == 0) {
			if(is.atomic(v)) {
				env$var[k] = v
			} else {
				env$var = list()
				env$var[[k]] = v
			}
		} else {
			if(is.list(env$var)) {
				env$var[[k]] = v
			} else {
				if(is.atomic(v)) {
					env$var[k] = v
				} else {
					env$var = lapply(1:n, function(i) env$var[i])
					env$var[[k]] = v
				}
			}
		}
		
	}

	if(length(as.list(formals(fun))) == 1) {
		formals(fun) = alist(dend = , index = NULL)
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

	return(env$var)
}

# breadth first
edit_node = function(dend, fun = function(dend, index) dend) {

	env = new.env()
	env$dend = dend

	if(length(as.list(formals(fun))) == 1) {
		formals(fun) = alist(dend = , index = NULL)
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