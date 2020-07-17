
# # for e.g. kmeans
# binary_cut_stability = function(mat, nrun = 100, ...) {
# 	partition_list = list()
	
# 	for(i in seq_len(nrun)) {
# 		qqcat("@{i}/@{nrun} runs...\n")
# 		partition_list[[i]] = as.cl_hard_partition(binary_cut(mat, ...))
# 	}

# 	partition_list = cl_ensemble(list = partition_list)
#     partition_consensus = cl_consensus(partition_list)
    
#     class_ids = as.vector(cl_class_ids(partition_consensus))

# 	membership_each = do.call("cbind", lapply(seq_along(partition_list), function(i) {
# 		x = partition_list[[i]]
# 		class = as.vector(cl_class_ids(x))
# 		map = relabel_class(class, class_ids)
# 		class = as.numeric(map[as.character(class)])
# 		as.integer(class)
# 	}))

# 	membership_each2 = membership_each
# 	membership_each2[is.na(membership_each2)] = as.integer(0)
# 	consensus_mat = cola:::get_consensus_matrix(membership_each2)

# 	concordance = apply(membership_each2, 2, function(x) sum(x == class_ids)/length(x))
 	
#  	list(class_ids = class_ids, consensus_mat = consensus_mat, membership_mat = membership_each2,
#  		concordance = concordance)
# }
