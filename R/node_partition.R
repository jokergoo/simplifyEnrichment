
consensus_kmeans = function(mat, centers, km_repeats = 10) {
    partition_list = lapply(seq_len(km_repeats), function(i) {
        as.cl_hard_partition(kmeans(mat, centers))
    })
    partition_list = cl_ensemble(list = partition_list)
    partition_consensus = cl_consensus(partition_list)
    cl = as.vector(cl_class_ids(partition_consensus))
    if(length(unique(cl)) == 1) {
    	cl = partition_list[[1]]$.Data$cluster
    }
    cl
}

partition_by_pam = function(mat, k) {
    cluster::pam(mat, k)$clustering
}

partition_by_hclust = function(mat, k) {
    cutree(hclust(dist(mat)), k)
}

partition_by_skmeans = function(mat, k) {
    skmeans::skmeans(x = mat, k = k)$cluster
}
