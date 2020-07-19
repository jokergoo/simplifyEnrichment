
# == title
# Partition by kmeans
#
# == param
# -mat The similarity matrix.
# -k Number of clusters.
# -n_repeats Number of repeated run of k-means.
#
# == details
# Since k-means clustering brings randomness, this function performs
# k-means clustering several times and uses the final consensus clustering.
partition_by_kmeans = function(mat, k, n_repeats = 10) {
    partition_list = lapply(seq_len(n_repeats), function(i) {
        as.cl_hard_partition(kmeans(mat, k))
    })
    partition_list = cl_ensemble(list = partition_list)
    partition_consensus = cl_consensus(partition_list)
    cl = as.vector(cl_class_ids(partition_consensus))
    if(length(unique(cl)) == 1) {
    	cl = partition_list[[1]]$.Data$cluster
    }
    cl
}

# == title
# Partition by PAM
#
# == param
# -mat The similarity matrix.
# -k Number of clusters.
#
# == details
# The clustering is performed by `cluster::pam`.
partition_by_pam = function(mat, k) {
    pam(mat, k)$clustering
}

# == title
# Partition by hclust
#
# == param
# -mat The similarity matrix.
# -k Number of clusters.
#
partition_by_hclust = function(mat, k) {
    cutree(hclust(dist(mat)), k)
}

# partition_by_skmeans = function(mat, k) {
#     skmeans::skmeans(x = mat, k = k)$cluster
# }
