
# == title
# Partition by kmeans
#
# == param
# -mat The similarity matrix.
# -k Number of clusters.
# -n_repeats Number of repeated runs of k-means.
#
# == details
# Since k-means clustering brings randomness, this function performs
# k-means clustering several times and uses the final consensus partitioning.
#
# This function is used to set to the ``partition_fun`` argument in `binary_cut`.
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
# The clustering is performed by `cluster::pam` with setting ``pamonce`` argument to 5.
#
# This function is used to set to the ``partition_fun`` argument in `binary_cut`.
partition_by_pam = function(mat, k) {
    pam(mat, k, pamonce = 5)$clustering
}

# == title
# Partition by hclust
#
# == param
# -mat The similarity matrix.
# -k Number of clusters.
#
# == details
# The clusering method and the distance method both take the defaults of `stats::hclust` and `stats::dist`.
#
# This function is used to set to the ``partition_fun`` argument in `binary_cut`.
partition_by_hclust = function(mat, k) {
    cutree(stats::hclust(stats::dist(mat)), k)
}

# partition_by_skmeans = function(mat, k) {
#     skmeans::skmeans(x = mat, k = k)$cluster
# }

# == title
# Partition by multiple methods
#
# == param
# -mat The similarity matrix.
# -k Number of clusters.
#
# == details
# It runs three partitionning methods: `partition_by_pam`, `partition_by_kmeans`
# and `partition_by_hclust` and pick the classification with the highest
# difference scores (calculated by `difference_score`).
#
# This function is used to set to the ``partition_fun`` argument in `binary_cut`.
partition_by_max_ds = function(mat, k) {
    cl_list = list(by_pam = partition_by_pam(mat, k),
                   by_kmeans = partition_by_kmeans(mat, k),
                   by_hclust = partition_by_hclust(mat, k))
    ds = sapply(cl_list, function(cl) difference_score(mat, cl))
    ds[is.na(ds)] = -Inf
    i = which.max(ds)
    # qqcat("select partition method '@{names(cl_list)[i]}'.\n")
    cl_list[[i]]
}
