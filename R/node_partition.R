
# == title
# Partition by kmeans
#
# == param
# -mat The similarity matrix.
# -n_repeats Number of repeated runs of k-means.
#
# == details
# Since k-means clustering brings randomness, this function performs
# k-means clustering several times and uses the final consensus partitioning.
#
# This function is used to set to the ``partition_fun`` argument in `binary_cut`.
partition_by_kmeans = function(mat, n_repeats = 10) {
    partition_list = lapply(seq_len(n_repeats), function(i) {
        as.cl_hard_partition(kmeans(mat, 2))
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
#
# == details
# The clustering is performed by `cluster::pam` with setting ``pamonce`` argument to 5.
#
# This function is used to set to the ``partition_fun`` argument in `binary_cut`.
partition_by_pam = function(mat) {
    fit = pam(mat, 2, pamonce = 5)
    fit$clustering
}


# == title
# Partition by hclust
#
# == param
# -mat The similarity matrix.
#
# == details
# The "ward.D2" clusering method was used.
#
# This function is used to set to the ``partition_fun`` argument in `binary_cut`.
partition_by_hclust = function(mat) {
    cutree(stats::hclust(stats::dist(mat), method = "ward.D2"), 2)
}

# partition_by_skmeans = function(mat) {
#     skmeans::skmeans(x = mat, k = 2)$cluster
# }


# == title
# Partition by kmeans++
#
# == param
# -mat The similarity matrix.
#
# == details
# This function is used to set to the ``partition_fun`` argument in `binary_cut`.
partition_by_kmeanspp = function(mat) {
    check_pkg("flexclust", bioc = FALSE)
    
    cl = flexclust::kcca(mat, k = 2, 
        family = flexclust::kccaFamily("kmeans"),
        control = list(initcent = "kmeanspp"))@cluster
    cl
}
