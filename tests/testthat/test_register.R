context("test register_methods")

test_that("test register_methods", {
    expect_that(length(all_clustering_methods()), equals(9))
    remove_clustering_methods(c("binary_cut"))
    expect_that(length(all_clustering_methods()), equals(8))
    reset_clustering_methods()
	expect_that(length(all_clustering_methods()), equals(9))

	m = matrix(rnorm(100), 10)
	cluster_terms(m, method = "kmeans", control = list(max_k = 4))
})
