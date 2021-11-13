context("test register_methods")

test_that("test register_methods", {
    expect_that(length(all_clustering_methods()), equals(12))
    remove_clustering_methods(c("binary_cut"))
    expect_that(length(all_clustering_methods()), equals(11))
    reset_clustering_methods()
	expect_that(length(all_clustering_methods()), equals(12))
})
