context("test register_methods")

test_that("test register_methods", {
    expect_that(length(all_clustering_methods()), equals(11))
    remove_clustering_methods(c("binary_cut"))
    expect_that(length(all_clustering_methods()), equals(10))
    reset_clustering_methods()
	expect_that(length(all_clustering_methods()), equals(11))
})
