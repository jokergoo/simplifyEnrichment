context("test dend_utils")

dend = as.dendrogram(hclust(dist(1:5)))
test_that("test dend_node_apply()", {
    h = dend_node_apply(dend, function(d, index) {
        h = attr(d, "height")
        if(is.null(index)) {
            names(h) = "top"
        } else {
            names(h) = paste(index, collapse = "")
        }
        h
    })

    expect_that(h[["top"]], equals(4))
    expect_that(h[["1"]], equals(1))
    expect_that(h[["2"]], equals(2))
   
})


test_that("test edit_node()", {
    dend2 = edit_node(dend, function(d, index) {
        attr(d, "index") = index
        d
    })

    expect_that(attr(dend2[[c(1, 2)]], "index"), equals(c(1, 2)))
})
