context("test test_sim")

library(Matrix)
library(proxyC)

m = matrix(0, ncol = 100, nrow = 4)
m[sample(400, 50)] = 1

sm = as(m, "sparseMatrix")

jm = proxyC::simil(sm, method = "jaccard")

test_that("test jaccard similarity", {
	expect_that(jm[1, 2], equals(sum(m[1, ] & m[2, ])/sum(m[1, ] | m[2, ])))
	expect_that(jm[1, 1], equals(1))
	expect_that(jm[2, 1], equals(sum(m[2, ] & m[1, ])/sum(m[2, ] | m[1, ])))
	expect_that(jm[2, 2], equals(1))
})

# x and y are logical
kappa = function(x, y) {
	tab = length(x)
	oab = sum(x == y)/tab
	aab = (sum(x)*sum(y) + sum(!x)*sum(!y))/tab/tab
	k = (oab - aab)/(1 - aab)
	if(k < 0) k = 0
	return(k)
}

# by rows
kappa_dist = function(m) {
	tab = ncol(m)
	oab = proxyC::simil(m, method = "simple matching")
	m1 = rowSums(m)
	m2 = abs(rowSums(m - 1))
	aab = (outer(m1, m1) + outer(m2, m2))/tab/tab
	k = (oab - aab)/(1 - aab)
	k[k < 0] = 0
	return(k)
}

km = kappa_dist(sm)
test_that("test jaccard similarity", {
	expect_that(km[1, 2], equals(kappa(m[1, ], m[2, ])))
	expect_that(km[2, 1], equals(kappa(m[2, ], m[1, ])))
}