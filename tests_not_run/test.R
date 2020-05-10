

set.seed(123)
go_id = random_GO(500)
mat = GO_similarity(go_id)
simplifyGO(mat)

cl = cluster_terms(mat)
ht_clusters(mat, cl)

dend = cluster_mat(mat)

render_dend(dend)

cl = cut_dend(dend)
plot_heatmap(mat, cl)


for(i in 1:25) {
	go_id = random_GO(500)
	mat = GO_similarity(go_id)
	lt = simplifyGO(mat)

	compare_methods(mat)


lt = list()
n = list()
v = list()
for(i in 1:20) {
	print(i)
	go_id = random_GO(500)
	mat = GO_similarity(go_id)
	clt = compare_methods_make_clusters(mat, "all")
	x = sapply(clt, function(x) difference_score(mat, x))
	lt[[i]] = x
	n[[i]] = sapply(clt, function(x) length(unique(x)))
	v[[i]] = sapply(clt, function(x) block_mean(mat, x))
}

df = do.call(rbind, lapply(v, function(x) data.frame(method = names(x), value = x)))
ggplot(df, aes(x=method, y=value)) + 
  geom_violin(trim = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

bp_list = readRDS("../simplifyGO_figures/bp_list.rds")

go_list = lapply(bp_list, function(x) {
	names(x)[p.adjust(x, "BH") < 0.05]
})

n = sapply(go_list, length)
go_list = go_list[n > 100 & n < 2000]

lt = list()
n = list()
v = list()
for(i in seq_along(go_list)[1:20]) {
	go_id = go_list[[i]]
	mat = GO_similarity(go_id)
	# compare_methods(mat)
	clt = compare_methods_make_clusters(mat, "all")
	x = sapply(clt, function(x) difference_score(mat, x))
	lt[[i]] = x
	n[[i]] = sapply(clt, function(x) length(unique(x)))
	v[[i]] = sapply(clt, function(x) block_mean(mat, x))
}
df = do.call(rbind, lapply(v, function(x) data.frame(method = names(x), value = x)))
ggplot(df, aes(x=method, y=value)) + 
  geom_violin(trim = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


for(i in seq_along(go_list)[1:20]) {
	go_id = go_list[[i]]
	mat = GO_similarity(go_id)
	dend = cluster_mat(mat)
	x = dend_node_apply(dend, function(d) attr(d, "members"))
	y = dend_node_apply(dend, function(d) attr(d, "score2"))
	plot(x, y, log = "x")
	readline("enter: ")
}

d = NULL
cutoff = seq(0.6, 0.95, by = 0.01)
for(i in seq_along(cutoff)) {
	cl = binary_cut(mat, cutoff = cutoff[i])
	d[i] = difference_score(mat, cl)
}
plot(cutoff, d)



go_id = random_GO(500)
mat = GO_similarity(go_id)
clt = compare_methods_make_clusters(mat)
g = graph_from_adjacency_matrix(mat, mode = "upper", weighted = TRUE)

sapply(clt, function(x) modularity(g, x, weights = E(g)$weight))


##########

go_id = random_GO(500)
semData <- godata('org.Hs.eg.db', ont = "BP")
system.time(termSim(go_id, go_id, method = "Wang", semData = semData))
system.time({
	m = matrix(nrow = 500, ncol = 500)
	for(i in 1:500) {
		for(j in 1:500) {
			m[i, j] = termSim(go_id[i], go_id[j], method = "Wang", semData = semData)
		}
	}
})

split_by_block = function(n, size) {
	size = min(c(n, size))
	REST <- n%%size
    LARGE <- n - REST
    NBLOCKS <- n%/%size
    GROUP <- rep(1:NBLOCKS, each = size)
    if (REST > 0)
        GROUP <- c(GROUP, rep(NBLOCKS + 1, REST))
    split(1:n, GROUP)
}

NCOL <- length(go_id)
SPLIT = split_by_block(NCOL, size)
COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
COMBS <- t(apply(COMBS, 1, sort))
COMBS <- unique(COMBS)

system.time(lt <- mclapply(seq_len(nrow(COMBS)), function(i) {
	ind1 = SPLIT[[ COMBS[i, 1] ]]
	ind2 = SPLIT[[ COMBS[i, 2] ]]
	termSim(go_id[ind1], go_id[ind2], method = "Wang", semData = semData)
}, mc.cores = 2))

m = matrix(nrow = 500, ncol = 500)
for(i in seq_len(nrow(COMBS))) {
	ind1 = SPLIT[[ COMBS[i, 1] ]]
	ind2 = SPLIT[[ COMBS[i, 2] ]]
	if(COMBS[i, 1] == COMBS[i, 2]) {
		m[ind1, ind2] = lt[[i]]
	} else {
		m[ind1, ind2] = lt[[i]]
		m[ind2, ind1] = lt[[i]]
	}
}


