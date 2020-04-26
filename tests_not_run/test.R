

set.seed(123)
go_id = random_GO(500)
mat = GO_similarity(go_id)
simplifyGO(mat)

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
