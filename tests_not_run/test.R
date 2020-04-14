

set.seed(123)
go_id = random_GO(500)
mat = GO_similarity(go_id)
simplify(mat)

dend = cluster_mat(mat)

render_dend(dend)

cl = cut_dend(dend)
plot_heatmap(mat, cl)
