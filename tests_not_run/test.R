

set.seed(123)
go_id = random_GO(500)
mat = GO_similarity(go_id)
simplify(mat)

dend = cluster_mat(mat)

render_dend(dend)

cl = cut_dend(dend)
plot_heatmap(mat, cl)


go_id = random_GO(500)
mat = GO_similarity(go_id)
compare_methods(mat)


pdf("~/test.pdf")
for(method in list_seriation_methods("matrix")) {
	qqcat("cluster by @{method}\n")
	order = seriate(mat, method = "method", control = list(verbose = TRUE))
	od = get_order(order)

	ht = Heatmap(mat, col = colorRamp2(c(0, 1), c("white", "red")),
		name = "Similarity",
		show_row_names = FALSE, show_column_names = FALSE, 
		row_order = od, column_order = od,
		column_title = method)
	draw(ht)
}
for(method in list_seriation_methods("dist")) {
	if(method %in% c("BBURCG", "BBWRCG", "Identity", "Random", "SA", "Spectral", "Spectral_norm")) next
	qqcat("cluster by @{method}\n")
	order = seriate(as.dist(1-mat), method = method, control = list(verbose = TRUE))
	od = get_order(order)

	ht = Heatmap(mat, col = colorRamp2(c(0, 1), c("white", "red")),
		name = "Similarity",
		show_row_names = FALSE, show_column_names = FALSE, 
		row_order = od, column_order = od,
		column_title = method)
	draw(ht)
}
dev.off()



bp_list = readRDS("bp_list.rds")

go_list = lapply(bp_list, function(x) {
	names(x)[p.adjust(x, "BH") < 0.05]
})

n = sapply(go_list, length)
go_list = go_list[n > 100 & n < 2000]

for(i in seq_along(go_list)) {
	go_id = go_list[[i]]
	mat = GO_similarity(go_id)
	compare_methods(mat)
}