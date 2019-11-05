hsGO = godata('org.Hs.eg.db', ont = "BP")
all_go_id = hsGO@geneAnno$GO

go_id = sample(all_go_id, 500)

mat = get_GO_sim_mat(go_id)

dend = cluster_mat(mat)
make_rule(dend, T)
plot_dend(dend, mat)

size = dend_node_apply(dend, function(x) attr(x, "member"))
score = dend_node_apply(dend, function(x) attr(x, "score"))

env = new.env()
env$dend = dend
parent_score = dend_node_apply(dend, function(x) {
	index = attr(x, "index")
	if(is.null(index)) {
		return(NA)
	} else if(length(index) == 1) {
		return(attr(env$dend, "score") - attr(dend, "score"))
	} else {
		return(attr(env$dend[[ index[-length(index)] ]], "score") - attr(dend, "score"))
	}
})


