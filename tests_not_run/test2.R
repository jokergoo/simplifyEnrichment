
################# test go id ###############3

set.seed(123)
go_id = random_GO(500)
mat = GO_similarity(go_id)
simplifyGO(mat)
simplifyGO(mat, word_cloud_grob_param=list(col = "white"), bg_gp = gpar(fill = "red", col = "blue"))


cutoff = seq(0.6, 0.9, by = 0.01)
s = NULL
for(i in seq_along(cutoff)) {
	qqcat("@{i}/@{length(cutoff)}...\n")
	cl = binary_cut(mat, cutoff = cutoff[i])
	s[i] = difference_score(mat, cl)
}

ht_clusters(mat, cl)

################# test Do id ###############3
do_id = random_DOSE(500)
mat = DOSE_similarity(do_id)


################# test reactome
reactome_sim = readRDS("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/rds/reactome_sim_kappa.rds")
mat = reactome_sim[[4]]

############ random term Msigdb
library(msigdbr)
m_df = msigdbr(category = "C5", subcategory = "BP")

gl = split(m_df$entrez_gene, m_df$gs_name)

gl = gl[sample(length(gl), 500)]

mat = term_similarity(gl, method = "kappa")


### test multiple times
go_id = random_GO(500)
mat = GO_similarity(go_id)
plot_binary_cut(mat)


### anno_word_cloud
gm = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds", package = "simplifyEnrichment"))
go_id = rownames(gm)
go_term = select(GO.db::GO.db, keys = go_id, columns = "TERM")$TERM


split = sample(letters[1:4], 100, replace = TRUE)
align_to = split(1:100, split)
term = lapply(letters[1:4], function(x) sample(go_term, sample(100:400, 1)))
names(term) = letters[1:4]

mat = matrix(rnorm(100*10), nrow = 100)
Heatmap(mat, cluster_rows = FALSE, row_split = split, 
	right_annotation = rowAnnotation(foo = anno_word_cloud(align_to, term)))
