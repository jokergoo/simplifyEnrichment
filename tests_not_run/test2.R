
################# test go id ###############3

set.seed(123)
go_id = random_GO(500)
mat = GO_similarity(go_id)
simplifyGO(mat)


cutoff = seq(0.6, 0.9, by = 0.01)
s = NULL
for(i in seq_along(cutoff)) {
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
