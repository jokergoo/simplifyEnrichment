# == title
# Similarity between terms based on the overlap of genes
#
# == param
# -gl A list of genes that are in the terms.
#
# == value
# A symmetric matrix
term_similarity = function(gl, method = c("jaccard", "kappa")) {
	all = unique(unlist(gl))
	gl = lapply(gl, function(x) as.numeric(factor(x, levels = all)))
	n = length(gl)

	mg = matrix(0, ncol = length(all), nrow = n)
	for(i in seq_len(n)) {
		mg[i, gl[[i]]] = 1
	}

	method = match.arg(method)[1]
	if(method == "jaccard") {
		mat = 1 - dist(mg, method = "binary")
		mat = as.matrix(mat)
	} else if(method == "kappa") {
		mat = matrix(1, nrow = n, ncol = n)
		for(i in seq(1, n - 1)) {
			for(j in seq(i+1, n)) {
				mat[i, j] = mat[j, i] = kappa(mg[i, ], mg[j, ])
			}
		}
	}
	diag(mat) = 1
	rownames(mat) = colnames(mat) = names(gl)
	return(mat)
}

# x and y are logical
kappa = function(x, y) {
	tab = length(x)
	oab = sum(x == y)/tab
	aab = (sum(x)*sum(y) + sum(!x)*sum(!y))/tab/tab
	(oab - aab)/(1 - aab)
}

#### similarity from enrichResult object ########

# == title
# Subset method of the enrichResult class
#
# == param
# -x A ``enrichResult`` object.
# -i Row indices.
# -j Column indices.
# -drop Please ignore.
#
"[.enrichResult" = function(x, i, j, drop = FALSE) {

	rownames(x@result) = x@result$ID

	if(missing(i) && missing(j)) {

	} else if(missing(i)) {
		x@result = x@result[, j, drop = FALSE]
	} else if(missing(j)) {
		x@result = x@result[i, , drop = FALSE]
	} else {
		x@result = x@result[i, j, drop = FALSE]
	}

	x@geneSets = x@geneSets[x@result$ID]
	return(x)
}

# == title
# Subset method of the enrichResult class
#
# == param
# -x A ``enrichResult`` object.
# -i Row indices.
#
subset_enrichResult = function(x, i) {
	rownames(x@result) = x@result$ID

	x@result = x@result[i, , drop = FALSE]

	x@geneSets = x@geneSets[x@result$ID]
	return(x)
}

# == title
# Similarity between terms in the enrichResult class
#
# == param
# -x A ``enrichResult`` object.
#
# == details
# The object is normally from the clusterProfiler, DOSE, meshes or ReactomePA package.
#
term_similarity_from_enrichResult = function(x) {
	term_similarity(x@geneSets[x@result$ID])
}

#### similarity directly from term IDs ######

# == title
# Similarity between KEGG terms
#
# == param
# -term_id A vector of KEGG IDs.
#
# == value
# A symmetric matrix
term_similarity_from_KEGG = function(term_id) {

	if(requireNamespace("clusterProfiler")) {
		stop_wrap("'clusterProfiler' package should be installed.")
	}

	species = gsub("^([a-zA-Z]+)(\\d+$)", "\\1", term_id[1])
	oe = try(KEGG_DATA <- getFromNamespace("prepare_KEGG", "clusterProfiler")(species, "KEGG", "kegg"), silent = TRUE)
	if(inherits(oe, "try-error")) {
		KEGG_DATA = getFromNamespace("get_data_from_KEGG_db", "clusterProfiler")(species)
	}

	gl = KEGG_DATA$PATHID2EXTID[term_id]

	term_similarity(gl)
}

# == title
# Similarity between Reactome terms
#
# == param
# -term_id A vector of Reactome IDs.
#
# == value
# A symmetric matrix
term_similarity_from_Reactome = function(term_id) {
	if(requireNamespace("reactome.db")) {
		stop_wrap("'reactome.db' package should be installed.")
	}

	all = as.list(reactome.db::reactomePATHID2EXTID)
	gl = all[term_id]

	term_similarity(gl)
}

# == title
# Similarity between MSigDB terms
#
# == param
# -term_id A vector of MSigDB gene set names.
# -category E.g., C1, C2, ...
#
# == value
# A symmetric matrix
term_similarity_from_MSigDB = function(term_id, category = NULL) {
	if(requireNamespace("msigdbr")) {
		stop_wrap("'msigdbr' package should be installed.")
	}

	m_df = msigdbr::msigdbr(species = "Homo sapiens", category = category)
	lt = split(m_df$entrez_gene, m_df$gs_name)
	gl = lt[term_id]

	term_similarity(gl)
}

# == title
# Similarity between terms from a gmt file
#
# == param
# -term_id A vector of terms.
# -gmt The path of the gmt file.
# -extract_term_id If the term ID in contained in the first column only as a substring,
#      setting a function to extract this substring.
#
# == value
# A symmetric matrix
term_similarity_from_gmt = function(term_id, gmt, extract_term_id = NULL) {
	ln = readLines(gmt)

	term_id = sapply(ln, function(x) x[1])
	if(!is.null(extract_term_id)) {
		term_id = extract_term_id(term_id)
		if(any(duplicated(term_id))) {
			stop_wrap("Duplicated term IDs generated after applying `extract_term_id`.")
		}
	}

	gl = lapply(ln, function(x) x[-(1:2)])
	names(gl) = term_id
	term_similarity(gl)
}
