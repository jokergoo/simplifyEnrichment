# == title
# Similarity between terms based on the overlap of genes
#
# == param
# -gl A list of genes that are in the terms.
# -method The similarity measurement.
# -all The universe set.
# -remove_negative If the value is negative, reset to zero
#
# == details
# The definition of the four similarity measurements can be found at https://jokergoo.github.io/simplifyEnrichment_supplementary/supplS01_coefficient_definition/supplS01_coefficient_definition.html .
#
# == value
# A symmetric matrix.
term_similarity = function(gl, method = c("kappa", "jaccard", "dice", "overlap"), all = NULL, remove_negative = TRUE) {
	if(is.null(all)) {
		all = unique(unlist(gl))
	} else {
		gl = lapply(gl, intersect, all)
	}
	gl = lapply(gl, function(x) as.numeric(factor(x, levels = all)))
	n = length(gl)

	mg = matrix(0, ncol = length(all), nrow = n)
	for(i in seq_len(n)) {
		mg[i, gl[[i]]] = 1
	}
	mg = as(mg, "sparseMatrix")

	method = match.arg(method)[1]
	if(method == "kappa") {
		mat = kappa_dist(mg, remove_negative = remove_negative)
	} else if(method == "overlap") {
		mat = overlap_dist(mg)
	} else {
		mat = proxyC::simil(mg, method = method)
	}

	mat = as.matrix(mat)
	diag(mat) = 1
	rownames(mat) = colnames(mat) = names(gl)
	return(mat)
}

# x and y are logical
kappa = function(x, y) {
	tab = length(x)
	oab = sum(x == y)/tab
	aab = (sum(x)*sum(y) + sum(!x)*sum(!y))/tab/tab
	k = (oab - aab)/(1 - aab)
	# if(k < 0) k = 0
	return(k)
}

# by rows
kappa_dist = function(m, remove_negative = TRUE) {
	tab = ncol(m)
	oab = proxyC::simil(m, method = "simple matching")
	m1 = rowSums(m)
	m2 = abs(rowSums(m - 1))
	aab = (outer(m1, m1) + outer(m2, m2))/tab/tab
	k = (oab - aab)/(1 - aab)
	if(remove_negative) k[k < 0] = 0
	return(k)
}

overlap_dist = function(m) {
	n = rowSums(m)
	proxyC::simil(m, method = "dice")*outer(n, n, "+")/2/outer(n, n, pmin)
}

# only for testing
overlap_single = function(x, y) {
	sum(x & y)/min(sum(x), sum(y))
}

#### similarity from enrichResult object ########

# title
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
# -x A ``enrichResult`` object from 'clusterProfiler' or other related packages.
# -i Row indices.
#
# == value
# Still a ``enrichResult`` object but with the selected subset of rows.
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
# -x A ``enrichResult`` object from 'clusterProfiler' or other related packages.
# -... Pass to `term_similarity`.
#
# == details
# The object is normally from the 'clusterProfiler', 'DOSE', 'meshes' or 'ReactomePA' package.
#
# == value
# A symmetric matrix.
term_similarity_from_enrichResult = function(x, ...) {
	term_similarity(x@geneSets[x@result$ID], ...)
}

#### similarity directly from term IDs ######

# == title
# Similarity between KEGG terms
#
# == param
# -term_id A vector of KEGG IDs, e.g., hsa001.
# -... Pass to `term_similarity`.
#
# == value
# A symmetric matrix.
term_similarity_from_KEGG = function(term_id, ...) {

	check_pkg("clusterProfiler", bioc = TRUE)

	species = gsub("^([a-zA-Z]+)(\\d+$)", "\\1", term_id[1])
	oe = try(KEGG_DATA <- getFromNamespace("prepare_KEGG", "clusterProfiler")(species, "KEGG", "kegg"), silent = TRUE)
	if(inherits(oe, "try-error")) {
		KEGG_DATA = getFromNamespace("get_data_from_KEGG_db", "clusterProfiler")(species)
	}

	gl = KEGG_DATA$PATHID2EXTID[term_id]

	term_similarity(gl, ...)
}

# == title
# Similarity between Reactome terms
#
# == param
# -term_id A vector of Reactome IDs.
# -... Pass to `term_similarity`.
#
# == value
# A symmetric matrix.
term_similarity_from_Reactome = function(term_id, ...) {

	check_pkg("reactome.db", bioc = TRUE)

	all = as.list(reactome.db::reactomePATHID2EXTID)
	gl = all[term_id]

	term_similarity(gl, ...)
}


# == title
# Similarity between MSigDB terms
#
# == param
# -term_id A vector of MSigDB gene set names.
# -category E.g., 'C1', 'C2', pass to `msigdbr::msigdbr`.
# -subcategory E.g., 'CGP', 'BP', pass to `msigdbr::msigdbr`.
# -... Pass to `term_similarity`.
#
# == value
# A symmetric matrix.
term_similarity_from_MSigDB = function(term_id, category = NULL, subcategory = NULL, ...) {

	check_pkg("msigdbr", bioc = FALSE)

	m_df = msigdbr::msigdbr(species = "Homo sapiens", category = category, subcategory = subcategory)
	if(all(grepl("^M\\d+$", term_id))) {
		lt = split(m_df$entrez_gene, m_df$gs_id)
	} else {
		lt = split(m_df$entrez_gene, m_df$gs_name)
	}
	gl = lt[term_id]

	term_similarity(gl, ...)
}

# == title
# Similarity between terms from a gmt file
#
# == param
# -term_id A vector of terms.
# -gmt The path of the gmt file.
# -extract_term_id If the term ID in the first column only as a substring,
#      setting a function to extract this substring.
# -... Pass to `term_similarity`.
#
# == value
# A symmetric matrix.
term_similarity_from_gmt = function(term_id, gmt, extract_term_id = NULL, ...) {
	ln = readLines(gmt)
	ln = strsplit(ln, "\t")
	
	term_id = sapply(ln, function(x) x[1])
	if(!is.null(extract_term_id)) {
		term_id = extract_term_id(term_id)
		if(any(duplicated(term_id))) {
			stop_wrap("Duplicated term IDs generated after applying `extract_term_id`.")
		}
	}

	gl = lapply(ln, function(x) x[-(1:2)])
	names(gl) = term_id
	term_similarity(gl, ...)
}


gene_weight = function(gl) {
	all = unique(unlist(gl))
	gl = lapply(gl, function(x) as.numeric(factor(x, levels = all)))
	n = length(gl)

	mg = Matrix(0, ncol = length(all), nrow = n)
	for(i in seq_len(n)) {
		mg[i, gl[[i]]] = 1
	}

	log(1 + nrow(mg)/colSums(mg))
}
