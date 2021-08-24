

# gl = list()
# for(db in c("reactome", "pid", "panther", "kegg", "pathbank")) {
# 	gl[[db]] = import_pathway_common(db)
# }


get_pathway_list = function(db) {
	# pl = readRDS(system.file("extdata", qq("pathway_network_@{db}.rds"), package = "simplifyEnrichment"))
	pl = readRDS(qq("~/project/development/simplifyEnrichment/inst/extdata/pathway_network_@{db}.rds"))
	map = pl[[2]]
	pl = pl[[1]]
	names(pl) = map[names(pl)]
	# a list of igraph object
	return(pl)
}

guess_pathway_database = function(id) {
	id = sample(id, 1)

	if(grepl("^hsa\\d+$", id)) {
		return("kegg")
	} else if(grepl("^R-\\w+-\\d+$", id)) {
		return("reactome")
	} else if(grepl("^P\\d+$", id)) {
		return("panther")
	} else if(grepl("^SMP\\d+$", id)) {
		return("pathbank")
	} else {
		return(NA)
	}
}

is_pathway_id = function(id) {
	!is.na(guess_pathway_database(id))
}

random_KEGG = function(n) {

	pl = get_pathway_list("kegg")
	id = names(pl)
	sample(id, min(n, length(id)))
}

random_Reactome = function(n) {

	pl = get_pathway_list("reactome")
	id = names(pl)
	sample(id, min(n, length(id)))
}

random_PANTHER = function(n) {
	pl = get_pathway_list("panther")
	id = names(pl)
	sample(id, min(n, length(id)))
}

random_PathBank= function(n) {
	pl = get_pathway_list("pathbank")
	id = names(pl)
	sample(id, min(n, length(id)))
}

# -pathway_id A vector of pathway IDs
pathway_similarity = function(pathway_id, centrality = "degree") {

	database = guess_pathway_database(pathway_id)
	if(is.na(database)) {
		stop_wrap("Cannot decide which pathway database the IDs are from.")
	}

	pl = get_pathway_list(database)
	cn = intersect(names(pl), pathway_id)

	if(length(cn) == 0) {
		stop_wrap("Cannot find proper pathways.")
	}

	if(length(cn) != length(pl)) {
		message_wrap(qq("@{length(pathway_id) - length(cn)}/@{length(pathway_id)} @{database} pathways removed."))
	}

	pl = pl[cn]
	pathway_id = cn

	centrality = switch(centrality, 
		degree = igraph::degree,
		closeness = igraph::closeness,
		reach = reach,
		spread = spread,
		betweenness = igraph::betweenness
	)

	lt = lapply(pl, function(g) {
		centrality(g)
	})

	n = length(lt)
	m = matrix(1, nrow = n, ncol = n)

	for(i in seq(1, n-1)) {
		for(j in seq(i+1, n)) {
			m[i, j] = weighted_jaccard(lt[[i]], lt[[j]])
			m[j, i] = m[i, j]
		}
	}

	rownames(m) = colnames(m) = cn
	return(m)
}

reach = function(graph, weights = E(graph)$weight, mode = c("all", "in", "out")) {
    mode = match.arg(mode)[1]
    sp = shortest.paths(graph, weights = weights, mode = mode)
    s = apply(sp, 1, function(x) {
        if (all(x == Inf)) {
            return(0)
        } else {
            return(max(x[x != Inf]))
        }
    })
    return(s)
}

spread = function(graph, mode = c("all", "in", "out"), weights = E(graph)$weight, f = function(x) 1/x) {
    mode = match.arg(mode)[1]
    sp = shortest.paths(graph, mode = mode, weights = weights)
    s = apply(sp, 1, function(x) {
        return(sum(f(x[x > 0])))
    })
    return(s)
}

weighted_jaccard = function(x, y) {
	nm = union(names(x), names(y))
	n = length(nm)
	y2 = x2 = numeric(n)
	names(y2) = names(x2) = nm
	x2[names(x)] = x
	y2[names(y)] = y

	sum(pmin(x2, y2))/sum(pmax(x2, y2))
}


simplifyPathway = function(mat, method = "binary_cut", control = list(), 
	plot = TRUE, term = NULL, verbose = TRUE, 
	column_title = qq("@{nrow(mat)} GO terms clustered by '@{method}'"),
	ht_list = NULL, ...) {

	if(is.atomic(mat) && !is.matrix(mat)) {
		pathway_id = mat
		mat = pathway_similarity(pathway_id)
	}
	
	cl = do.call(cluster_terms, list(mat = mat, method = method, verbose = verbose, control = control))
	pathway_id = rownames(mat)

	db = guess_pathway_database(pathway_id)

	if(db == "kegg") {
		term = get_gene_desc_from_KEGG(term_id)
	} else if(db == "reactome") {
		term = get_gene_desc_from_Reactome(term_id)
	} else if(db == "panther") {
		term = get_gene_desc_from_PANTHER(term_id)
	} else if(db == "pathbank") {
		term = get_gene_desc_from_PathBank(term_id)
	}
	term_type = "gene_description"

	if(plot) ht_clusters(mat, cl, term = term, term_type = term_type, column_title = column_title, ht_list = ht_list, ...)

	return(invisible(data.frame(id = go_id, term = term, cluster = cl, stringsAsFactors = FALSE)))
}

import_pathway_common = function(db) {

	db = tolower(db)

	f1 = qq("PathwayCommons12.@{db}.hgnc.txt.gz")
	if(!file.exists(f1)) {
		download.file(qq("https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.@{db}.hgnc.txt.gz"), f1)
	}

	f2 = qq("PathwayCommons12.@{db}.hgnc.gmt.gz")
	if(!file.exists(f2)) {
		download.file(qq("https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.@{db}.hgnc.gmt.gz"), f2)
	}

	
	qqcat("reading @{f1}...\n")	
	pathway_common = read.table(pipe(qq("gzip -d -c @{f1} | awk 'BEGIN{FS=\"\\t\";OFS=\"\\t\"} {if($6!=\"\") print $1,$2,$3,$4,$6}'")), header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char="")
	
	# http://www.pathwaycommons.org/pc2/formats#sif_relations
	l = pathway_common[, 2] %in% c("controls-state-change-of",
		                           "controls-transport-of",
		                           "controls-phosphorylation-of",
		                           "controls-expression-of",
		                           "catalysis-precedes",
		                           "in-complex-with",
		                           "interacts-with",
		                           "neighbor-of",
		                           "consumption-controled-by",
		                           "controls-production-of",
		                           "controls-transport-of-chemical",
		                           "chemical-affects",
		                           "reacts-with",
		                           "used-to-produce")
	pathway_common = pathway_common[l, ]
	if(nrow(pathway_common) == 0) {
		return(NULL)
	}

	# for the following interaction type:
	#   in-complex-with	
	#   interacts-with	
	#   neighbor-of	
	#   reacts-with	
	# we need to assign both directions
	l = pathway_common[, 2] %in% c("in-complex-with",
		                           "interacts-with",
		                           "neighbor-of",
		                           "reacts-with")
	if(any(l)) {
		qqcat("  assign directions to @{sum(l)} interactions.\n")
		pathway_common = rbind(pathway_common, pathway_common[l, c(3, 2, 1, 4, 5)])
	}

	all_pathway = strsplit(pathway_common[, "PATHWAY_NAMES"], ";")
	unique_pathway = unique(unlist(all_pathway))
	qqcat("  there are @{length(unique_pathway)} pathways.\n")
	pathway_list = lapply(seq_along(unique_pathway), function(i) {
		pa = unique_pathway[i]
		l = sapply(all_pathway, function(x) any(x == pa))
		qqcat("  [@{db} @{i}/@{length(unique_pathway)}] @{pa}: @{sum(l)} interactions.\n")
		tb = pathway_common[l, ]
		graph.edgelist(as.matrix(tb[, c(1, 3)]))
	})
	names(pathway_list) = unique_pathway

	ln = readLines(gzfile(f2))
	ln = strsplit(ln, "\t")
	map = structure(
		sapply(ln, function(x) gsub("^.*/(.*)$", "\\1", x[1])),
		names = sapply(ln, function(x) gsub("^name: (.*?);.*$", "\\1", x[2]))
	)

	cn = intersect(unique_pathway, names(map))
	pathway_list = pathway_list[cn]
	map = map[cn]

	return(list(pathway_list = pathway_list, map = map))
}

