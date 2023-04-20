

get_gene_desc_from_GO = function(go_id, ont = "BP") {
	n = length(go_id)
	desc = character(n)
	names(desc) = go_id

	GO_DATA = getFromNamespace("get_GO_data", ns = "clusterProfiler")("org.Hs.eg.db", ont, "ENTREZID")
	gl = GO_DATA$PATHID2EXTID[which(names(GO_DATA$PATHID2EXTID) %in% go_id)]

	# tb = readRDS(system.file("extdata", "refseq_gene_desc_human.rds", package = "simplifyEnrichment"))
	
	map = structure(tb[, 2], names = tb[, 1])
	tb = readRDS("~/project/development/simplifyEnrichment/inst/extdata/refseq_gene_desc_human.rds")
	
	desc[names(gl)] = sapply(gl, function(x) paste(map[as.character(x)], collapse = " "))

}

get_gene_desc_from_KEGG = function(pathway_id) {

	n = length(pathway_id)
	desc = character(n)
	names(desc) = pathway_id

	species = gsub("^([a-zA-Z]+)(\\d+$)", "\\1", pathway_id[1])
	oe = try(KEGG_DATA <- getFromNamespace("prepare_KEGG", "clusterProfiler")(species, "KEGG", "kegg"), silent = TRUE)
	if(inherits(oe, "try-error")) {
		KEGG_DATA = getFromNamespace("get_data_from_KEGG_db", "clusterProfiler")(species)
	}

	gl = KEGG_DATA$PATHID2EXTID

	# tb = readRDS(system.file("extdata", "refseq_gene_desc_human.rds", package = "simplifyEnrichment"))
	tb = readRDS("~/project/development/simplifyEnrichment/inst/extdata/refseq_gene_desc_human.rds")
	map = structure(tb[, 2], names = tb[, 1])
	
	desc[names(gl)] = sapply(gl, function(x) paste(map[as.character(x)], collapse = " "))
}

get_gene_desc_from_Reactome = function(pathway_id) {
	n = length(pathway_id)
	desc = character(n)
	names(desc) = pathway_id

	all = as.list(reactome.db::reactomePATHID2EXTID)
	gl = all[pathway_id]

	# tb = readRDS(system.file("extdata", "refseq_gene_desc_human.rds", package = "simplifyEnrichment"))
	tb = readRDS("~/project/development/simplifyEnrichment/inst/extdata/refseq_gene_desc_human.rds")
	map = structure(tb[, 2], names = tb[, 1])
	
	desc[names(gl)] = sapply(gl, function(x) paste(map[as.character(x)], collapse = " "))
}

