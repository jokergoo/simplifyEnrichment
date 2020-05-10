
env = new.env()
env$semData_hash = ""

# == title
# Calculate GO similarity matrix
#
# == param
# -go_id A vector of GO IDs.
# -ont GO ontology. Value should be one of "BP", "CC" or "MF". If it is not specified,
#      the function automatically identifies it by random sampling 10 IDs from ``go_id`` (see `guess_ont`).
# -db Annotation database. It should be from https://bioconductor.org/packages/3.10/BiocViews.html#___OrgDb
# -measure Measurement for the GO similarity, pass to `GOSemSim::termSim`.
#
# == details
# This function is basically a wrapper on `GOSemSim::termSim`.
#
# == value
# A symmetric matrix.
#
GO_similarity = function(go_id, ont, db = 'org.Hs.eg.db', measure = "Rel") {

	if(missing(ont)) {
		
		ont = guess_ont(go_id, db)

		if(is.null(ont)) {
			stop("You need to specify the ontology by `ont`.")
		} else {
			message(qq("You haven't provided value for `ont`, guess it as `@{ont}`."))
		}
	}

	hash = digest::digest(list(ont = ont, db = db))
	if(hash == env$semData_hash) {
		semData = env$semData
	} else {
		suppressMessages(semData <- godata(db, ont = ont))
		env$semData_hash = hash
		env$semData = semData
	}
	go_removed = setdiff(go_id, Lkeys(getFromNamespace("getAncestors", "GOSemSim")(semData@ont)))

	if(length(go_removed)) {
		message(qq("@{length(go_removed)}/@{length(go_id)} GO ID@{ifelse(length(go_removed) == 1, ' is', 's are')} removed."))
	}
	go_id = setdiff(go_id, go_removed)
	# go_sim = calc_similarity(go_id, measure = measure, semData = semData, mc.cores = mc.cores)
	go_sim = termSim(go_id, go_id, method = measure, semData = semData)
	go_sim[is.na(go_sim)] = 0

	go_sim[lower.tri(go_sim)]  = t(go_sim)[lower.tri(go_sim)]
	diag(go_sim) = 1

	attr(go_sim, "measure") = measure
	attr(go_sim, "ontology") = "GO"
	return(go_sim)
}


split_by_block = function(n, size) {
	size = min(c(n, size))
	REST = n %% size
    LARGE = n - REST
    NBLOCKS = n %/% size
    GROUP = rep(1:NBLOCKS, each = size)
    if (REST > 0) GROUP = c(GROUP, rep(NBLOCKS + 1, REST))
    split(1:n, GROUP)
}

# calc_similarity = function(go_id, measure, semData, mc.cores = 1) {

# 	n = length(go_id)
# 	SPLIT = split_by_block(n, floor(sqrt(n)))
# 	COMBS = expand.grid(1:length(SPLIT), 1:length(SPLIT))
# 	COMBS = t(apply(COMBS, 1, sort))
# 	COMBS = unique(COMBS)

# 	lt = mclapply(seq_len(nrow(COMBS)), function(i) {
# 		ind1 = SPLIT[[ COMBS[i, 1] ]]
# 		ind2 = SPLIT[[ COMBS[i, 2] ]]
# 		termSim(go_id[ind1], go_id[ind2], method = measure, semData = semData)
# 	}, mc.cores = mc.cores)

# 	m = matrix(nrow = n, ncol = n)
# 	dimnames(m) = list(go_id, go_id)
# 	dimnames(m) = list(go_id, go_id)
# 	for(i in seq_len(nrow(COMBS))) {
# 		ind1 = SPLIT[[ COMBS[i, 1] ]]
# 		ind2 = SPLIT[[ COMBS[i, 2] ]]
# 		if(COMBS[i, 1] == COMBS[i, 2]) {
# 			m[ind1, ind2] = lt[[i]]
# 		} else {
# 			m[ind1, ind2] = lt[[i]]
# 			m[ind2, ind1] = lt[[i]]
# 		}
# 	}
# 	return(m)
# }

# == title
# Guess the ontology of the input GO IDs
#
# == param
# -go_id A vector of GO IDs.
# -db Annotation database. It should be from https://bioconductor.org/packages/3.10/BiocViews.html#___OrgDb
#
# == details
# 10 GO IDs are randomly sampled and checked.
#
# == value
# A single character scalar of "BP", "CC" or "MF".
#
# If there are more than one ontologies detected. It returns ``NULL``.
guess_ont = function(go_id, db = 'org.Hs.eg.db') {
	test_go_id = sample(go_id, min(c(length(go_id), 10)))
	suppressMessages(df <- select(get(db, asNamespace(db)), keys = test_go_id, columns = "ONTOLOGY", keytype = "GO"))
	guess_ont = unique(df$ONTOLOGY)
	guess_ont = guess_ont[!is.na(guess_ont)]
	if(length(guess_ont) != 1) {
		return(NULL)
	} else {
		return(guess_ont)
	}
}

# == title
# Generate random GO IDs
#
# == param
# -n Number of GO IDs.
# -ont GO ontology. Value should be one of "BP", "CC" or "MF".
# -db Annotation database. It should be from https://bioconductor.org/packages/3.10/BiocViews.html#___OrgDb
#
# == value
# A vector of GO IDs.
#
random_GO = function(n, ont = "BP", db = 'org.Hs.eg.db') {
	hash = digest::digest(list(ont = ont, db = db))
	if(hash == env$semData_hash) {
		semData = env$semData
	} else {
		suppressMessages(semData <- godata(db, ont = ont))
		env$semData_hash = hash
		env$semData = semData
	}

	all_go_id = unique(semData@geneAnno$GO)

	sample(all_go_id, min(n, length(all_go_id)))
}
