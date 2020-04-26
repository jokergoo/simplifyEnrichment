
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
# -measure Measurement for the GO similarity, pass to `GOSemSim::mgoSim`.
#
# == details
# This function is basically a wrapper on `GOSemSim::mgoSim`.
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
	go_sim = mgoSim(go_id, go_id, measure = measure, semData = semData, combine = NULL)
	go_sim[is.na(go_sim)] = 0

	go_sim[lower.tri(go_sim)]  = t(go_sim)[lower.tri(go_sim)]
	diag(go_sim) = 1

	attr(go_sim, "measure") = measure
	return(go_sim)
}

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
