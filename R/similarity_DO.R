# == title
# Calculate Disease Ontology (DO) semantic similarity matrix
#
# == param
# -do_id A vector of DO IDs.
# -measure Semantic measure for the DO similarity, pass to `DOSE::doSim`.
# -remove_orphan_terms Whether to remove terms that have zero similarity to all other terms?
#
# == details
# This function is basically a wrapper on `DOSE::doSim`.
#
# == value
# A symmetric matrix.
#
# == example
# \donttest{
# require(DOSE)
# do_id = random_DO(10)
# DO_similarity(do_id)
# }
DO_similarity = function(do_id, measure = "Rel", remove_orphan_terms = FALSE) {

	check_pkg("DOSE", bioc = FALSE)

	if(!exists(".DOSEEnv", envir = .GlobalEnv)) {
		message_wrap("'DOSE' package requires a '.DOSEEnv' variable stored in the '.GlobalEnv' environment. Please manualy create one by `.GlobalEnv$.DOSEEnv = new.env()` or simply running `library(DOSE)`.")
	}

	do_sim = DOSE::doSim(do_id, do_id, measure = measure)
	do_sim[is.na(do_sim)] = 0

	do_sim[lower.tri(do_sim)]  = t(do_sim)[lower.tri(do_sim)]

	if(remove_orphan_terms) {
		do_sim_tmp = do_sim
		diag(do_sim_tmp) = 0
		l = rowSums(do_sim_tmp) == 0
		if(any(l)) {
			message(qq("@{sum(l)} DO term@{ifelse(sum(l) == 1, ' is', 's are')} removed because @{ifelse(sum(l) == 1, 'it has', 'they have')} zero similarity to all other terms."))
			do_sim = do_sim[l, l, drop = FALSE]
		}
	}

	attr(do_sim, "measure") = measure
	attr(do_sim, "ontology") = "DO"
	return(do_sim)
}

# == title
# Generate random Disease Ontology (DO) IDs
#
# == param
# -n Number of DO IDs.
#
# == details
# ``DO.db`` package should be installed.
#
# == value
# A vector of DO IDs.
#
# == example
# \donttest{
# random_DO(100)
# }
random_DO = function(n) {

	check_pkg("DO.db", bioc = TRUE)

	all_do_id = names(as.list(DO.db::DOTERM))
	sample(all_do_id, min(n, length(all_do_id)))
}

