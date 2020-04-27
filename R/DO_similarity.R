# DO_similarity = function(do_id, measure = "Rel") {

# 	do_sim = DOSE::doSim(do_id, do_id, measure = measure)
# 	do_sim[is.na(do_sim)] = 0

# 	do_sim[lower.tri(do_sim)]  = t(do_sim)[lower.tri(do_sim)]
# 	diag(do_sim) = 1

# 	attr(do_sim, "measure") = measure
# 	attr(do_sim, "ontology") = "DO"
# 	return(do_sim)
# }


# random_DO = function(n) {
# 	all_do_id = names(as.list(DO.db::DOTERM))
# 	sample(all_do_id, min(n, length(all_do_id)))
# }
