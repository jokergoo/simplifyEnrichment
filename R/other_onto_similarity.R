# == title
# Calculate Disease Ontology (DO) similarity matrix
#
# == param
# -do_id A vector of DO IDs.
# -measure Measurement for the DO similarity, pass to `GOSemSim::termSim`.
#
# == details
# This function is basically a wrapper on `GOSemSim::termSim`.
#
# == value
# A symmetric matrix.
#
DO_similarity = function(do_id, measure = "Rel") {

	if(!requireNamespace("DOSE")) {
		stop_wrap("'DOSE' package should be installed.")
	}

	if(!exist(".DOSEEnv", envir = .GlobalEnv)) {
		message_wrap("'DOSE' package requires a '.DOSEEnv' variable stored in the '.GlobalEnv' environment. Please manualy create one by `.GlobalEnv$.DOSEEnv = new.env()` or simply `library(DOSE)`.")
	}

	do_sim = DOSE::doSim(do_id, do_id, measure = measure)
	do_sim[is.na(do_sim)] = 0

	do_sim[lower.tri(do_sim)]  = t(do_sim)[lower.tri(do_sim)]
	diag(do_sim) = 1

	attr(do_sim, "measure") = measure
	attr(do_sim, "ontology") = "DO"
	return(do_sim)
}

# == title
# Generate random DO IDs
#
# == param
# -n Number of DO IDs.
#
# == value
# A vector of DO IDs.
#
random_DO = function(n) {
	if(!requireNamespace("DO.db")) {
		stop_wrap("'DO.db' package should be installed.")
	}
	all_do_id = names(as.list(DO.db::DOTERM))
	sample(all_do_id, min(n, length(all_do_id)))
}


# == title
# Calculate mesh similarity matrix
#
# == param
# -mesh_id A vector of mesh IDs.
# -measure Measurement for the mesh similarity, pass to `GOSemSim::termSim`.
#
# == details
# This function is basically a wrapper on `GOSemSim::termSim`.
#
# == value
# A symmetric matrix.
#
mesh_similarity = function(mesh_id, measure = "Rel") {
	if(!requireNamespace("meshes")) {
		stop_wrap("'meshes' package should be installed.")
	}

	mesh_sim = meshes::meshSim(mesh_id, mesh_id, measure = measure)
	mesh_sim[is.na(mesh_sim)] = 0

	mesh_sim[lower.tri(mesh_sim)]  = t(mesh_sim)[lower.tri(mesh_sim)]
	diag(mesh_sim) = 1

	attr(mesh_sim, "measure") = measure
	attr(mesh_sim, "ontology") = "mesh"
	return(mesh_sim)
}

# == title
# Generate random mesh IDs
#
# == param
# -n Number of mesh IDs.
#
# == value
# A vector of mesh IDs.
#
random_mesh = function(n) {
	if(!requireNamespace("MeSH.db")) {
		stop_wrap("'MeSH.db' package should be installed.")
	}
	all_mesh_id = keys(MeSH.db::MeSH.db, keytype = "MESHID")
	sample(all_mesh_id, min(n, length(all_mesh_id)))
}
