# == title
# Calculate Disease Ontology (DO) semantic similarity matrix
#
# == param
# -do_id A vector of DO IDs.
# -measure Semantic measure for the DO similarity, pass to `DOSE::doSim`.
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
DO_similarity = function(do_id, measure = "Rel") {

	check_pkg("DOSE", bioc = FALSE)

	if(!exists(".DOSEEnv", envir = .GlobalEnv)) {
		message_wrap("'DOSE' package requires a '.DOSEEnv' variable stored in the '.GlobalEnv' environment. Please manualy create one by `.GlobalEnv$.DOSEEnv = new.env()` or simply running `library(DOSE)`.")
	}

	do_sim = DOSE::doSim(do_id, do_id, measure = measure)
	do_sim[is.na(do_sim)] = 0

	do_sim[lower.tri(do_sim)]  = t(do_sim)[lower.tri(do_sim)]

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


# # == title
# # Calculate MeSH similarity matrix
# #
# # == param
# # -mesh_id A vector of MeSH IDs.
# # -measure Measurement for the MeSH similarity, pass to `meshes::meshSim`.
# #
# # == details
# # This function is basically a wrapper on `meshes::meshSim`.
# #
# # == value
# # A symmetric matrix.
# #
# # == example
# # \donttest{
# # mesh_id = random_mesh(10)
# # mesh_similarity(mesh_id)
# # }
# mesh_similarity = function(mesh_id, measure = "Rel") {
# 	if(!requireNamespace("meshes")) {
# 		stop_wrap("'meshes' package should be installed.")
# 	}

# 	mesh_sim = meshes::meshSim(mesh_id, mesh_id, measure = measure)
# 	mesh_sim[is.na(mesh_sim)] = 0

# 	mesh_sim[lower.tri(mesh_sim)]  = t(mesh_sim)[lower.tri(mesh_sim)]
# 	diag(mesh_sim) = 1

# 	attr(mesh_sim, "measure") = measure
# 	attr(mesh_sim, "ontology") = "mesh"
# 	return(mesh_sim)
# }

# # == title
# # Generate random MeSH IDs
# #
# # == param
# # -n Number of MeSH IDs.
# #
# # == value
# # A vector of MeSH IDs.
# #
# # == example
# # \donttest{
# # random_mesh(10)
# # }
# random_mesh = function(n) {
# 	if(!requireNamespace("MeSH.db")) {
# 		stop_wrap("'MeSH.db' package should be installed.")
# 	}
# 	all_mesh_id = keys(MeSH.db::MeSH.db, keytype = "MESHID")
# 	sample(all_mesh_id, min(n, length(all_mesh_id)))
# }
