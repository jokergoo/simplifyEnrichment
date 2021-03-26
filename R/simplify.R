

# == title
# Simplify Gene Ontology (GO) enrichment results
#
# == param
# -mat A GO similarity matrix.
# -method Method for clustering the matrix. See `cluster_terms`.
# -control A list of parameters for controlling the clustering method, passed to `cluster_terms`.
# -plot Whether to make the heatmap.
# -term The full name or the description of the corresponding GO IDs. The values are automatically
#      extracted if it is not provided.
# -column_title Column title for the heatmap.
# -verbose Whether to print messages.
# -ht_list A list of additional heatmaps added to the left of the similarity heatmap.
# -... Arguments passed to `ht_clusters`.
#
# == details
# This is basically a wrapper function that it first runs `cluster_terms` to cluster
# GO terms and then runs `ht_clusters` to visualize the clustering.
#
# The arguments in `simplifyGO` passed to `ht_clusters` are:
#
# -``draw_word_cloud`` Whether to draw the word clouds.
# -``min_term`` Minimal number of GO terms in a cluster. All the clusters
#     with size less than ``min_term`` are all merged into one single cluster in the heatmap.
# -``order_by_size`` Whether to reorder GO clusters by their sizes. The cluster
#      that is merged from small clusters (size < ``min_term``) is always put to the bottom of the heatmap.
# -``exclude_words`` Words that are excluded in the word cloud.
# -``max_words`` Maximal number of words visualized in the word cloud.
# -``word_cloud_grob_param`` A list of graphic parameters passed to `word_cloud_grob`.
# -``fontsize_range`` The range of the font size. The value should be a numeric vector with length two.
#       The minimal font size is mapped to word frequency value of 1 and the maximal font size is mapped
#       to the maximal word frequency. The font size interlopation is linear.
# -``bg_gp`` Graphic parameters for controlling the background of word cloud annotations.
#
# == value
# A data frame with three columns: GO IDs, GO term names and cluster labels.
#
# == example
# \donttest{
# set.seed(123)
# go_id = random_GO(500)
# mat = GO_similarity(go_id)
# df = simplifyGO(mat, word_cloud_grob_param = list(max_width = 80))
# head(df)
# }
simplifyGO = function(mat, method = "binary_cut", control = list(), 
	plot = TRUE, term = NULL, verbose = TRUE, 
	column_title = qq("@{nrow(mat)} GO terms clustered by '@{method}'"),
	ht_list = NULL, ...) {

	if(is.atomic(mat) && !is.matrix(mat)) {
		go_id = mat
		if(!all(grepl("^GO:\\d+$", go_id))) {
			stop_wrap("If you specify a vector, it should contain all valid GO IDs.")
		}

		mat = GO_similarity(go_id)
	}
	
	cl = do.call(cluster_terms, list(mat = mat, method = method, verbose = verbose, control = control))
	go_id = rownames(mat)

	if(!all(grepl("^GO:\\d+$", go_id))) {
		stop_wrap("Please ensure GO IDs are the row names of the similarity matrix and should be matched to '^GO:\\d+$'.")
	}

	if(is.null(term)) {
		suppressMessages(term <- select(GO.db::GO.db, keys = go_id, columns = "TERM")$TERM)
	}

	if(plot) ht_clusters(mat, cl, term = term, column_title = column_title, ht_list = ht_list, ...)

	return(invisible(data.frame(id = go_id, term = term, cluster = cl, stringsAsFactors = FALSE)))
}

# == title
# Simplify functional enrichment results
#
# == param
# -mat A similarity matrix.
# -method Method for clustering the matrix. See `cluster_terms`.
# -control A list of parameters for controlling the clustering method, passed to `cluster_terms`.
# -plot Whether to make the heatmap.
# -term The full name or the description of the corresponding terms. 
# -column_title Column title for the heatmap.
# -verbose Whether to print messages.
# -ht_list A list of additional heatmaps added to the left of the similarity heatmap.
# -... Arguments passed to `ht_clusters`.
#
# == details
# The usage is the same as `simplifyGO`, except you need to manually provide the term names by ``term`` argument
# if you want to draw the word clouds.
#
simplifyEnrichment = function(mat, method = "binary_cut", control = list(), 
	plot = TRUE, term = NULL, verbose = TRUE, 
	column_title = qq("@{nrow(mat)} terms clustered by '@{method}'"),
	ht_list = NULL, ...) {
	
	cl = do.call(cluster_terms, list(mat = mat, method = method, verbose = verbose, control = control))
	term_id = rownames(mat)
	if(is.null(term_id)) {
		term_id = paste0("row_", 1:nrow(mat))
	}
	
	if(plot) ht_clusters(mat, cl, term = term, column_title = column_title, ht_list = ht_list, ...)

	if(is.null(term)) {
		return(data.frame(id = term_id, cluster = cl, stringsAsFactors = FALSE))
	} else {
		return(data.frame(id = term_id, term = term, cluster = cl, stringsAsFactors = FALSE))
	}
}

# == title
# Apply simplifyGO analysis with several lists of GO IDs
#
# == param
# -lt Preferable a list of numeric vectors where each numeric vector has GO IDs as names. It also accepts other format, see Details.
# -filter A self-defined function for filtering GO IDs, see Details.
# -default The default value, see Details.
# -heatmap_param Parameters that control the heatmap, see Details.
# -method Pass to `simplifyGO`.
# -control Pass to `simplifyGO`.
# -min_term Pass to `simplifyGO`.
# -verbose Pass to `simplifyGO`.
# -column_title Pass to `simplifyGO`.
# -... Pass to `simplifyGO`.
#
# == Details
# The input data is preferable to be a list of nmeric vectors where each numeric vector has GO IDs as names. Nevertheless, it also allows
# two other formats where they will be converted to a list of numeric vectors internally:
#
# - If ``lt`` is specified as a list of character vectors of GO IDs. Each character vector is changed to a numeric vector where
#   all values take 1 and the original GO IDs are used as names of the vector.
# - If ``lt`` is a list of data frames and each data frame is from clusterProfiler analysis, the column "p.adjust" in the data frame is
#   taken as the numeric vector.
#
# Now let's assume there are ``n`` GO lists, we first construct a global matrix where columns correspond to the GO lists and rows correspond
# to the "union" of all GO IDs in the lists. The value for the ith GO ID and in the jth list are taken from the corresponding numeric vector
# in ``lt``. If the jth vector in ``lt`` does not contain the ith GO ID, the value defined by ``default`` argument is taken there (e.g. if the numeric
# values are p-values or FDRs, ``default`` can be set to 1). Let's call this matrix as ``M0``.
#
# Next step is to filter ``M0`` so that we only take a subset of GO IDs of interest. We define a proper function via argument ``filter`` to remove
# GO IDs that are not important for the analysis. Functions for ``filter`` is applied to every row in ``M0`` and ``filter`` function needs
# to return a logical value to decide whether to remove the current GO ID. For example, if the values in ``lt`` are p-values or FDRs, the ``filter`` function
# can be set as ``function(x) any(x < 0.01)`` so that the GO ID is kept as long as it is signfiicant in at least one list. After the filter, let's call
# the filtered matrix ``M1``.
#
# GO IDs in ``M1`` (row names of ``M1``) are used for clustering. A heatmap of ``M1`` is attached to the left of the GO similarity heatmap so that
# the group-specific (or list-specific) patterns can be easily observed and to correspond to GO functions.
#
# Argument ``heatmap_param`` controls several parameters for heatmap ``M1``:
#
# - ``transform``: A self-defined function to transform the data for heatmap visualization. The most typical case is to transform FDRs by ``-log10(x)``.
# - ``breaks``: break values for color interpolation.
# - ``col``: The corresponding values for ``breaks``.
# - ``labels``: The corresponding labels.
# - ``name``: Legend title.
#
# == example
# \donttest{
# # perform functional enrichment on the signatures genes from cola anlaysis 
# require(cola)
# data(golub_cola) 
# res = golub_cola["ATC:skmeans"]
# require(hu6800.db)
# x = hu6800ENTREZID
# mapped_probes = mappedkeys(x)
# id_mapping = unlist(as.list(x[mapped_probes]))
# lt = functional_enrichment(res, k = 3, id_mapping = id_mapping) # you can check the value of `lt`
#
# # a list of data frames
# simplifyGOFromMultipleLists(lt)
#
# # a list of numeric values
# lt2 = lapply(lt, function(x) structure(x$p.adjust, names = x$ID))
# simplifyGOFromMultipleLists(lt2, filter = function(x) any(x <  0.001), default = 1,
#     heatmap_param = list(transform = function(x) -log10(x), name = "FDR"))
#
# # a list of GO IDS
# lt3 = lapply(lt, function(x) x$ID[x$p.adjust < 0.001])
# simplifyGOFromMultipleLists(lt3)
# }
simplifyGOFromMultipleLists = function(lt, filter = function(x) TRUE, default = NA, 
	heatmap_param = list(transform = function(x) x, 
		breaks = NULL, col = NULL, labels = NULL, name = "Value"
	), 
	method = "binary_cut", control = list(partial = TRUE), 
	min_term = NULL, verbose = TRUE, column_title = NULL, ...) {

	n = length(lt)

	if(is.data.frame(lt[[1]])) {
		if(identical(colnames(lt[[1]])[1:4], c("ID", "Description", "GeneRatio", "BgRatio"))) {
			lt = lapply(lt, function(x) structure(x$p.adjust, names = rownames(x)))
			return(simplifyGOFromMultipleLists(lt, filter = function(x) any(x < 1e-3), default = 1, 
				heatmap_param = list(transform = function(x) -log10(x), breaks = c(1, 1e-3, 1e-6), labels = gt_render(c("1", "1x10<sup>-3</sup>", "1x10<sup>-6</sup>")),
					col = c("green", "white", "red"), name = "FDR"), ...))
		}
	} else if(is.character(lt[[1]])) {
		lt = lapply(lt, function(x) structure(rep(1, length(x)), names = x))
		return(simplifyGOFromMultipleLists(lt, default = 0, 
			heatmap_param = list(breaks = c(0, 1), col = c("transparent", "red"), name = "", labels = c("not available", "available")), ...))
	}

	heatmap_param2 = list(transform = function(x) x, 
		breaks = NULL, col = NULL, labels = NULL, name = NULL
	)
	for(nm in names(heatmap_param)) {
		heatmap_param2[[nm]] = heatmap_param[[nm]]
	}

	transform = heatmap_param2$transform
	breaks = heatmap_param2$breaks
	col = heatmap_param2$col
	labels = heatmap_param2$labels
	name = heatmap_param2$name
	if(is.null(name)) name = ""

	if(is.null(breaks) && is.null(col)) {
		if(is_p_value(unlist(lt))) {
			breaks = c(1, 1e-2, 1e-4)
			transform = function(x) -log10(x)
			labels = gt_render(c("1", "1x10<sup>-2</sup>", "1x10<sup>-4</sup>"))
			col = c("green", "white", "red")
		}
	}

	all_go_id = unique(unlist(lapply(lt, names)))
	if(!all(grepl("^GO:\\d+$", all_go_id))) {
		stop_wrap("Only GO ID is allowed.")
	}

	m = matrix(default, nrow = length(all_go_id), ncol = n)
	rownames(m) = all_go_id
	colnames(m) = names(lt)
	if(is.null(colnames)) colnames = paste0("Group", 1:n)

	for(i in 1:n) {
		m[names(lt[[i]]), i] = lt[[i]]
	}

	l = apply(m, 1, function(x) {
		if(all(is.na(x))) {
			FALSE
		} else {
			l = filter(x[!is.na(x)])
			if(length(l) == 1) {
				return(l)
			} else {
				return(any(l))
			}
		}
	})
	m = m[l, , drop = FALSE]
	m = t(apply(m, 1, transform))

	if(verbose) qqcat("@{nrow(m)}/@{length(all_go_id)} GO IDs left for clustering.\n")

	if(length(unique(m[!is.na(m)])) <= 2) {
		col = structure(col, names = breaks)
	} else {
		if(is.null(breaks) && is.null(col)) {
			col = NULL
		} else if(!is.null(breaks) && !is.null(col)) {
			if(length(breaks) != length(col)) {
				stop_wrap("Length of `breaks` and `col` should be the same.")
			}
			col = colorRamp2(transform(breaks), col)
		} else {
			stop_wrap("Arguments `breaks` and `col` should be set at the same time.")
		}
	}

	heatmap_legend_param = list()
	heatmap_legend_param$at = transform(breaks)
	heatmap_legend_param$labels = if(is.null(labels)) breaks else labels
	heatmap_legend_param$title = name
	ht = Heatmap(m, col = col, name = if(name == "") NULL else name,
		show_row_names = FALSE, cluster_columns = FALSE,
		border = "black",
		heatmap_legend_param = heatmap_legend_param,
		width = unit(0.5, "cm")*n, use_raster = TRUE)

	all_go_id = rownames(m)
	sim_mat = GO_similarity(all_go_id)

	if(is.null(min_term)) min_term = round(nrow(sim_mat)*0.02)
	if(is.null(column_title)) column_title = qq("@{length(all_go_id)} GO terms clustered by '@{method}'")
	
	simplifyGO(sim_mat, ht_list = ht, method = method, 
		verbose = verbose, min_term = min_term, control= control, column_title = column_title, ...)
}

is_p_value = function(x) {
	if(!all(x <= 1 & x >= 0)) {
		return(FALSE)
	}
	v = -log10(x)
	if(sum(v > 2)/length(v) > 0.05) {
		TRUE
	} else {
		FALSE
	}
}
