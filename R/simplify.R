

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
# == seealso
# `simplifyGOFromMultipleLists` which performs simplifyGO analysis with multiple lists of GO IDs.
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
	plot = TRUE, term = "TERM", verbose = TRUE, 
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
		stop_wrap("Please ensure GO IDs are the row names of the similarity matrix and should be matched to '^GO:\\\\d+$'.")
	}

	if(is.null(term)) term = "term"

	if(length(term) == 1) {
		term = tolower(term)

		if(term == "term") {
			suppressMessages(term <- select(GO.db::GO.db, keys = go_id, columns = "TERM")$TERM)
			term_type = "term"
		} else if(term == "definition") {
			suppressMessages(term <- select(GO.db::GO.db, keys = go_id, columns = "DEFINITION")$DEFINITION)
			term_type = "definition"
		} else if(term == "gene_description" || term == "gene description") {
			ont = guess_ont(go_id)
			term = get_gene_desc_from_GO(go_id, ont)
			term_type = "gene_description"
		}
	} else {
		term_type = "text"
	}

	if(plot) ht_clusters(mat, cl, term = term, term_type = term_type, column_title = column_title, ht_list = ht_list, ...)

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

	term_type = ""
	if(is.null(term)) {
		if(is_GO_id(term_id)) {
			term = get_gene_desc_from_GO(term_id)
			term_type = "gene_description"
		} else {
			db = guess_pathway_database(term_id)
			if(!is.null(db)) {
				if(db == "kegg") {
					term = get_gene_desc_from_KEGG(term_id)
					term_type = "gene_description"
				} else if(db == "reactome") {
					term = get_gene_desc_from_Reactome(term_id)
					term_type = "gene_description"
				} else if(db == "panther") {
					term = get_gene_desc_from_PANTHER(term_id)
					term_type = "gene_description"
				} else if(db == "pathbank") {
					term = get_gene_desc_from_PathBank(term_id)
					term_type = "gene_description"
				}
			}
		}
	}
	
	if(plot) ht_clusters(mat, cl, term = term, term_type = term_type, column_title = column_title, ht_list = ht_list, ...)

	if(is.null(term)) {
		return(data.frame(id = term_id, cluster = cl, stringsAsFactors = FALSE))
	} else {
		return(data.frame(id = term_id, term = term, cluster = cl, stringsAsFactors = FALSE))
	}
}

# == title
# Perform simplifyGO analysis with multiple lists of GO IDs
#
# == param
# -lt A data frame, a list of numeric vectors (e.g. adjusted p-values) where each numeric vector has GO IDs as names, or a list of GO IDs.
# -go_id_column Column index of GO ID if ``lt`` contains a list of data frames.
# -padj_column Column index of adjusted p-values if ``lt`` contains a list of data frames.
# -padj_cutoff Cut off for adjusted p-values
# -filter A self-defined function for filtering GO IDs. By default it requires GO IDs should be significant in at least one list.
# -default The default value for the adjusted p-values. See Details.
# -ont GO ontology. Value should be one of "BP", "CC" or "MF". If it is not specified,
#      the function automatically identifies it by random sampling 10 IDs from ``go_id`` (see `guess_ont`).
# -db Annotation database. It should be from https://bioconductor.org/packages/3.10/BiocViews.html#___OrgDb
# -measure Semantic measure for the GO similarity, pass to `GOSemSim::termSim`.
# -heatmap_param Parameters for controlling the heatmap, see Details.
# -method Pass to `simplifyGO`.
# -control Pass to `simplifyGO`.
# -min_term Pass to `simplifyGO`.
# -verbose Pass to `simplifyGO`.
# -column_title Pass to `simplifyGO`.
# -... Pass to `simplifyGO`.
#
# == Details
# The input data can have three types of formats:
#
# - A list of numeric vectors of adjusted p-values where each vector has the GO IDs as names.
# - A data frame. The column of the GO IDs can be specified with ``go_id_column`` argument and the column of the adjusted p-values can be
#      specified with ``padj_column`` argument. If these columns are not specified, they are automatically identified. The GO ID column
#      is found by checking whether a column contains all GO IDs. The adjusted p-value column is found by comparing the column names of the 
#      data frame to see whether it might be a column for adjusted p-values. These two columns are used to construct a numeric vector
#      with GO IDs as names.
# - A list of character vectors of GO IDs. In this case, each character vector is changed to a numeric vector where
#   all values take 1 and the original GO IDs are used as names of the vector.
#
# Now let's assume there are ``n`` GO lists, we first construct a global matrix where columns correspond to the ``n`` GO lists and rows correspond
# to the "union" of all GO IDs in the lists. The value for the ith GO ID and in the jth list are taken from the corresponding numeric vector
# in ``lt``. If the jth vector in ``lt`` does not contain the ith GO ID, the value defined by ``default`` argument is taken there (e.g. in most cases the numeric
# values are adjusted p-values, ``default`` is set to 1). Let's call this matrix as ``M0``.
#
# Next step is to filter ``M0`` so that we only take a subset of GO IDs of interest. We define a proper function via argument ``filter`` to remove
# GO IDs that are not important for the analysis. Functions for ``filter`` is applied to every row in ``M0`` and ``filter`` function needs
# to return a logical value to decide whether to remove the current GO ID. For example, if the values in ``lt`` are adjusted p-values, the ``filter`` function
# can be set as ``function(x) any(x < padj_cutoff)`` so that the GO ID is kept as long as it is signfiicant in at least one list. After the filter, let's call
# the filtered matrix ``M1``.
#
# GO IDs in ``M1`` (row names of ``M1``) are used for clustering. A heatmap of ``M1`` is attached to the left of the GO similarity heatmap so that
# the group-specific (or list-specific) patterns can be easily observed and to corresponded to GO functions.
#
# Argument ``heatmap_param`` controls several parameters for heatmap ``M1``:
#
# - ``transform``: A self-defined function to transform the data for heatmap visualization. The most typical case is to transform adjusted p-values by ``-log10(x)``.
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
# simplifyGOFromMultipleLists(lt, padj_cutoff = 0.001)
#
# # a list of numeric values
# lt2 = lapply(lt, function(x) structure(x$p.adjust, names = x$ID))
# simplifyGOFromMultipleLists(lt2, padj_cutoff = 0.001)
#
# # a list of GO IDS
# lt3 = lapply(lt, function(x) x$ID[x$p.adjust < 0.001])
# simplifyGOFromMultipleLists(lt3)
# }
simplifyGOFromMultipleLists = function(lt, go_id_column = NULL, padj_column = NULL, padj_cutoff = 1e-2,
	filter = function(x) any(x < padj_cutoff), default = 1, 
	ont = NULL, db = 'org.Hs.eg.db', measure = "Rel",
	heatmap_param = list(NULL), 
	method = "binary_cut", control = list(partial = TRUE), 
	min_term = NULL, verbose = TRUE, column_title = NULL, ...) {

	n = length(lt)

	if(is.data.frame(lt[[1]])) {

		if(is.null(go_id_column)) {
			go_id_column = which(sapply(lt[[1]], function(x) all(grepl("^GO:\\d+$", x))))[1]
			if(length(go_id_column) == 0) {
				if(!is.null(rownames(lt[[1]]))) {
					go_id_column = rownames
					if(is.null(rownames(lt[[1]]))) {
						stop_wrap("Cannot find the GO ID column in the data frames. Please explicitly set argument `go_id_column`.")
					}
					if(verbose) {
						qqcat("Use row names of the data frame as `go_id_column`.\n")
					}
				} else {
					stop_wrap("Cannot find the GO ID column in the data frames. Please explicitly set argument `go_id_column`.")
				}
			} else {
				if(verbose) {
					qqcat("Use column '@{colnames(lt[[1]])[go_id_column]}' as `go_id_column`.\n")
				}
			}
		}
		if(is.null(padj_column)) {
			cn = colnames(lt[[1]])
			ind = test_padj_column(cn)
			if(length(ind)) {
				padj_column = ind
				if(verbose) {
					qqcat("Use column '@{colnames(lt[[1]])[padj_column]}' as `padj_column`.\n")
				}
			} else {
				stop_wrap("Cannot find the column the contains adjusted p-values in the data frames. Please explicitly set argument `padj_column`.")
			}
		}

		lt = lapply(lt, function(x) {
			if(is.function(go_id_column)) {
				structure(x[, padj_column], names = go_id_column(x))
			} else {
				structure(x[, padj_column], names = x[, go_id_column])
			}
		})
		return(simplifyGOFromMultipleLists(lt, padj_cutoff = padj_cutoff, filter = filter, default = default, heatmap_param = heatmap_param, method = method, 
			control = control, min_term = min_term, verbose = verbose, column_title = column_title, ...))
		
	} else if(is.character(lt[[1]])) {
		lt = lapply(lt, function(x) structure(rep(1, length(x)), names = x))
		return(simplifyGOFromMultipleLists(lt, default = 0, filter = function(x) TRUE,
			heatmap_param = list(transform = function(x) x, breaks = c(0, 1), col = c("transparent", "red"), name = "", labels = c("not available", "available")), ...))
	}

	heatmap_param2 = list(transform = NULL, 
		breaks = NULL, col = NULL, labels = NULL, name = "padj"
	)
	for(nm in names(heatmap_param)) {
		heatmap_param2[[nm]] = heatmap_param[[nm]]
	}

	transform = heatmap_param2$transform
	if(is.null(transform)) transform = function(x) -log10(x)
	breaks = heatmap_param2$breaks
	col = heatmap_param2$col
	labels = heatmap_param2$labels
	name = heatmap_param2$name
	if(is.null(name)) name = ""

	if(is.null(breaks) && is.null(col)) {
		digit = ceiling(-log10(padj_cutoff))
		base = padj_cutoff*10^digit
		breaks = c(1, padj_cutoff, base*10^(-digit*2))
		col = c("green", "white", "red")
		labels = gt_render(c("1", qq("@{base}x10<sup>-@{digit}</sup>"), qq("@{base}x10<sup>-@{digit*2}</sup>")))
	} else if(!is.null(breaks) && !is.null(col)) {
		if(length(breaks) != length(col)) {
			stop_wrap("Length of `breaks` must be the same as the length of `col`.")
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
	sim_mat = GO_similarity(all_go_id, ont = ont, db = db, measure = measure)

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

test_padj_column = function(cn) {
	test_cn = c("p.adjust", "p_adjust", "padjust", "padj", "fdr", "FDR", "BH", "p.value", "p-value", "pvalue", "p_value")
	for(x in test_cn) {
		ind = which(cn %in% x)
		if(length(ind)) {
			return(ind[1])
		}
	}
	return(NULL)
}
