
# == title
# Visualize the similarity matrix and the clustering
#
# == param
# -mat A similarity matrix.
# -cl Cluster labels inferred from the similarity matrix, e.g. from `cluster_terms` or `binary_cut`.
# -dend Used internally.
# -col A vector of colors that map from 0 to the 95^th percentile of the similarity values.
# -draw_word_cloud Whether to draw the word clouds.
# -term The full name or the description of the corresponding GO IDs. 
# -min_term Minimal number of functional terms in a cluster. All the clusters
#     with size less than ``min_term`` are all merged into one separated cluster in the heatmap.
# -order_by_size Whether to reorder clusters by their sizes. The cluster
#      that is merged from small clusters (size < ``min_term``) is always put to the bottom of the heatmap.
# -exclude_words Words that are excluded in the word cloud.
# -max_words Maximal number of words visualized in the word cloud.
# -word_cloud_grob_param A list of graphic parameters passed to `word_cloud_grob`.
# -fontsize_range The range of the font size. The value should be a numeric vector with length two.
#       The minimal font size is mapped to word frequency value of 1 and the maximal font size is mapped
#       to the maximal word frequency. The font size interlopation is linear.
# -bg_gp Graphics parameters for controlling word cloud annotation background.
# -column_title Column title for the heatmap.
# -ht_list A list of additional heatmaps added to the left of the similarity heatmap.
# -use_raster Whether to write the heatmap as a raster image.
# -run_draw Internally used.
# -... Other arguments passed to `ComplexHeatmap::draw,HeatmapList-method`.
#
# == value
# A `ComplexHeatmap::HeatmapList-class` object.
#
# == example
# mat = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds",
#     package = "simplifyEnrichment"))
# cl = binary_cut(mat)
# ht_clusters(mat, cl, word_cloud_grob_param = list(max_width = 80))
# ht_clusters(mat, cl, word_cloud_grob_param = list(max_width = 80),
#     order_by_size = TRUE)
ht_clusters = function(
	mat, 
	cl, 
	dend = NULL, 
	col = c("white", "red"),

	# arguments that control the word cloud annotation
	draw_word_cloud = is_GO_id(rownames(mat)[1]) || !is.null(term), 
	term = NULL, 
	min_term = round(nrow(mat)*0.01), 
	order_by_size = FALSE, 
	exclude_words = character(0), 
	max_words = 10,
	word_cloud_grob_param = list(), 
	fontsize_range = c(4, 16), 
	bg_gp = gpar(fill = "#DDDDDD", col = "#AAAAAA"),

	# arguments that control the heatmaps
	column_title = NULL, 
	ht_list = NULL, 
	use_raster = TRUE, 
	run_draw = TRUE,
	...) {

	if(length(col) == 1) col = c("white", rgb(t(col2rgb(col)), maxColorValue = 255))
	col_fun = colorRamp2(seq(0, quantile(mat, 0.95), length = length(col)), col)
	if(!is.null(dend)) {
		ht = Heatmap(mat, col = col_fun,
			name = "Similarity", column_title = column_title,
			show_row_names = FALSE, show_column_names = FALSE,
			cluster_rows = dend, cluster_columns = dend, 
			show_row_dend = TRUE, show_column_dend = FALSE,
			row_dend_width = unit(4, "cm"),
			border = "#404040", row_title = NULL,
			use_raster = use_raster)
		draw(ht)
		return(invisible(NULL))
	} else {
		if(inherits(cl, "try-error")) {
			grid.newpage()
			pushViewport(viewport())
			grid.text("Clustering has an error.")
			popViewport()
			return(invisible(NULL))
		}

		# if(!is.factor(cl)) cl = factor(cl, levels = unique(cl))
		cl = as.vector(cl)
		cl_tb = table(cl)
		cl[as.character(cl) %in% names(cl_tb[cl_tb < min_term])] = 0
		cl = factor(cl, levels = c(setdiff(sort(cl), 0), 0))

		if(order_by_size) {
			cl = factor(cl, levels = c(setdiff(names(sort(table(cl), decreasing = TRUE)), 0), 0))
		}
		# od2 = order.dendrogram(dend_env$dend)
		od2 = unlist(lapply(levels(cl), function(le) {
			l = cl == le
			if(sum(l) <= 1) {
				return(which(l))
			} else {
				mm = mat[l, l, drop = FALSE]
				which(l)[hclust(stats::dist(mm))$order]
			}
		}))
		ht = Heatmap(mat, col = col_fun,
			name = "Similarity", column_title = column_title,
			show_row_names = FALSE, show_column_names = FALSE,
			show_row_dend = FALSE, show_column_dend = FALSE,
			row_order = od2, column_order = od2,
			border = "#404040", row_title = NULL,
			use_raster = use_raster) + NULL

		if(is.null(term)) {
			if(is.null(rownames(mat))) {
				draw_word_cloud = FALSE
			} else if(!grepl("^GO:[0-9]+$", rownames(mat)[1])) {
				draw_word_cloud = FALSE
			}
		}

		if(draw_word_cloud) {
			go_id = rownames(mat)

			if(!is.null(term)) {
				if(length(term) != length(go_id)) {
					stop_wrap("Length of `term` should be the same as the nrow of `mat`.")
				}
			}

			align_to = split(seq_along(cl), cl)
			go_id = split(go_id, cl)
			if(!is.null(term)) term = split(term, cl)

			align_to = align_to[names(align_to) != "0"]
			go_id = go_id[names(go_id) != "0"]
			if(!is.null(term)) term = term[names(term) != 0]

			if(length(align_to)) {
				ht = ht + rowAnnotation(keywords = anno_word_cloud_from_GO(align_to, go_id, term,
					exclude_words = exclude_words, max_words = max_words, word_cloud_grob_param = word_cloud_grob_param, 
					fontsize_range = fontsize_range, bg_gp = bg_gp))
			} else {
				ht = ht + Heatmap(ifelse(cl == "0", "< 5", ">= 5"), col = c("< 5" = "darkgreen", ">= 5" = "white"), width = unit(1, "mm"),
					heatmap_legend_param = list(title = "", at = "< 5", labels = "Small clusters"),
					show_column_names = FALSE)
			}
		} else {
			if(any(cl == "0")) {
				ht = ht + Heatmap(ifelse(cl == "0", "< 5", ">= 5"), col = c("< 5" = "darkgreen", ">= 5" = "white"), width = unit(1, "mm"),
					heatmap_legend_param = list(title = "", at = "< 5", labels = "Small clusters"),
					show_column_names = FALSE)
			}
		}
	}

	gap = unit(2, "pt")
	if(!is.null(ht_list)) {
		n = length(ht_list)
		ht = ht_list + ht
		gap = unit.c(unit(rep(2, n), "mm"), gap)
	}
	
	ht@ht_list[[1]]@heatmap_param$post_fun = function(ht) {

		decorate_heatmap_body("Similarity", {
			grid.rect(gp = gpar(fill = NA, col = "#404040"))
			cl = factor(cl, levels = unique(cl[od2]))
			tbcl = table(cl)
			ncl = length(cl)
			x = cumsum(c(0, tbcl))/ncl
			grid.segments(x, 0, x, 1, default.units = "npc", gp = gpar(col = "#404040"))
			grid.segments(0, 1 - x, 1, 1 - x, default.units = "npc", gp = gpar(col = "#404040"))
		})
	}

	if(run_draw) {
		ht = draw(ht, main_heatmap = "Similarity", gap = gap, ...)
	}
	return(invisible(ht))
}

# == title
# Scale font size
#
# == param
# -x A numeric vector.
# -rg The range.
# -fs Range of the font size.
#
# == detaisl
# It is a linear interpolation.
#
# == value
# A numeric vector.
#
# == example
# x = runif(10, min = 1, max = 20)
# # scale x to fontsize 4 to 16.
# scale_fontsize(x)
scale_fontsize = function(x, rg = c(1, 30), fs = c(4, 16)) {
	k = (fs[2] - fs[1])/(rg[2] - rg[1]) 
	b = fs[2] - k*rg[2]
	y = k*x + b
	y[y < fs[1]] = fs[1]
	y[y > fs[2]] = fs[2]
	round(y)
}

