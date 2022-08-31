
# == title
# Calculate word frequency
#
# == param
# -term A vector of description texts.
# -exclude_words The words that should be excluded.
# -stop_words The stop words that should be be removed.
# -min_word_length Minimum length of the word to be counted.
# -tokenizer The tokenizer function, one of the values accepted by ``tm::termFreq``.
# -transform_case The function normalizing lettercase of the words.
# -remove_numbers Whether to remove numbers.
# -remove_punctuation Whether to remove punctuation.
# -custom_transformer Custom function that transforms words.
# -stemming Whether to only keep the roots of inflected words.
# -dictionary A vector of words to be counted (if given all other words will be excluded).
#
# == details
# The text preprocessing followings the instructions from http://www.sthda.com/english/wiki/word-cloud-generator-in-r-one-killer-function-to-do-everything-you-need .
#
# == value
# A data frame with words and frequencies.
#
# == example
# gm = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds", package = "simplifyEnrichment"))
# go_id = rownames(gm)
# go_term = AnnotationDbi::select(GO.db::GO.db, keys = go_id, columns = "TERM")$TERM
# count_words(go_term)
count_words = function(term,
	exclude_words = NULL, stop_words = stopwords(),
	min_word_length = 1, tokenizer = 'words', transform_case = tolower,
	remove_numbers = TRUE, remove_punctuation = TRUE, custom_transformer = NULL,
	stemming = FALSE, dictionary = NULL) {
	
	# http://www.sthda.com/english/wiki/word-cloud-generator-in-r-one-killer-function-to-do-everything-you-need

	# Load the text as a corpus
	suppressWarnings({
		docs = VCorpus(VectorSource(term))
		# Convert the text to lower case
		docs = tm_map(docs, content_transformer(transform_case))
		if (remove_numbers) {
			# Remove numbers
			docs = tm_map(docs, removeNumbers)
		}
		# Remove stopwords for the language
		docs = tm_map(docs, removeWords, stop_words)
		if (remove_punctuation) {
			# Remove punctuations
			docs = tm_map(docs, removePunctuation)
		}
		# Eliminate extra white spaces
		docs = tm_map(docs, stripWhitespace)
		# Remove your own stopwords
		docs = tm_map(docs, removeWords, exclude_words)
		# Apply any user-provided transformer
		if (!is.null(custom_transformer)) {
			docs = tm_map(docs, content_transformer(custom_transformer))
		}

		# Create term-document matrix
		tdm = TermDocumentMatrix(
			docs,
			control = list(
				wordLengths = c(min_word_length, Inf),
				tokenize = tokenizer,
				stemming = stemming,
				dictionary = dictionary,
				# letter case transformation is handled earlier (above), let's not overwrite the results
				tolower = FALSE
			)
		)
	})
	v = sort(slam::row_sums(tdm), decreasing = TRUE)
	d = data.frame(word = names(v), freq = v, stringsAsFactors = FALSE)
	d
}


# == title
# A simple grob for the word cloud
#
# == param
# -text A vector of words.
# -fontsize The corresponding font size. With the frequency of the words known, `scale_fontsize` can be used to linearly interpolate frequencies to font sizes.
# -line_space Space between lines. The value can be a `grid::unit` object or a numeric scalar which is measured in mm.
# -word_space Space between words. The value can be a `grid::unit` object or a numeric scalar which is measured in mm.
# -max_width The maximal width of the viewport to put the word cloud. The value can be a `grid::unit` object or a numeric scalar which is measured in mm.
#        Note this might be larger than the final width of the returned grob object.
# -col Colors for the words. The value can be a vector, in numeric or character, which should have the same
#      length as ``text``. Or it is a self-defined function that takes the font size vector as 
#      the only argument. The function should return a color vector. See Examples.
# -add_new_line Whether to add new line after every word? If ``TRUE``, each word will be in a separated line.
# -test Internally used. It basically adds borders to the words and the viewport.
#
# == value
# A `grid::grob` object. The width and height of the grob can be get by `grid::grobWidth` and `grid::grobHeight`.
#
# == example
# # very old R versions do not have strrep() function
# if(!exists("strrep")) {
#     strrep = function(x, i) paste(rep(x, i), collapse = "")
# }
# words = sapply(1:30, function(x) strrep(sample(letters, 1), sample(3:10, 1)))
# require(grid)
# gb = word_cloud_grob(words, fontsize = runif(30, min = 5, max = 30), 
#     max_width = 100)
# grid.newpage(); grid.draw(gb)
#
# # color as a single scalar
# gb = word_cloud_grob(words, fontsize = runif(30, min = 5, max = 30), 
#     max_width = 100, col = 1)
# grid.newpage(); grid.draw(gb)
#
# # color as a vector
# gb = word_cloud_grob(words, fontsize = runif(30, min = 5, max = 30), 
#     max_width = 100, col = 1:30)
# grid.newpage(); grid.draw(gb)
#
# # color as a function
# require(circlize)
# col_fun = colorRamp2(c(5, 17, 30), c("blue", "black", "red"))
# gb = word_cloud_grob(words, fontsize = runif(30, min = 5, max = 30), 
#     max_width = 100, col = function(fs) col_fun(fs))
# grid.newpage(); grid.draw(gb)
#
word_cloud_grob = function(text, fontsize, 
	line_space = unit(4, "pt"), word_space = unit(4, "pt"), max_width = unit(80, "mm"), 
	col = function(fs) circlize::rand_color(length(fs), luminosity = "dark"),
	add_new_line = FALSE, test = FALSE) { # width in mm

	if(length(text) != length(fontsize)) {
		stop_wrap("`text` and `fontsize` should the same length.")
	}

	n = length(text)

	if(is.character(col) || is.numeric(col)) {
		if(length(col) == 1) col = rep(col, n)
	} else if(is.function(col)) {
		col = col(fontsize)
	} else {
		stop_wrap("`col` can only be a function or a character vector.")
	}
	
	od = order(fontsize, decreasing = TRUE)
	text = text[od]
	fontsize = fontsize[od]
	col = col[od]

	# if(Sys.info()["sysname"] == "Darwin" && dev.interactive()) {
	# 	ComplexHeatmap:::dev.null()
	# 	on.exit(ComplexHeatmap:::dev.off2())
	# }

	text_gb_lt = lapply(seq_len(n), function(i) textGrob(text[i], gp = gpar(fontsize = fontsize[i])))
	text_width = vapply(text_gb_lt, function(gb) convertWidth(grobWidth(gb), "mm", valueOnly = TRUE), 0)
	text_height = vapply(text_gb_lt, function(gb) convertHeight(grobHeight(gb), "mm", valueOnly = TRUE), 0)

	if(is.unit(line_space)) line_space = convertHeight(line_space, "mm", valueOnly = TRUE)
	if(is.unit(word_space)) word_space = convertWidth(word_space, "mm", valueOnly = TRUE)

	x = numeric(n)
	y = numeric(n)
	current_line_height = 0
	current_line_width = 0

	# the first text
	current_line_height = text_height[1]
	current_line_width = text_width[1]
	x[1] = 0
	y[1] = 0

	w = text_width[1]
	h = text_height[1]

	if(is.unit(max_width)) {
		max_width = convertWidth(max_width, "mm", valueOnly = TRUE)
	} 

	for(i in seq_len(n)[-1]) {
		# the next text can be put on the same line
		if(current_line_width + text_width[i] <= max_width && !add_new_line) {
			x[i] = current_line_width + word_space
			y[i] = y[i-1] # same as previous one
			current_line_width = x[i] + text_width[i]
			w = max(w, current_line_width)
			h = max(h, y[i] + text_height[i])
		} else { # the next text need to be put on the next line
			x[i] = 0
			y[i] = current_line_height + line_space
			current_line_width = text_width[i]
			current_line_height = y[i] + text_height[i]
			w = max(w, current_line_width)
			h = max(h, current_line_height)
		}
	}

	if(test) {
		gl = gList(
			rectGrob(),
			textGrob(text, x = x, y = y, gp = gpar(fontsize = fontsize, col = col), 
				default.units = "mm", just = c(0, 0)),
			rectGrob(x = x, y = y, width = text_width, height = text_height, default.units = "mm", just = c(0, 0))

		)
	} else {
		gl = gList(
			textGrob(text, x = x, y = y, gp = gpar(fontsize = fontsize, col = col), 
				default.units = "mm", just = c(0, 0))
		)
	}

	gb = gTree(children = gl, cl = "word_cloud", vp = viewport(width = unit(w, "mm"), height = unit(h, "mm")))
	return(gb)
}

# == title
# Width for word_cloud grob
#
# == param
# -x The ``word_cloud`` grob returned by `word_cloud_grob`.
#
# == value
# A `grid::unit` object.
widthDetails.word_cloud = function(x) {
	x$vp$width
}

# == title
# Height for word_cloud grob
#
# == param
# -x The ``word_cloud`` grob returned by `word_cloud_grob`.
#
# == value
# A `grid::unit` object.
heightDetails.word_cloud = function(x) {
	x$vp$height
}


# == title
# Word cloud annotations
#
# == param
# -align_to How to align the annotations to the heatmap. Similar as in `ComplexHeatmap::anno_link`, the value of ``align_to``
#           can be a list of row indices or a categorical vector where each vector in the list corresponds to a word cloud. 
#           If it is a categorical vector, rows with the same level correspond to a same word cloud. 
#           If ``align_to`` is a categorical vector and ``term`` is a list, names of ``term`` should have overlap to the levels in ``align_to``.
#           When ``align_to`` is set as a categorical vector, normally the same value is set to ``row_split`` in the main heatmap so that each row slice
#           can correspond to a word cloud.
# -term The description text used for constructing the word clouds. The value should have the same format as ``align_to``. If ``align_to``
#          is a list, ``term`` should also be a list. In this case, the length of vectors in ``term`` is not necessarily the same
#          as in ``align_to``. E.g. ``length(term[[1]])`` is not necessarily equal to ``length(align_to[[1]]``. If ``align_to``
#          is a categorical vector, ``term`` should also be a character vector with the same length as ``align_to``.
#          To make it more genrall, when ``align_to`` is a list, ``term`` can also be a list of data frames where the first column contains 
#          keywords and the second column contains numeric values that will be mapped to font sizes in the word clouds.
# -exclude_words The words excluced for construcing word cloud.
# -max_words Maximal number of words visualized in the word cloud.
# -word_cloud_grob_param A list of graphics parameters passed to `word_cloud_grob`.
# -fontsize_range The range of the font size. The value should be a numeric vector with length two.
#       The font size interpolation is linear.
# -value_range The range of values to map to font sizes.
# -bg_gp Graphics parameters for controlling the background.
# -side Side of the annotation relative to the heatmap.
# -add_new_line Whether to add new line after every word? If ``TRUE``, each word will be in a separated line.
# -count_words_param A list of parameters passed to `count_words`.
# -... Other parameters.
#
# == details
# The word cloud annotation is constructed by `ComplexHeatmap::anno_link`.
#
# If the annotation is failed to construct or no keyword is found, the function returns a `ComplexHeatmap::anno_empty` with 1px width.
#
# English stop words, punctuation and numbers are removed by default when counting words. As specific stop words might
# coincide with gene or pathway names, and numbers in genes names might be meaningful it is recommended to adjust this
# behaviour by passing appropriate arguments to the `count_words` function using ``count_words_param``.
#
# == example
# gm = readRDS(system.file("extdata", "random_GO_BP_sim_mat.rds", package = "simplifyEnrichment"))
# go_id = rownames(gm)
# go_term = AnnotationDbi::select(GO.db::GO.db, keys = go_id, columns = "TERM")$TERM
#
# split = sample(letters[1:4], 100, replace = TRUE)
# align_to = split(1:100, split)
# term = lapply(letters[1:4], function(x) sample(go_term, sample(100:400, 1)))
# names(term) = letters[1:4]
#
# require(ComplexHeatmap)
# mat = matrix(rnorm(100*10), nrow = 100)
# Heatmap(mat, cluster_rows = FALSE, row_split = split, 
# 	right_annotation = rowAnnotation(foo = anno_word_cloud(align_to, term)))
#
anno_word_cloud = function(align_to, term, exclude_words = NULL, max_words = 10,
	word_cloud_grob_param = list(), fontsize_range = c(4, 16), value_range = NULL,
	bg_gp = gpar(fill = "#DDDDDD", col = "#AAAAAA"), side = c("right", "left"),
	add_new_line = FALSE, count_words_param = list(), ...) {

	if(is.atomic(align_to) && is.list(term)) {
		align_to = split(seq_along(align_to), align_to)

		cn = intersect(names(align_to), names(term))
		if(length(cn) == 0) {
			stop_wrap("names of `term` should have overlap to the `align_to`.")
		} else {
			align_to = align_to[cn]
			term = term[cn]
		}

	} else if(is.list(align_to)) {
		if(!is.list(term)) {
			stop_wrap("`term` should have the same format as `align_to`, which is a list of term descriptions. Note, e.g. `length(term[[1]])` are not necessarily equal to `length(align_to[[1]]`.")
		} else {
			if(length(align_to) != length(term)) {
				stop_wrap("`align_to` and `term` should be two list with the same length.")
			}
		}
		if(!is.null(names(align_to)) && !is.null(names(term))) {
			if(length(setdiff(names(align_to), names(term))) == 0) {
				term = term[names(align_to)]
			}
		}
	} else if(is.atomic(align_to)) {
		if(!is.atomic(term)) {
			stop_wrap("`term` should have the same format as `align_to`, which is a vector of term descriptions with the same length as `align_to`.")
		}
		term = split(term, align_to)
		align_to = split(seq_along(align_to), align_to)
	}

	flag = 0
	if(is.list(term)) {
		is_term_a_data_frame = sapply(term, is.data.frame)
		if(any(is_term_a_data_frame)) {
			if(!all(is_term_a_data_frame)) {
				stop_wrap("Elements in `term` should all be data.frames")
			}

			sapply(term, function(x) {
				if(ncol(x) < 2) {
					stop_wrap("Data frame in `term` should have more than one column.")
				}
				if(!any(sapply(x[, -1, drop = FALSE], is.numeric))) {
					stop_wrap("There should be a numeric column in the data frame.")
				}
			})
			flag = 1
		}
	}

	if(flag) {
		keywords = lapply(term, function(df) {
			if(nrow(df) > max_words) {
				df = df[order(df[, 2], decreasing = TRUE)[seq_len(max_words)], ]
			}
			df[!df[, 1] %in% exclude_words, , drop = FALSE]
		})
	} else {
		keywords = lapply(term, function(desc) {
			combined_count_params = c(list(term = desc, exclude_words = exclude_words), count_words_param)
			suppressMessages(suppressWarnings(df <- do.call(count_words, combined_count_params)))
			# df = df[df$freq > 1, , drop = FALSE]
			if(nrow(df) > max_words) {
				df = df[order(df$freq, decreasing = TRUE)[seq_len(max_words)], ]
			}
			df
		})
	}
	keywords = keywords[vapply(keywords, nrow, 0) > 0]

	align_to = align_to[names(keywords)]

	if(length(keywords) == 0) {
		return(anno_empty(border = FALSE, width = unit(1, "pt")))
	}

	if(length(word_cloud_grob_param)) {
		if(is.atomic(word_cloud_grob_param)) {
			stop_wrap("`word_cloud_grob_param` should be a named list.")
		}
	}

	all_keywords = unique(unlist(lapply(keywords, function(x) x[, 1])))
	all_keywords_col = circlize::rand_color(length(all_keywords), luminosity = "dark")
	i_try = 1
	while(1) {
		hsv = coords(as(hex2RGB(all_keywords_col), "HSV"))
		l = hsv[, 3] > 0.85 | (hsv[, 1] > 40 & hsv[, 1] < 65)
		if(!any(l)) break
		i_try = i_try + 1
		if(i_try > 50) break
		all_keywords_col[l] = circlize::rand_color(sum(l), luminosity = "dark")
	}
	all_keywords_col = structure(all_keywords_col, names = all_keywords)

	word_cloud_grob_param = word_cloud_grob_param[setdiff(names(word_cloud_grob_param), c("text", "fontsize"))]
	ComplexHeatmap:::dev.null()
	oe = try({
		gbl <- lapply(names(align_to), function(nm) {
			kw = rev(keywords[[nm]][, 1])
			freq = rev(keywords[[nm]][, 2])
			if(is.null(value_range)) {
				fontsize = scale_fontsize(freq, rg = c(1, max(10, freq)), fs = fontsize_range)
			} else {
				fontsize = scale_fontsize(freq, rg = value_range, fs = fontsize_range)
			}

			if(!"col" %in% names(word_cloud_grob_param)) {
				word_cloud_grob_param$col = function(fs) all_keywords_col[kw]
			}

			lt = c(list(text = kw, fontsize = fontsize, add_new_line = add_new_line), word_cloud_grob_param)
			do.call(word_cloud_grob, lt)
		})
		names(gbl) = names(align_to)
	
		margin = unit(8, "pt")
		gbl_h = lapply(gbl, function(x) convertHeight(grobHeight(x), "cm") + margin)
		gbl_h = do.call(unit.c, gbl_h)

		gbl_w = lapply(gbl, function(x) convertWidth(grobWidth(x), "cm"))
		gbl_w = do.call(unit.c, gbl_w)
		gbl_w = max(gbl_w) + margin

	}, silent = TRUE)
	ComplexHeatmap:::dev.off2()
	if(inherits(oe, "try-error")) {
		stop(oe)
	}

	if(is.null(bg_gp$fill)) bg_gp$fill = "#DDDDDD"
	if(is.null(bg_gp$col)) bg_gp$col = "#AAAAAA"
	if(is.null(bg_gp$lty)) bg_gp$lty = 1
	if(is.null(bg_gp$lwd)) bg_gp$lwd = 1

	side = match.arg(side)[1]

	panel_fun = function(index, nm) {
		pushViewport(viewport())
		grid.rect(gp = gpar(fill = bg_gp$fill, col = bg_gp$fill, lty = bg_gp$lty, lwd = bg_gp$lwd))
		if(side == "right") {
			grid.lines(c(0, 1, 1, 0), c(0, 0, 1, 1), gp = gpar(col = bg_gp$col, lty = bg_gp$lty, lwd = bg_gp$lwd), default.units = "npc")
		} else {
			grid.lines(c(1, 0, 0, 1), c(0, 0, 1, 1), gp = gpar(col = bg_gp$col, lty = bg_gp$lty, lwd = bg_gp$lwd), default.units = "npc")
		}
	    pushViewport(viewport(width = unit(1, "npc") - margin, height = unit(1, "npc") - margin))
	    gb = gbl[[nm]]
	    gb$vp$x = gb$vp$width*0.5
	    gb$vp$y = gb$vp$height*0.5
	    grid.draw(gb)
	    popViewport()
	    popViewport()
	}

	anno = anno_link(align_to = align_to, which = "row", panel_fun = panel_fun, 
    	size = gbl_h, gap = unit(2, "mm"), width = gbl_w + unit(5, "mm"),
    	link_gp = bg_gp, internal_line = FALSE, side = side, ...)
	anno
}

# == title
# Word cloud annotations from GO
#
# == param
# -align_to The same format as in `anno_word_cloud`.
# -go_id The value should be in the same format as ``align_to``. If ``go_id`` is a vector, it should have the
#       same length as ``align_to``, and if ``go_id`` is a list, note, e.g. ``length(go_id[[1]])`` is not necessarily equal to ``length(align_to[[1]]``.
#       If ``align_to`` is a categorical vector and ``go_id`` is a list, names of ``go_id`` should have overlap to the levels in ``align_to``.
# -min_stat Minimal value for ``stat`` for selecting keywords.
# -stat What type of value to map to font sizes of the keywords. There are two possible values. "pvalue": enrichment is applied to keywords and -log10(p-value)
#       is used to map to font size; "count": simply word frequency of keywords.
# -term Alternatively the GO description can be set via the ``term`` argument. The same format as in `anno_word_cloud`.
# -exclude_words The words excluced for construcing word cloud. Some words are internally exclucded: ``c("via", "protein", "factor", "side", "type", "specific")``.
# -... All other arguments passed to `anno_word_cloud`.
#
anno_word_cloud_from_GO = function(align_to, go_id, stat = c("pvalue", "count"), 
	min_stat = ifelse(stat == "count", 5, 0.05),
	term = NULL, exclude_words = NULL, ...) {

	if(is.null(term)) {
		if(is.atomic(align_to) && is.list(go_id)) {
			align_to = split(seq_along(align_to), align_to)

			cn = intersect(names(align_to), names(go_id))
			if(length(cn) == 0) {
				stop_wrap("names of `go_id` should have overlap to the `align_to`.")
			} else {
				align_to = align_to[cn]
				go_id = go_id[cn]
			}

		} else if(is.list(align_to)) {
			if(!is.list(go_id)) {
				stop_wrap("`go_id` should have the same format as `align_to`, which is a list of term descriptions. Note, e.g. `length(go_id[[1]])` are not necessarily equal to `length(align_to[[1]]`.")
			} else {
				if(length(align_to) != length(go_id)) {
					stop_wrap("`align_to` and `go_id` should be two list with the same length.")
				}
			}
			if(!is.null(names(align_to)) && !is.null(names(go_id))) {
				if(length(setdiff(names(align_to), names(go_id))) == 0) {
					go_id = go_id[names(align_to)]
				}
			}
		} else if(is.atomic(align_to)) {
			if(!is.atomic(go_id)) {
				stop_wrap("`go_id` should have the same format as `align_to`, which is a vector of GO IDs with the same length as `align_to`.")
			}
			go_id = split(go_id, align_to)
			align_to = split(seq_along(align_to), align_to)
		}

		stat = match.arg(stat)[1]
		if(is.null(env$tdm_GO)) {
			# env$tdm_GO = readRDS("~/project/development/simplifyEnrichment/inst/extdata/tdm_GO.rds")
			env$tdm_GO = readRDS(system.file("extdata", "tdm_GO.rds", package = "simplifyEnrichment"))
		}

		if(stat == "pvalue") message(qq("Perform keywords enrichment for @{length(go_id)} GO lists..."))
		term = lapply(go_id, function(x) {
			if(is_GO_id(x[1])) {
				
				if(stat == "count") {
					v = row_sums(env$tdm_GO[, intersect(x, colnames(env$tdm_GO))])
					v = v[v >= min_stat]
					data.frame(names(v), v)
				} else if(stat == "pvalue") {
					df = keyword_enrichment(x, env$tdm_GO)
					df = df[df$p <= min_stat, , drop = FALSE]
					data.frame(df[, 1], -log10(df$p))
				}
			} else {
				stop_wrap("Cannot automatically retrieve the term names by the input ID, please set values for `term` argument manually.")
			}
			
		})
	} 

	anno_word_cloud(align_to, term, exclude_words = c(exclude_words, GO_EXCLUDE_WORDS), ...)
	
}
