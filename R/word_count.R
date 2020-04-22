
# == title
# Calculate word cloud
#
# == param
# -go_id A vector of GO IDs.
# -term The corresponding GO terms.
# -exclude_words The words that should be excluded.
#
# == value
# A data frame with words and frequencies.
#
count_word = function(go_id, term = NULL, exclude_words = NULL) {
	
	if(is.null(term)) term = select(GO.db::GO.db, keys = go_id, columns = "TERM")$TERM
	
	# http://www.sthda.com/english/wiki/word-cloud-generator-in-r-one-killer-function-to-do-everything-you-need

	# Load the text as a corpus
	docs = Corpus(VectorSource(term))
	# Convert the text to lower case
	docs = tm_map(docs, content_transformer(tolower))
	# Remove numbers
	docs = tm_map(docs, removeNumbers)
	# Remove stopwords for the language 
	docs = tm_map(docs, removeWords, stopwords())
	# Remove punctuations
	docs = tm_map(docs, removePunctuation)
	# Eliminate extra white spaces
	docs = tm_map(docs, stripWhitespace)
	# Remove your own stopwords
	docs = tm_map(docs, removeWords, c(exclude_words, EXCLUDE_WORDS))
	
	# Create term-document matrix
	tdm = TermDocumentMatrix(docs)

	v = sort(slam::row_sums(tdm), decreasing = TRUE)
	d = data.frame(word = names(v), freq = v, stringsAsFactors = FALSE)
	d
}


# generate excluded words that are too general
all_word_count = function() {
	all_go = as.list(GO.db::GOTERM)

	ontology = sapply(all_go, slot, "Ontology")
	term = sapply(all_go, slot, "Term")

	lt = tapply(term, ontology, function(x) {
		df = count_word(term = x)
		df = df[order(df$freq, decreasing = TRUE)[1:min(50, nrow(df))], ]
	})
	lt[c("BP", "CC", "MF")]
}

EXCLUDE_WORDS = c("via", "protein", "factor", "side", "type", "specific")

# == title
# A simple grob for word cloud
#
# == param
# -text A vector of words.
# -fontsize The corresponding font size.
# -max_width The maximal width of the viewport to put the word cloud. 
#            The value should be numeric. It is measured in mm.
#
# == value
# A `grid::grob` object. The width and height of the grob can be get by `grid::grobWidth` and `grid::grobHeight`.
#
# == example
# if(!exists("strrep")) {
# 	strrep = function(x, i) paste(rep(x, i), collapse = "")
# }
# words = sapply(1:30, function(x) strrep(sample(letters, 1), sample(3:10, 1)))
# require(grid)
# gb = simple_word_cloud_grob(words, fontsize = runif(30, min = 5, max = 30), max_width = 80)
# grid.newpage()
# w = grobWidth(gb)
# h = grobHeight(gb)
# pushViewport(viewport(width = w, height = grobHeight(gb)))
# grid.draw(gb)
# grid.rect()
# popViewport()
simple_word_cloud_grob = function(text, fontsize, max_width = 40) { # width in mm
	
	od = order(fontsize, decreasing = TRUE)
	text = text[od]
	fontsize = fontsize[od]

	n = length(text)
	text_gb_lt = lapply(seq_len(n), function(i) textGrob(text[i], gp = gpar(fontsize = fontsize[i])))
	text_width = sapply(text_gb_lt, function(gb) convertWidth(grobWidth(gb), "mm", valueOnly = TRUE))
	text_height = sapply(text_gb_lt, function(gb) convertHeight(grobHeight(gb), "mm", valueOnly = TRUE))

	margin = c(0.5, 0.5)

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

	for(i in seq_len(n)[-1]) {
		# the next text can be put on the same line
		if(current_line_width + text_width[i] <= max_width) {
			x[i] = current_line_width + margin[1]
			y[i] = y[i-1] # same as previous one
			current_line_width = x[i] + text_width[i]
			w = max(w, current_line_width)
			h = max(h, y[i] + text_height[i])
		} else { # the next text need to be put on the next line
			x[i] = 0
			y[i] = current_line_height + margin[2]
			current_line_width = text_width[i]
			current_line_height = y[i] + text_height[i]
			w = max(w, current_line_width)
			h = max(h, current_line_height)
		}
	}

	gl = gList(
		textGrob(text, x = x, y = y, gp = gpar(fontsize = fontsize, col = rand_color(n, luminosity = "dark")), 
			default.units = "mm", just = c(0, 0))
		# rectGrob(x = x, y = y, width = text_width, height = text_height, default.units = "mm", just = c(0, 0))
	)
	gb = gTree(children = gl, cl = "simple_word_cloud")
	attr(gb, "size") = c(w, h)
	return(gb)
}

# == title
# Width for simple_word_cloud Grob
#
# == param
# -x The ``simple_word_cloud`` grob returned by `simple_word_cloud_grob`.
#
widthDetails.simple_word_cloud = function(x) {
	unit(attr(x, "size")[1], "mm")
}

# == title
# Height for simple_word_cloud Grob
#
# == param
# -x The ``simple_word_cloud`` grob returned by `simple_word_cloud_grob`.
#
heightDetails.simple_word_cloud = function(x) {
	unit(attr(x, "size")[2], "mm")
}

	