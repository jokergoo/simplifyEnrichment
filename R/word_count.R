
count_word = function(go_id, term = NULL) {
	
	if(is.null(term)) term = select(GO.db::GO.db, keys = go_id, columns = "TERM")$TERM
	
	# http://www.sthda.com/english/wiki/word-cloud-generator-in-r-one-killer-function-to-do-everything-you-need

	# Load the text as a corpus
	docs <- Corpus(VectorSource(term))
	# Convert the text to lower case
	docs <- tm_map(docs, content_transformer(tolower))
	# Remove numbers
	docs <- tm_map(docs, removeNumbers)
	# Remove stopwords for the language 
	docs <- tm_map(docs, removeWords, stopwords())
	# Remove punctuations
	docs <- tm_map(docs, removePunctuation)
	# Eliminate extra white spaces
	docs <- tm_map(docs, stripWhitespace)
	# Remove your own stopwords
	# docs <- tm_map(docs, removeWords, .excludeWords) 
	
	# Create term-document matrix
	tdm <- TermDocumentMatrix(docs)

	v <- sort(slam::row_sums(tdm),decreasing = TRUE)
	d <- data.frame(word = names(v), freq = v, stringsAsFactors = FALSE)
	d
}

# k <- 30
# SEED <- 2010
# jss_TM <- list(
# VEM = LDA(tdm, k = k, control = list(seed = SEED)),
# VEM_fixed = LDA(tdm, k = k, control = list(estimate.alpha = FALSE,
#   seed = SEED)),
# Gibbs = LDA(tdm, k = k, method = "Gibbs", control = list(
#   seed = SEED, burnin = 1000, thin = 100, iter = 1000)),
# CTM = CTM(tdm, k = k, control = list(seed = SEED,
#   var = list(tol = 10^-4), em = list(tol = 10^-3))))


# # generate excluded words that are too general
# all_go = as.list(GO.db::GOTERM)

# ontology = sapply(all_go, slot, "Ontology")
# term = sapply(all_go, slot, "Term")

# lt = tapply(term, ontology, function(x) {
# 	df = count_word(term = x)
# 	df = df[df$freq > 1, ]
# })


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

	for(i in seq_len(n)[-1]) {
		# the next text can be put on the same line
		if(current_line_width + text_width[i] <= max_width) {
			x[i] = current_line_width + margin[1]
			y[i] = y[i-1] # same as previous one
			current_line_width = x[i] + text_width[i]
		} else { # the next text need to be put on the next line
			x[i] = 0
			y[i] = current_line_height + margin[2]
			current_line_width = text_width[i]
			current_line_height = y[i] + text_height[i]
		}
	}

	textGrob(text, x = x, y = y, gp = gpar(fontsize = fontsize), 
		default.units = "mm", just = c(0, 0))
}

	