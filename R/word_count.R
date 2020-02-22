library(tm)

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
	d <- data.frame(word = names(v), freq = v)
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

