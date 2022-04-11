# # generate excluded words that are too general
# all_GO_word_count = function() {
# 	all_go = as.list(GO.db::GOTERM)

# 	ontology = sapply(all_go, slot, "Ontology")
# 	term = sapply(all_go, slot, "Term")
# 	l = ontology %in% c("BP", "CC", "MF")
# 	ontology = ontology[l]
# 	term = term[l]

# 	lt = tapply(term, ontology, function(x) {
# 		df = count_word(term = x)
# 		structure(df[, 2], names = df[, 1])
# 	})
# }


# n_grams = function(sequence, n) {
#     if (n <= 0) {
#         stop('n must be greater than 0')
#     }
#     if (length(sequence) < n) {
#         return(character())
#     }
#     indices = seq(0, n - 1)
#     sapply(
#         seq(1, length(sequence) - n + 1),
#         function(start) {
#             paste(sequence[indices + start], collapse=' ')
#         }
#     )
# }


# count_word2 = function(term, n_gram = 1) {

# 	n = length(term)

# 	#### analyze text ###########
# 	docs = VCorpus(VectorSource(term))
# 	# Convert the text to lower case
# 	docs = tm_map(docs, content_transformer(tolower))
# 	# Remove stopwords for the language
# 	# docs = tm_map(docs, removeWords, stopwords())
# 	# docs = tm_map(docs, removePunctuation)
# 	# docs = tm_map(docs, stripWhitespace)

# 	if(n_gram == 1) {
# 		tokenizer = scan_tokenizer
# 		docs = tm_map(docs, removeWords, stopwords::stopwords("en", source = "stopwords-iso"))
# 	} else if(n_gram == 2) {
# 		tokenizer = function(x) n_grams(scan_tokenizer(x), 2)
# 	} else if(n_gram == 3) {
# 		tokenizer = function(x) n_grams(scan_tokenizer(x), 3)
# 	}

# 	# Create term-document matrix
# 	tdm = TermDocumentMatrix(
# 		docs,
# 		control = list(
# 			tokenize = tokenizer,
# 			wordLengths = c(2, Inf),
# 			tolower = FALSE
# 		)
# 	)

# 	v = sort(slam::row_sums(tdm), decreasing = TRUE)

# 	if(n_gram >= 2) {
# 		v1 = count_word2(term, 1)
# 		p2 = sapply(strsplit(names(v), " "), function(x) prod((v1[x]+5)/n))
# 		v = cbind(v = v, p1 = (v+5)/n, p2 = p2)
# 	}
# 	v
# }

