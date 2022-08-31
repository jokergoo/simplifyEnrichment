# # generate excluded words that are too general
# all_GO_word_count = function() {
# 	all_go = as.list(GO.db::GOTERM)

# 	ontology = sapply(all_go, slot, "Ontology")
# 	term = sapply(all_go, slot, "Term")
# 	l = ontology %in% c("BP", "CC", "MF")
# 	ontology = ontology[l]
# 	term = term[l]

# 	lt = tapply(term, ontology, function(x) {
# 		df = count_words(term = x, remove_punctuation = FALSE, remove_numbers = FALSE, transform_case = function(x) x)
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


# tokenizer = function(n_gram = 1) {
# 	if(n_gram == 1) {
# 		function(x) {
# 			w = x$content
# 			w = strsplit(w, "[[:punct:]]|[[:space:]]+", perl = TRUE)[[1]]
# 			w = w[!grepl("^\\d*$", w)]
# 			w
# 		}
# 	} else if(n_gram == 2) {
# 		function(x) {
# 			w = x$content
# 			w = strsplit(w, "[[:punct:]]|[[:space:]]+", perl = TRUE)[[1]]
# 			w = w[!grepl("^\\d*$", w)]
# 			n_grams(w, 2)
# 		}
# 	} else if(n_gram == 3) {
# 		function(x) {
# 			w = x$content
# 			w = strsplit(w, "[[:punct:]]|[[:space:]]+", perl = TRUE)[[1]]
# 			w = w[!grepl("^\\d*$", w)]
# 			n_grams(w, 3)
# 		}
# 	}
# }

# count_word2 = function(term, n_gram = 1) {

# 	n = length(term)

# 	docs = VCorpus(VectorSource(term))

# 	docs = tm_map(docs, removeWords, stopwords())

# 	# Create term-document matrix
# 	tdm = TermDocumentMatrix(
# 		docs,
# 		control = list(
# 			tokenize = tokenizer(n_gram),
# 			wordLengths = c(1, Inf),
# 			tolower = FALSE
# 		)
# 	)

# 	v = sort(slam::row_sums(tdm), decreasing = TRUE)

# 	data.frame(word = names(v), freq = v, stringsAsFactors = FALSE)
# }

# n_words = length(unlist(lapply(term, function(x) {
# 	x = list(content = x)
# 	setdiff(tokenizer(1)(x), stopwords())
# })))

# d1 = count_word2(term, 1)
# d2 = count_word2(term, 2)
# d3 = count_word2(term, 3)


# d2$p_expected = sapply(strsplit(d2$word, " "), function(w) {
# 	prod(d1[w, "p"])
# })


# d2 = d2[d2$freq >= 5, ]
# d2$ratio = d2$p/d2$p_expected




