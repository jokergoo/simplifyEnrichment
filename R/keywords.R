

GO_EXCLUDE_WORDS = c("via", "protein", "factor", "side", "type", "specific")


make_term_document_matrix = function(term,
	exclude_words = NULL, stop_words = stopwords(),
	min_word_length = 2, tokenizer = 'words', transform_case = tolower,
	remove_numbers = TRUE, remove_punctuation = TRUE, custom_transformer = NULL,
	stemming = FALSE, dictionary = NULL) {
	
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
				tolower = FALSE
			)
		)
	})

	tdm
}

prepare_keywords_tdm_for_GO = function(use_desc = FALSE, ...) {

	all_go = as.list(GO.db::GOTERM)

	ontology = sapply(all_go, slot, "Ontology")
	term = sapply(all_go, slot, "Term")
	l = ontology %in% c("BP", "CC", "MF")
	term = term[l]
	term_id = names(term)

	if(use_desc) {
		suppressMessages(term <- select(GO.db::GO.db, keys = term_id, columns = "DEFINITION")$DEFINITION)
		tdm = make_term_document_matrix(term, stop_words = c(stopwords(), GO_EXCLUDE_WORDS), ...)
	} else {
		tdm = make_term_document_matrix(term, stop_words = c(stopwords(), GO_EXCLUDE_WORDS), ...)
	}

	colnames(tdm) = term_id
	attr(tdm, "GO_dbInfo") = GO.db::GO_dbInfo()
	
	tdm
}

keywords_enrichment = function(term_id, tdm, min_bg = 5, min_set = 2) {
	tdm2 = tdm[row_sums(tdm) >= min_bg, ]

	l = colnames(tdm2) %in% term_id

	n = nrow(tdm2)
	n_term = numeric(n)
	n_bg = numeric(n)
	p = numeric(n)
	for(i in seq_len(n)) {
		if(interactive() && se_opt$verbose) {
			if(i %% 100 == 0 || i == n) {
				cat(strrep("\r", 100))
				qqcat("keyword enrichment, @{i}/@{n}...")
			}
		}
		v = as.vector(tdm2[i, ])
		s11 = sum(v & l)
		if(s11 < min_set) {
			next
		}
		s12 = sum(!v & l)
		s21 = sum(v & !l)
		s22 = sum(!v & !l)

		n_term[i] = s11
		n_bg[i] = s11 + s21

		p[i] = fisher.test(cbind(c(s11, s21), c(s12, s22)), alternative = "greater")$p.value

	}
	if(interactive() && se_opt$verbose) {
		cat("\n")
	}

	df = data.frame(keyword = rownames(tdm2), n_term = n_term, n_bg = n_bg, p = p)
	df = df[df$n_term >= min_set, , drop = FALSE]
	df$padj = p.adjust(df$p)
	df[order(df$padj, df$p), , drop = FALSE]
}

