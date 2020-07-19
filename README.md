# Simplify functional enrichment Results

[![Build Status](https://travis-ci.org/jokergoo/simplifyEnrichment.svg)](https://travis-ci.org/jokergoo/simplifyEnrichment)
[![bioc](http://www.bioconductor.org/shields/downloads/devel/simplifyEnrichment.svg)](https://bioconductor.org/packages/stats/bioc/simplifyEnrichment/) 
[![bioc](http://www.bioconductor.org/shields/years-in-bioc/simplifyEnrichment.svg)](http://bioconductor.org/packages/devel/bioc/html/simplifyEnrichment.html)

### Features

- A new method (binary cut) is proposed to effectively cluster functional terms (e.g. GO terms) into groups from the semantic similarity matrix.
- Summaries of GO terms in each cluster are visualized by word clouds.


### Install

`simplifyEnrichment` is available on [Bioconductor](http://www.bioconductor.org/packages/devel/bioc/html/simplifyEnrichment.html), you can install it by:

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("simplifyEnrichment")
```

If you want the latest version, install it directly from GitHub:

```r
library(devtools)
install_github("jokergoo/simplifyEnrichment")
```

### Vignette

- [Simplify Gene Ontology Enrichment Results](https://jokergoo.github.io/simplifyEnrichment/articles/simplifyGO.html)

### Usage

As an example, we first generate a list of random GO IDs.

```r
library(simplifyEnrichment)
set.seed(88)
go_id = random_GO(500)
head(go_id)
# [1] "GO:0042981" "GO:0000338" "GO:0006929" "GO:0043161" "GO:0006353" "GO:0046101"
```

Then generate the GO similarity matrix, split GO terms into clusters and visualize it.

```r
mat = GO_similarity(go_id)
simplifyGO(mat)
```

![image](https://user-images.githubusercontent.com/449218/81970861-78eae000-9620-11ea-9ea2-c9d51d36ff4c.png)

### Examples

- [Randomly generated GO terms](https://jokergoo.github.io/simplifyGO_figures/random_BP.html)
- [Enriched GO terms from EBI Expression Atlas datasets](https://jokergoo.github.io/simplifyGO_figures/EBI_Expression_Atlas.html)

### License

MIT @ Zuguang Gu
