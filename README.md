# Simplify Functional Enrichment Results

[![Build Status](https://travis-ci.org/jokergoo/simplifyEnrichment.svg)](https://travis-ci.org/jokergoo/simplifyEnrichment)
[![bioc](http://www.bioconductor.org/shields/downloads/devel/simplifyEnrichment.svg)](https://bioconductor.org/packages/stats/bioc/simplifyEnrichment/) 
[![bioc](http://www.bioconductor.org/shields/years-in-bioc/simplifyEnrichment.svg)](http://bioconductor.org/packages/devel/bioc/html/simplifyEnrichment.html)

### Features

- A new method (binary cut) is proposed to efficiently cluster functional terms (_e.g._ GO terms) into groups from the semantic similarity matrix.
- Summaries of functional terms in each cluster are visualized by word clouds.


### Install

`simplifyEnrichment` is available on [Bioconductor](http://www.bioconductor.org/packages/devel/bioc/html/simplifyEnrichment.html), you can install it by:

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("simplifyEnrichment")
```

If you want to try the latest version, install it directly from GitHub:

```r
library(devtools)
install_github("jokergoo/simplifyEnrichment")
```

### Vignette

- [Simplify Functional Enrichment Results](https://jokergoo.github.io/simplifyEnrichment/articles/simplifyEnrichment.html)
- [Word Cloud Annotation](https://jokergoo.github.io/simplifyEnrichment/articles/word_cloud_anno.html)

### Usage

As an example, I first generate a list of random GO IDs.

```r
library(simplifyEnrichment)
set.seed(888)
go_id = random_GO(500)
head(go_id)
# [1] "GO:0003283" "GO:0060032" "GO:0031334" "GO:0097476" "GO:1901222"
# [6] "GO:0018216"
```

Then generate the GO similarity matrix, split GO terms into clusters and visualize it.

```r
mat = GO_similarity(go_id)
simplifyGO(mat)
```

![image](https://user-images.githubusercontent.com/449218/89673686-133c8600-d8e7-11ea-89fe-5221cb64d819.png)

### Examples

- [Examples of simplifyEnrichment](https://simplifyenrichment.github.io/examples/)
- [Compare different similarity measures for functional terms](https://simplifyenrichment.github.io/compare_similarity/)
- [Compare different partitioning methods in binary cut clustering](https://simplifyenrichment.github.io/test_partition_methods/)

### License

MIT @ Zuguang Gu
