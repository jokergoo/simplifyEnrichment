# Simplify GO enrichment Results

[![Build Status](https://travis-ci.org/jokergoo/simplifyGO.svg)](https://travis-ci.org/jokergoo/simplifyGO)


### Usage

As an example, we first generate a list of random GO IDs.

```r
set.seed(88)
go_id = random_GO(500)
head(go_id)
# [1] "GO:0042981" "GO:0000338" "GO:0006929" "GO:0043161" "GO:0006353"
# [6] "GO:0046101"
```

Then generate the GO similarity matrix and split GO terms into clusters.

```r
mat = GO_similarity(go_id)
simplify(mat)
```

![image](https://user-images.githubusercontent.com/449218/79051702-027a4d00-7c32-11ea-887e-ed3e171a03a0.png)

