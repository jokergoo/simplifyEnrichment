.onAttach = function(libname, pkgname) {
    version = packageDescription(pkgname, fields = "Version")

  	msg = paste0("========================================
", pkgname, " version ", version, "
Bioconductor page: https://bioconductor.org/packages/simplifyEnrichment/
Github page: https://github.com/jokergoo/simplifyEnrichment
Documentation: https://jokergoo.github.io/simplifyEnrichment/
Examples: https://simplifyenrichment.github.io/

This message can be suppressed by:
  suppressPackageStartupMessages(library(simplifyEnrichment))
========================================
")	

    packageStartupMessage(msg)

    ComplexHeatmap::ht_opt(save_last = TRUE)
}


GO_global_word_freq = readRDS(system.file("extdata", "GO_global_word_freq.rds", package = "simplifyEnrichment"))
