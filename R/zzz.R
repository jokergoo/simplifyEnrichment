.onAttach = function(libname, pkgname) {
    version = packageDescription(pkgname, fields = "Version")

  	msg = paste0("========================================
", pkgname, " version ", version, "
Github page: https://github.com/jokergoo/simplifyEnrichment
Documentation: https://jokergoo.github.io/simplifyEnrichment/
Examples: https://jokergoo.github.io/simplifyGO_figures/

This message can be suppressed by:
  suppressPackageStartupMessages(library(simplifyEnrichment))
========================================
")	

    packageStartupMessage(msg)
}
