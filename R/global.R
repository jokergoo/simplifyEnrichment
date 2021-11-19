

# == title
# Global parameters
#
# == param
# -... Arguments for the parameters, see "details" section.
# -RESET Whether to reset to default values.
# -READ.ONLY Please ignore.
# -LOCAL Please ignore.
# -ADD Please ignore.
# 
# == details
# There are the following global options:
#
# -``verobse`` Whether to print messages.
#
se_opt = function(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE, ADD = FALSE) {}
se_opt = setGlobalOptions(
	verbose = TRUE
)
