
relabel_class = function (class, ref, full_set = union(class, ref), return_map = TRUE) {
    if(!requireNamespace("clue")) {
        stop_wrap("Package 'clue' should be installed.")
    }

    md = mode(class)
    class = as.character(class)
    ref = as.character(ref)
    full_set = as.character(full_set)
    all_class = union(class, ref)
    n = length(all_class)
    m = matrix(0, nrow = n, ncol = n, dimnames = list(all_class, all_class))
    tb = table(class, ref)
    m[rownames(tb), colnames(tb)] = tb
    imap = clue::solve_LSAP(m, maximum = TRUE)
    map = structure(rownames(m)[imap], names = rownames(m))
    map = map[unique(class)]
    unmapped = setdiff(setdiff(full_set, class), names(map))
    if (length(unmapped)) {
        map = c(map, structure(unmapped, names = unmapped))
    }
    df = data.frame(class = class, adjusted = map[class], ref = ref, 
        stringsAsFactors = FALSE)
    attr(map, "df") = df
    if (return_map) {
        return(map)
    }
    else {
        return(as(df$adjusted, md))
    }
}

stop_wrap = function (...) {
    x = paste0(...)
    x = paste(strwrap(x), collapse = "\n")
    stop(x, call. = FALSE)
}

message_wrap = function (...)  {
    x = paste0(...)
    x = paste(strwrap(x), collapse = "\n")
    message(x)
}
