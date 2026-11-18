decisions <- jsonlite::fromJSON("tools/param-advanced-decisions.json", simplifyVector = FALSE)$decisions
keys <- vapply(decisions, function(x) paste(x$file, x$class, x$parameter, sep = "::"), character(1))
files <- unique(vapply(decisions, `[[`, character(1), "file"))

# Locate parameter list entries in each config-definition call while respecting quotes/nesting.
parse_file <- function(path) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    starts <- grep("paramConfigDefsFact\\s*\\(\\s*list\\s*\\(", lines, perl = TRUE)
    out <- list()
    for (start in starts) {
        text <- paste(lines[start:length(lines)], collapse = "\n")
        call_start <- regexpr("paramConfigDefsFact\\s*\\(\\s*list\\s*\\(", text, perl = TRUE)[[1]]
        text <- substring(text, call_start)
        ch <- strsplit(text, "", fixed = TRUE)[[1]]
        depth <- 0L; quote <- ""; esc <- FALSE; call_end <- NA_integer_
        for (i in seq_along(ch)) {
            c <- ch[[i]]
            if (nzchar(quote)) {
                if (esc) esc <- FALSE else if (c == "\\") esc <- TRUE else if (c == quote) quote <- ""
                next
            }
            if (c %in% c("'", '"', "`")) quote <- c else if (c == "(") depth <- depth + 1L else if (c == ")") {
                depth <- depth - 1L
                if (depth == 0L) { call_end <- i; break }
            }
        }
        if (is.na(call_end)) next
        call <- substring(text, 1L, call_end)
        tail <- paste(lines[start:min(length(lines), start + 180L)], collapse = "\n")
        cm <- regexpr('setClass\\("[^"]+",\\s*contains\\s*=\\s*"param"', tail, perl = TRUE)
        cls <- if (cm[[1]] > 0) sub('^setClass\\("([^"]+)".*$', "\\1", regmatches(tail, cm)) else basename(path)
        locs <- gregexpr("[A-Za-z][A-Za-z0-9._]*\\s*=\\s*list\\s*\\(", call, perl = TRUE)[[1]]
        if (locs[[1]] < 0) next
        for (pos in locs) {
            prefix <- substring(call, 1L, pos - 1L)
            pch <- strsplit(prefix, "", fixed = TRUE)[[1]]
            d <- 0L; q <- ""; e <- FALSE
            for (c in pch) {
                if (nzchar(q)) { if (e) e <- FALSE else if (c == "\\") e <- TRUE else if (c == q) q <- ""; next }
                if (c %in% c("'", '"', "`")) q <- c else if (c == "(") d <- d + 1L else if (c == ")") d <- d - 1L
            }
            if (d != 2L) next
            name <- sub("^([A-Za-z][A-Za-z0-9._]*).*", "\\1", substring(call, pos))
            out[[length(out) + 1L]] <- data.frame(file = path, class = cls, parameter = name,
                line = start + sum(pch == "\n"), stringsAsFactors = FALSE)
        }
    }
    if (length(out)) do.call(rbind, out) else data.frame(file=character(), class=character(), parameter=character(), line=integer())
}

entries <- do.call(rbind, lapply(files, parse_file))
entry_keys <- paste(entries$file, entries$class, entries$parameter, sep = "::")
if (anyDuplicated(keys)) stop("Duplicate decision keys in JSON")
missing <- setdiff(keys, entry_keys)
if (length(missing)) stop("Decisions do not map uniquely to source entries: ", paste(missing, collapse = ", "))
if (anyDuplicated(entry_keys[match(keys, entry_keys)])) stop("Ambiguous source parameter definition")

for (path in files) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    es <- entries[entries$file == path, , drop = FALSE]
    # Process backwards, so inserted lines do not invalidate later source locations.
    es <- es[order(es$line, decreasing = TRUE), , drop = FALSE]
    for (j in seq_len(nrow(es))) {
        key <- paste(path, es$class[[j]], es$parameter[[j]], sep = "::")
        decision <- decisions[[match(key, keys)]]
        line <- es$line[[j]]
        # Entry extent ends before the next sibling parameter assignment.
        stop <- min(length(lines), line + 40L)
        following <- if (line < stop) seq.int(line + 1L, stop) else integer()
        siblings <- following[grepl("^[[:space:]]*[A-Za-z][A-Za-z0-9._]*[[:space:]]*=[[:space:]]*list[[:space:]]*\\(", lines[following], perl = TRUE)]
        if (length(siblings)) stop <- siblings[[1]] - 1L
        block <- seq.int(line, stop)
        flag <- paste0("advanced = ", if (isTRUE(decision$advanced)) "TRUE" else "FALSE", ",")
        existing <- block[grepl("^[[:space:]]*advanced[[:space:]]*=[[:space:]]*(TRUE|FALSE)[[:space:]]*,?", lines[block], perl = TRUE)]
        if (length(existing)) {
            lines[existing[[1]]] <- sub("advanced[[:space:]]*=[[:space:]]*(TRUE|FALSE)[[:space:]]*,?", flag, lines[existing[[1]]], perl = TRUE)
        } else {
            # Definitions may omit a default; add the flag at the beginning of the entry.
            di <- block[[1]] + 1L
            indent <- sub("^(\\s*).*", "\\1", lines[[di]])
            lines <- append(lines, paste0(indent, flag), after = di - 1L)
        }
    }
    writeLines(lines, path, useBytes = TRUE)
}

# Reparse definitions and confirm every decision's literal flag.
verified <- do.call(rbind, lapply(files, parse_file))
verified_keys <- paste(verified$file, verified$class, verified$parameter, sep = "::")
for (i in seq_along(decisions)) {
    decision <- decisions[[i]]
    key <- keys[[i]]
    row <- verified[verified_keys == key, , drop = FALSE]
    if (nrow(row) != 1L) stop("Validation did not find exactly one definition for ", key)
    source_lines <- readLines(decision$file, warn = FALSE, encoding = "UTF-8")
    start <- row$line[[1]]
    end <- min(length(source_lines), start + 40L)
    block <- source_lines[start:end]
    sibling <- which(grepl("^[[:space:]]*[A-Za-z][A-Za-z0-9._]*[[:space:]]*=[[:space:]]*list[[:space:]]*\\(", block, perl = TRUE))
    if (length(sibling) > 1L) block <- block[seq_len(sibling[[2]] - 1L)]
    flags <- grep("^[[:space:]]*advanced[[:space:]]*=[[:space:]]*(TRUE|FALSE)", block, value = TRUE)
    expected <- if (isTRUE(decision$advanced)) "TRUE" else "FALSE"
    if (length(flags) != 1L || !grepl(paste0("advanced[[:space:]]*=[[:space:]]*", expected), flags[[1]], perl = TRUE)) stop("Advanced-flag mismatch for ", key)
}
cat(sprintf("Updated and verified %d decisions in %d source files.\n", length(decisions), length(files)))
