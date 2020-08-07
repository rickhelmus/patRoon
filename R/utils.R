#' @include main.R
NULL

printf <- function(...) cat(sprintf(...), sep = "")
fprintf <- function(file, ..., append = FALSE) cat(sprintf(...), sep = "", file = file, append = append)

if (!exists("hasName")) # should be defined in latest R versions
    hasName <- function(x, name) return(name %in% names(x))

# backslashes are converted to forwardslashes for *nix compat
baseName <- function(...) basename(gsub("\\\\", "/", ...))

simplifyAnalysisNames <- function(slist)
{
    # simplify sample file names: remove extension and path
    return(as.character(sapply(slist, function(s) baseName(tools::file_path_sans_ext(s)), USE.NAMES = F)))
}

getMzMLAnalysisPath <- function(file, path) file.path(path, paste0(file, ".mzML"))
getMzXMLAnalysisPath <- function(file, path) file.path(path, paste0(file, ".mzXML"))
getMzDataAnalysisPath <- function(file, path) file.path(path, paste0(file, ".mzData"))
getAnalysisPath <- function(file, path, ext) file.path(path, paste0(file, ".", ext))

getMzMLOrMzXMLAnalysisPath <- function(file, path)
{
    ret <- getMzMLAnalysisPath(file, path)
    if (!file.exists(ret))
        ret <- getMzXMLAnalysisPath(file, path)
    return(ret)
}

mergeAnaSubsetArgWithSets <- function(i, sets, anaInfo)
{
    setAna <- anaInfo$analysis[anaInfo$set %in% sets]
    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, anaInfo$analysis)
        return(union(i, setAna))
    }
    return(setAna)
}

mkdirp <- function(path)
{
    for (p in path)
    {
        if (!dir.exists(p))
            dir.create(p, recursive = TRUE)
    }
}

insertDTColumn <- function(dt, col, d, before)
{
    dt <- copy(dt)

    dt[, (col) := d]

    if (before < ncol(dt))
    {
        ind <- c()

        if (before > 1)
            ind <- append(ind, 1:(before-1))

        ind <- append(ind, c(ncol(dt), before:(ncol(dt)-1)))
        setcolorder(dt, ind)
    }

    return(dt)
}

checkPackage <- function(pkg, gh = NULL)
{
    # from http://stackoverflow.com/a/20333756
    if (!requireNamespace(pkg, quietly = TRUE))
    {
        if (!is.null(gh))
            stop(sprintf("Please install %s from github: remotes::install_github('%s')", pkg, gh))
        else
            stop(sprintf("Please install %s: install.packages('%s')", pkg, pkg))
    }
}

executeCommand <- function(cmd, args = character(), ...)
{
    return(system2(cmd, sapply(args, shQuote), ...))
}

# NOTE: keep in sync with install-patRoon version
getCommandWithOptPath <- function(cmd, opt, verify = TRUE)
{
    if (Sys.info()[["sysname"]] == "Windows")
        cmd <- paste0(cmd, ".exe") # add file extension for Windows

    opt <- paste0("patRoon.path.", opt)
    path <- getOption(opt)
    if (!is.null(path) && nzchar(path))
    {
        ret <- file.path(path.expand(path), cmd)
        if (!file.exists(ret))
        {
            if (verify)
                stop(sprintf("Cannot find '%s'. Is the option '%s' set correctly?", ret, opt))
            return(NULL)
        }

        return(ret)
    }

    # assume command is in PATH --> no need to add path
    if (!nzchar(Sys.which(cmd)))
    {
        if (verify)
            stop(sprintf("Cannot find '%s'. Either add the correct file location to the PATH environment variable or set '%s' with options().", cmd, opt))
        return(NULL)
    }

    return(cmd)
}

# convert to unnamed character vector where previous names are followed by set values
OpenMSArgListToOpts <- function(args) as.vector(mapply(names(args), args, FUN = c, USE.NAMES = FALSE))

# NOTE: keep in sync with install-patRoon version
findPWizPath <- function()
{
    # try to find ProteoWizard
    # order: options --> win registry --> PATH
    # the PATH is searched last because OpenMS might have added its own old version.

    path <- getOption("patRoon.path.pwiz")
    if (!is.null(path) && nzchar(path))
        return(path)

    if (Sys.info()[["sysname"]] == "Windows")
    {
        # Inspired by scan_registry_for_rtools() from pkgload
        key <- "Directory\\shell\\Open with SeeMS\\command"
        reg <- tryCatch(utils::readRegistry(key, "HCR"), error = function(e) NULL)

        # not sure if this might occur
        if (is.null(reg))
            reg <- tryCatch(utils::readRegistry(key, "HLM"), error = function(e) NULL)

        if (!is.null(reg))
        {
            path <- tryCatch(dirname(sub("\"([^\"]*)\".*", "\\1", reg[[1]])), error = function(e) NULL)
            if (!is.null(path) && file.exists(file.path(path, "msconvert.exe"))) # extra check: see if msconvert is there
                return(path)
        }
    }

    # check PATH
    msc <- if (Sys.info()[["sysname"]] == "Windows") "msconvert.exe" else "msconvert"
    path <- dirname(Sys.which(msc))
    if (nzchar(path))
        return(path)

    return(NULL)
}

unFactorDF <- function(df)
{
    # from http://stackoverflow.com/a/2853231
    i <- sapply(df, is.factor)
    df[i] <- lapply(df[i], as.character)
    return(df)
}

# from http://stackoverflow.com/a/25455968
loadRData <- function(fileName)
{
    load(fileName)
    get(ls()[ls() != "fileName"])
}

# from https://stackoverflow.com/a/34788691
jgc <- function()
{
    gc()
    rJava::.jcall("java/lang/System", method = "gc")
}

# based on https://stackoverflow.com/a/17972280
recursiveApplyDT <- function(l, f, appl = lapply, ...)
{
    rec <- function(x)
    {
        if (isS4(x))
        {
            for (sn in slotNames(x))
                slot(x, sn) <- rec(slot(x, sn))
        }
        else if (is.list(x))
        {
            if (is.data.table(x))
                x <- f(x)
            else
            {
                # retain attributes: https://stackoverflow.com/a/48905113
                x <- "attributes<-"(appl(x, rec, ...), attributes(x))
            }
        }
        return(x)
    }

    return(rec(l))
}

prepareDTForComparison <- function(dt)
{
    setattr(dt, ".internal.selfref", NULL)
    setindex(dt, NULL)
}

readAllFile <- function(f) readChar(f, file.size(f))

getArgNames <- function(..., def = NULL)
{
    args <- sapply(substitute(list(...))[-1], deparse)
    n <- names(args)
    if (is.null(def))
        def <- args
    if (!is.null(n))
        args <- ifelse(nzchar(n), n, def)
    else if (!is.null(def))
        args <- def
    return(args)
}

getStrListWithMax <- function(l, m, collapse)
{
    if (length(l) > m)
    {
        l <- l[seq_len(m + 1)]
        l[m + 1] <- "..."
    }
    paste0(l, collapse = collapse)
}

showObjectSize <- function(object) printf("Object size (indication): %s\n", format(utils::object.size(object), "auto", "SI"))

allSame <- function(l)
{
    if (length(l) > 1)
    {
        if (all(is.na(l)))
            return(TRUE)
        if (any(is.na(l)))
            return(FALSE)

        return(all(sapply(l[-1], identical, l[[1]])))
    }

    return(TRUE)
}

# splitting a vector in chunks: https://stackoverflow.com/a/3321659
splitInBatches <- function(x, size) split(x, ceiling(seq_along(x) / size))
splitInNBatches <- function(x, n)
{
    if (n == 1)
        return(list(x))
    return(split(x, cut(seq_along(x), n, labels = FALSE)))
}

normalize <- function(x, minMax, xrange = range(x, na.rm = TRUE))
{
    xn <- x[!is.na(x)]
    if (length(xn) == 0 || all(xn == 0))
        return(x) # all NA or all values zero

    if (allSame(xn))
        xn <- rep(as(1, typeof(xn[[1]])), length(xn))
    else
    {
        minv <- xrange[1]
        if (!minMax)
            minv <- min(minv, 0) # force minMax if min <0

        xn <- (xn - minv) / (xrange[2] - minv)
    }

    x[!is.na(x)] <- xn

    return(x)
}

# from https://stackoverflow.com/a/38228840
countCharInStr <- function(str, ch) sum(charToRaw(str) == charToRaw(ch))

curTimeMS <- function() as.numeric(Sys.time()) * 1000

doDynamicTreeCut <- function(dendro, maxTreeHeight, deepSplit,
                             minModuleSize)
{
    if (minModuleSize == 1)
    {
        # workaround adapted from RAMClustR (ramclustR.R)
        ret <- dynamicTreeCut::cutreeDynamicTree(dendro = dendro, maxTreeHeight = maxTreeHeight,
                                                 deepSplit = deepSplit, minModuleSize = 2)
        single <- which(ret == 0) # all unassigned have length = 1
        ret[single] <- max(ret) + seq_len(length(single))
    }
    else
        ret <- dynamicTreeCut::cutreeDynamicTree(dendro = dendro, maxTreeHeight = maxTreeHeight,
                                                 deepSplit = deepSplit, minModuleSize = minModuleSize)

    return(ret)
}

# don't use all.equal here so functions are vectorized
# NOTE: actually all.equal seems to also take relative
# tolerances into account and thus may give different results.
numEQ <- function(x, y, tol = sqrt(.Machine$double.eps)) abs(x - y) <= tol
numGTE <- function(x, y, tol = sqrt(.Machine$double.eps)) numEQ(x, y, tol) | x > y
numLTE <- function(x, y, tol = sqrt(.Machine$double.eps)) numEQ(x, y, tol) | x < y

wrapStr <- function(s, width, sep = "\n") paste0(strwrap(s, width), collapse = sep)

pruneList <- function(l, checkEmptyElements = FALSE, checkZeroRows = FALSE)
{
    ret <- l[!sapply(l, is.null)]
    if (checkEmptyElements)
        ret <- ret[lengths(ret) > 0]
    if (checkZeroRows)
        ret <- ret[sapply(ret, nrow) > 0]
    return(ret)
}

makeEmptyListNamed <- function(li)
{
    if (length(li) == 0)
        names(li) <- character()
    return(li)
}

# based on tabular() from formatting vignette of roxygen
# nocov start
tabularRD <- function(df, ...)
{
    align <- function(x) if (is.numeric(x)) "r" else "l"
    col_align <- vapply(df, align, character(1))

    # add headers
    df <- rbind(sprintf("\\strong{%s}", names(df)), df)

    cols <- lapply(df, format, ...)

    contents <- do.call("paste",
                        c(cols, list(sep = " \\tab ", collapse = "\\cr\n  ")))

    paste("\\tabular{", paste(col_align, collapse = ""), "}{\n  ",
          contents, "\n}\n", sep = "")
}

# make a S4 class inheritance tree in a format compatible with data.tree::FromListSimple()
makeClassHierarchy <- function(class, showParents)
{
    cldef <- getClassDef(class)
    subcl <- list()
    if (length(cldef@subclasses) > 0)
        subcl <- cldef@subclasses[sapply(cldef@subclasses, function(sc) sc@distance == 1)]

    ret <- c(list(name = class),
             lapply(names(subcl), makeClassHierarchy, showParents = FALSE))

    if (showParents)
    {
        pars <- selectSuperClasses(class)
        if (length(pars) > 0)
            ret <- c(list(name = paste0(pars, collapse = ", ")), list(ret))
    }

    return(ret)
}

printClassHierarchy <- function(class, showParents = TRUE, RD = FALSE)
{
    doPrintTxt <- function(cl, level)
    {
        indent <- strrep(" ", (level + 1) * 2)

        if (level > 0)
            cat(paste0(indent, "|-- ", cl$name))
        else
            cat(cl$name)
        cat("\n")

        for (clsub in cl)
        {
            if (is.list(clsub))
                doPrintTxt(clsub, level + 1)
        }
    }

    doPrintRD <- function(cl, level)
    {
        indent <- strrep(" ", (level + 1) * 2)

        printf("%s%s\\item{%s}\n", if (level == 0) "\\itemize{\n" else "", indent, cl$name)

        more <- any(sapply(cl, is.list))
        if (more)
            cat(paste0(indent, "\\itemize{\n"))

        for (clsub in cl)
        {
            if (is.list(clsub))
                doPrintRD(clsub, level + 1)
        }

        if (more)
            cat(paste0(indent, "}\n"))
        if (level == 0)
            cat("}\n")
    }

    hier <- makeClassHierarchy(class, showParents)
    hasParents <- hier$name != class

    if (RD)
    {
        hier <- rapply(hier, function(h) sprintf("\\code{\\link{%s}}", h),
                       classes = "character", how = "replace")

        # if parents are shown make the second line (ie this class) bold
        if (hasParents)
            hier[[2]]$name <- sprintf("\\strong{%s}", hier[[2]]$name)
        else
            hier$name <- sprintf("\\strong{%s}", hier$name)

        doPrintRD(hier, 0)
    }
    else
        doPrintTxt(hier, 0)

    invisible(NULL)
}

getAllMethods <- function(gen)
{
    # automatically retrieve defined methods for a generic and create document
    # links. This only works if the arguments of the method are named obj, objX or x.

    cl <- showMethods(gen, where = "package:patRoon", inherited = FALSE, printTo = FALSE,
                      classes = getClasses(asNamespace("patRoon")))
    cl <- cl[grepl("obj.*|x=", cl)]
    cl <- gsub("[^\"]*\"([^\"]*)\"[^\"]*", "\\1,", cl)
    # cl <- cl[!grepl("ANY", cl)]
    cl <- gsub(" ", "", cl)
    cl <- gsub(",$", "", cl)

    return(cl[order(tolower(cl))])
}
# nocov end

NULLToZero <- function(x) if (is.null(x)) 0 else x
zeroToNULL <- function(x) if (is.numeric(x) && x == 0) NULL else x
NAToZero <- function(x) if (is.na(x)) 0 else x

# From https://stackoverflow.com/a/47955845
allArgs <- function(origValues = FALSE)
{
    # get formals for parent function
    parent_formals <- formals(sys.function(sys.parent(n = 1)))

    # Get names of implied arguments
    fnames <- names(parent_formals)

    # Remove '...' from list of parameter names if it exists
    fnames <- fnames[-which(fnames == '...')]

    # Get currently set values for named variables in the parent frame
    args <- evalq(as.list(environment()), envir = parent.frame())

    # Get the list of variables defined in '...'
    args <- c(args[fnames], evalq(list(...), envir = parent.frame()))

    if(origValues)
    {
        # get default values
        defargs <- as.list(parent_formals)
        defargs <- defargs[unlist(lapply(defargs, FUN = function(x) class(x) != "name"))]
        args[names(defargs)] <- defargs
        setargs <- evalq(as.list(match.call())[-1], envir = parent.frame())
        args[names(setargs)] <- setargs
    }

    return(args)
}

openProgBar <- function(min = 0, max, style = 3, ...)
{
    progOpts <- list(min = min, max = max, style = style, ...)
    progOpts <- modifyList(progOpts, getOption("patRoon.progress.opts", list()))
    return(do.call(txtProgressBar, progOpts))
}

verboseCall <- function(f, a, v) if (v) do.call(f, a) else suppressMessages(invisible(do.call(f, a)))

RUserDir <- function(...)
{
    if (getRversion() >= "4.0.0")
        return(tools::R_user_dir(...))
    return(backports:::R_user_dir(...))
}

assertAndGetMSPLSetsArgs <- function(fGroupsSet, MSPeakListsSet)
{
    checkmate::assertClass(MSPeakListsSet, "MSPeakListsSet")
    sd <- setdiff(sets(fGroupsSet), sets(MSPeakListsSet))
    if (length(sd) > 0)
        stop(paste("The following sets in fGroups are missing in MSPeakLists:"), paste0(sd, collapse = ", "))
    
    ionizedMSPeaksList <- lapply(sets(MSPeakListsSet), ionize, obj = MSPeakListsSet)
    return(lapply(ionizedMSPeaksList, function(x) list(MSPeakLists = x)))
}
