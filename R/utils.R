#' @include main.R
NULL

printf <- function(...) cat(sprintf(...), sep = "")
fprintf <- function(file, ..., append = FALSE) cat(sprintf(...), sep = "", file = file, append = append)
mdprintf <- function(...) cat(sprintf(...), sep = "", file = stderr()) # for rmarkdown

if (!exists("hasName")) # should be defined in latest R versions
    hasName <- function(x, name) return(name %in% names(x))

makeDT <- function(df) if (is.data.table(df)) copy(df) else as.data.table(df)

# backslashes are converted to forwardslashes for *nix compat
baseName <- function(...) basename(gsub("\\\\", "/", ...))

simplifyAnalysisNames <- function(slist)
{
    # simplify sample file names: remove extension and path
    return(as.character(sapply(slist, function(s) baseName(tools::file_path_sans_ext(s)), USE.NAMES = F)))
}

getMzMLAnalysisPath <- function(file, path, mustExist = FALSE) getMSFilePaths(file, path, "mzML", mustExist)
getMzXMLAnalysisPath <- function(file, path, mustExist = FALSE) getMSFilePaths(file, path, "mzXML", mustExist)

getMzMLOrMzXMLAnalysisPath <- function(file, path, mustExist = FALSE)
{
    ret <- getMzMLAnalysisPath(file, path)
    if (length(ret) == 0)
        ret <- getMzXMLAnalysisPath(file, path)
    if (mustExist && length(ret) == 0)
        stop(sprintf("Analysis %s is not found as mzML or mzXML file.", file), call. = FALSE)
    return(ret)
}

mergeAnaSubsetArgWithSets <- function(i, sets, anaInfo)
{
    if (length(sets) == 0)
        return(character())
    
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

moveLastDTColumn <- function(dt, before)
{
    if (before < ncol(dt))
    {
        dt <- copy(dt)
        
        ind <- if (before > 1) seq_len(before - 1) else NULL
        ind <- append(ind, c(ncol(dt), seq(before, ncol(dt) - 1)))
        setcolorder(dt, ind)
    }

    return(dt)
}

checkPackage <- function(pkg, gh = NULL, ghSubDir = NULL)
{
    # from http://stackoverflow.com/a/20333756
    if (!requireNamespace(pkg, quietly = TRUE))
    {
        if (!is.null(gh))
        {
            args <- if (!is.null(ghSubDir)) sprintf("'%s', subdir = '%s')", gh, ghSubDir) else sprintf("'%s'", gh)
            stop(sprintf("Please install %s from github: remotes::install_github(%s)", pkg, args), call. = FALSE)
        }
        else
            stop(sprintf("Please install %s: install.packages('%s')", pkg, pkg), call. = FALSE)
    }
}

executeCommand <- function(cmd, args = character(), ...)
{
    return(system2(cmd, sapply(args, shQuote), ...))
}

getExtDepPath <- function(what, subTool = NULL, verify = TRUE)
{
    # NOTE: pwiz is handled by findPWizPath()
    # NOTE: GenForm is handled by getGenFormBin()
    
    exts <- list(
        openms = list(
            name = "OpenMS",
            bin = subTool,
            opt = "OpenMS"
        ),
        sirius = list(
            name = "SIRIUS",
            bin = "sirius",
            opt = "SIRIUS"
        ),
        pngquant = list(
            name = "pngquant",
            bin = "pngquant",
            opt = "pngquant"
        ),
        openbabel = list(
            name = "OpenBabel",
            bin = "obabel",
            opt = "obabel"
        ),
        metfragcl = list(
            name = "MetFrag",
            opt = "MetFragCL"
        ),
        metfragct = list(
            name = "MetFrag CompTox DB",
            opt = "MetFragCompTox"
        ),
        metfragpcl = list(
            name = "MetFrag PubChemLite DB",
            opt = "MetFragPubChemLite"
        ),
        biotransformer = list(
            name = "BioTransformer",
            opt = "BioTransformer"
        )
    )
    
    ext <- exts[[what]]
    
    hasBin <- !is.null(ext[["bin"]])
    
    if (hasBin && Sys.info()[["sysname"]] == "Windows")
        ext$bin <- paste0(ext$bin, ".exe") # add file extension for Windows
    
    # order: options(), patRoonExt, PATH

    ext$opt <- paste0("patRoon.path.", ext$opt)
    path <- getOption(ext$opt)
    if (!is.null(path) && nzchar(path))
    {
        path <- path.expand(path)
        if (hasBin)
            path <- file.path(path, ext$bin)
        if (!file.exists(path))
        {
            if (verify)
                stop(sprintf("Cannot find '%s'. Is the option '%s' set correctly?", path, ext$opt))
            return(NULL)
        }
        
        return(path)
    }
    
    if (requireNamespace("patRoonExt", quietly = TRUE))
    {
        path <- patRoonExt::getExtPath(what, warn = FALSE)
        if (!is.null(path))
        {
            if (hasBin)
                path <- file.path(path, ext$bin)
            return(path)
        }
    }
    
    # only check binaries for PATH
    if (!hasBin || !nzchar(Sys.which(ext$bin)))
    {
        if (verify)
            stop(sprintf("Cannot find '%s'. You may need to install patRoonExt, set '%s' with options() or add the correct file location to the PATH environment variable.",
                         ext$name, ext$opt), call. = FALSE)
        return(NULL)
        
    }

    return(ext$bin) # if we're here it's a binary found in PATH
}

# convert to unnamed character vector where previous names are followed by set values
OpenMSArgListToOpts <- function(args) as.vector(mapply(names(args), args, FUN = c, USE.NAMES = FALSE))

OpenMSVersionAtLeast <- function(tool, minVersion)
{
    toolHelp <- suppressWarnings(executeCommand(getExtDepPath("openms", tool), stdout = TRUE, stderr = TRUE))
    ver <- grep("^Version\\:", toolHelp, value = TRUE)
    if (length(ver) == 1)
    {
        ver <- unlist(regmatches(ver, regexec("[0-9\\.]+", ver)))
        return(utils::compareVersion(ver, minVersion) >= 0)
    }
    NA
}

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

# as above
unFactorDT <- function(dt)
{
    cols <- names(which(sapply(dt, is.factor)))
    if (length(cols) > 0)
        dt[, (cols) := lapply(.SD, as.character), .SDcols = cols]
    return(dt[])
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

removeDTColumnsIfPresent <- function(dt, cols) dt[, setdiff(names(dt), cols), with = FALSE]
subsetDTColumnsIfPresent <- function(dt, cols) dt[, intersect(names(dt), cols), with = FALSE]

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

allSame <- function(l, func = identical)
{
    if (length(l) > 1)
    {
        if (all(is.na(l)))
            return(TRUE)
        if (any(is.na(l)))
            return(FALSE)

        return(all(sapply(l[-1], func, l[[1]])))
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
            ret <- c(list(name = paste0(pars, collapse = ",")), list(ret))
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

        # NOTE: don't go over four levels (0-3) since TeX cannot handle more 
        if (level < 3 && any(sapply(cl, is.list)))
        {
            cat(paste0(indent, "\\itemize{\n"))
            for (clsub in cl)
            {
                if (is.list(clsub))
                    doPrintRD(clsub, level + 1)
            }
            cat(paste0(indent, "}\n"))
        }
        if (level == 0)
            cat("}\n")
    }

    hier <- makeClassHierarchy(class, showParents)
    hasParents <- hier$name != class

    if (RD)
    {
        hier <- rapply(hier, function(h) paste0(sprintf("\\code{\\link{%s}}", unlist(strsplit(h, ","))), collapse = ","),
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

callAllNextMethods <- function(obj, f, ..., firstClass = NULL, startFrom = class(obj))
{
    supers <- selectSuperClasses(startFrom, directOnly = TRUE)
    if (!is.null(firstClass))
    {
        selectMethod(f, firstClass)(obj, ...)
        supers <- supers[supers != firstClass]
    }
    for (sc in supers)
        selectMethod(f, sc)(obj, ...)
}

setMethodMult <- function(f, signatures, definition)
{
    for (sig in signatures)
        setMethod(f, sig, definition)
}

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
    return(getFromNamespace("R_user_dir", "backports")(...))
}

ReduceWithArgs <- function(f, x, ..., fixedArgs = list())
{
    ret <- x[[1]]
    itArgs <- list(...)
    
    if (length(x) > 1)
    {
        x <- x[-1]
        for (i in seq_along(x))
            ret <- do.call(f, c(list(ret, x[[i]]),
                                lapply(itArgs, "[[", 1), lapply(itArgs, "[[", i+1),
                                fixedArgs))
    }
    
    return(ret)
}

# Based on https://www.inwt-statistics.com/read-blog/optimize-your-r-code-using-memoization.html
memoise <- function(fun)
{
    memory <- list()
    function(...)
    {
        valueName <- if (...length() == 1) as.character(..1) else makeHash(..., checkDT = FALSE)
        if (!is.null(memory[[valueName]]))
            return(memory[[valueName]])
        res <- fun(...)
        memory[[valueName]] <<- res
        res
    }
}

getHighestAbsValue <- function(abs, rel, size)
{
    abs <- NULLToZero(abs); rel <- NULLToZero(rel)
    return(max(abs, rel * size))
}

readYAML <- function(f) yaml::read_yaml(f, eval.expr = FALSE)
writeYAML <- function(...) yaml::write_yaml(..., indent = 4)

# get a vector of all (merged) columns
getAllMergedConsCols <- function(targetCols, allCols, mConsNames)
{
    if (length(mConsNames) > 0)
        targetCols <- c(targetCols, paste0(targetCols, "-", mConsNames))
    return(intersect(targetCols, allCols))
}

getDuplicatedStrings <- function(x) names(which(table(x) > 1))

doApply <- function(applyf, doPar, data, ...)
{
    if (doPar)
        applyf <- get(paste0("future_", applyf), envir = asNamespace("future.apply"))
    withProg(length(data), doPar, do.call(applyf, list(data, ...)))
}

getMS2QuantRes <- function(calibrants, unknowns, eluent, organicModifier, pHAq, allFPs)
{
    calibrants <- copy(calibrants)
    setnames(calibrants, c("name", "rt", "intensity", "conc"), c("identifier", "retention_time", "area", "conc_M"))
    
    # UNDONE: would be nice if we could just pass table directly
    quantFile <- tempfile(fileext = ".csv"); fwrite(rbind(calibrants, unknowns, fill = TRUE), quantFile)
    eluentFile <- tempfile(fileext = ".csv"); fwrite(eluent, eluentFile)
    suppressMessages(utils::capture.output(pr <- MS2Quant::MS2Quant_quantify(quantFile, eluentFile,
                                                                             organic_modifier = organicModifier,
                                                                             pH_aq = pHAq, fingerprints = allFPs)))
    
    RFs <- as.data.table(pr$suspects_concentrations)
    
    # the area we set to one, we calculated the concentration for that response, so the inverse is the RF (=area/conc=1/conc)
    RFs[, RF_pred := 1 / conc_M]
    
    RFs <- RFs[, c("identifier", "SMILES", "RF_pred"), with = FALSE]
    
    return(list(RFs = RFs[], MD = list(calibPlot = pr$logIE_logRF_calibration_plot,
                                       linModel = pr$calibration_linear_model_summary)))
}

aggrVec <- function(x, f)
{
    x <- x[!is.na(x)]
    if (length(x) == 0)
        is.na(x) <- TRUE # use is.na()<- to retain class of x
    else
        x <- f(x)
    return(x)
}

readIDLRules <- function(IDFile)
{
    ret <- readYAML(IDFile)
    
    if (!checkmate::test_named(ret))
        stop("No valid rules could be loaded")
    if (!all(grepl("^[[:digit:]]+[[:alpha:]]?$", names(ret))))
        stop("Levels should be defined as a number and may optionally followed by one character (e.g. 3, 2b etc)")
    
    ret <- ret[order(names(ret))] # sort to ensure lowest levels will be tested first
    
    return(ret)
}

estimateIdentificationLevel <- function(candidateName, candidateFGroup, candidateRTDev, candidateInChIKey1, candidateFormula,
                                        candidateAnnSimForm, candidateAnnSimComp, candidateAnnSimBoth,
                                        maxSuspFrags, maxFragMatches, formTable, formTableNorm, formRank, mFormNames,
                                        compTable, compTableNorm, compRank, mCompNames, absMzDev, IDLevelRules, logPath)
{
    formScores <- formScoreNames(FALSE)
    compScores <- compScoreNames(FALSE)
    fRow <- cRow <- NULL
    
    if (!is.null(formTable) && !is.null(candidateFormula) && !is.na(formRank))
    {
        fRow <- formTable[formRank]
        fRowNorm <- formTableNorm[formRank]
    }
    
    if (!is.null(compTable) && !is.null(candidateInChIKey1) && !is.na(compRank))
    {
        cRow <- compTable[compRank]
        cRowNorm <- compTableNorm[compRank]
    }
    
    getValType <- function(val, IDType)
    {
        if (!is.list(val) || is.null(val[["type"]]) || !val[["type"]] %in% c("formula", "compound"))
            stop(sprintf("Need to specify the type (formula/compound) for %s!", IDType))
        return(val[["type"]])
    }
    
    getVal <- function(val, IDType, valType)
    {
        if (is.list(val))
        {
            if (is.null(val[[valType]]))
                stop(sprintf("Need to specify %s for %s!", valType, IDType))
            return(val[[valType]])
        }
        return(val)
    }
    
    getOptVal <- function(val, valType, default)
    {
        if (is.list(val) && !is.null(val[[valType]]))
            return(val[[valType]])
        return(default)
    }
    
    checkAnnotationScore <- function(val, scType, rank, annRow, annTable, annRowNorm, annTableNorm, scCols)
    {
        scCols <- scCols[!is.na(unlist(annRow[, scCols, with = FALSE]))]
        if (length(scCols) == 0)
            return("score not available")
        
        minValue <- getVal(val, scType, "min")
        
        if (getOptVal(val, "relative", FALSE))
        {
            annRow <- annRowNorm
            annTable <- annTableNorm
        }
        
        scoreVal <- rowMeans(annRow[, scCols, with = FALSE], na.rm = TRUE)
        if (is.infinite(scoreVal)) # only NA values
            return("no value(s) available for this score")
        if (scoreVal < minValue)
            return(sprintf("(average) score too low: %f/%f", scoreVal, minValue))
        
        htn <- getOptVal(val, "higherThanNext", 0)
        if (htn > 0 && nrow(annTable) > 1)
        {
            otherVals <- rowMeans(annTable[-rank, scCols, with = FALSE], na.rm = TRUE)
            otherVals <- otherVals[is.finite(otherVals)] # remove results for NA rows
            if (length(otherVals) > 0)
            {
                otherHighest <- max(otherVals)
                if (is.infinite(htn)) # special case: should be highest
                {
                    if (otherHighest > 0)
                        return("not the highest score")
                }
                else if ((scoreVal - otherHighest) < htn)
                    return(sprintf("difference with highest score from other candidates too low: %f/%f", scoreVal - otherHighest, htn))
            }
        }
        
        return(TRUE)            
    }
    
    if (!is.null(logPath))
    {
        logDir <- R.utils::getAbsolutePath(logPath) # take absolute path for length calculation below
        logFile <- paste0(candidateName, "-", candidateFGroup, ".txt")
        
        # check if path would be too long for e.g Windows systems, which may happen with very long candidate names
        pathLen <- nchar(logDir) + nchar(logFile) + 1 # +1: path separator
        if (pathLen > 255)
        {
            # truncate end part of candidate name
            logFile <- paste0(substr(candidateName, 1, nchar(candidateName) - (pathLen - 255)), "-", candidateFGroup, ".txt")
        }
        logOut <- file.path(logDir, logFile)
        
        mkdirp(dirname(logOut))
        logFile <- withr::local_connection(file(logOut, "w"))
        doLog <- function(indent, s, ...) fprintf(logFile, paste0(strrep(" ", indent * 4), s), ...)
    }
    else
        doLog <- function(...) NULL
    
    formCompScores <- intersect(formScores, compScores)
    allScores <- union(formScores, compScores)
    
    checkLevelOK <- function(IDL, indent = 0)
    {
        indent <- indent + 1
        for (type in names(IDL))
        {
            levelOK <- NULL
            levelFailed <- NULL
            val <- IDL[[type]]
            
            if (type %in% c("rank", "annMSMSSim", formCompScores))
                doLog(indent, "Checking ID level type '%s' (for %s)\n", type,
                      getValType(val, type))
            else
                doLog(indent, "Checking ID level type '%s'\n", type)
            
            if (type %in% c("or", "and"))
            {
                if (!is.list(val) || checkmate::testNamed(val))
                    stop(sprintf("Specify a list with '%s'", type))
                checkf <- if (type == "or") any else all
                levelOK <- checkf(mapply(val, seq_along(val), FUN = function(IDL, i)
                {
                    doLog(indent + 1, "check %s condition %d/%d\n", toupper(type), i, length(val))
                    return(checkLevelOK(IDL, indent + 2))
                }))
            }
            else if (type == "all" && val == TRUE)
                levelOK <- TRUE # special case: this level is always valid
            else if (type == "suspectFragments")
            {
                minFrags <- min(val, maxSuspFrags, na.rm = TRUE)
                levelOK <- !is.na(maxFragMatches) && maxFragMatches >= minFrags
                if (!levelOK)
                {
                    if (is.na(maxFragMatches))
                        levelFailed <- "no fragments to match"
                    else
                        levelFailed <- sprintf("not enough fragments: %d/%d", maxFragMatches, minFrags)
                }
            }
            else if (type == "retention")
            {
                rtm <- getVal(val, type, "max")
                levelOK <- !is.null(candidateRTDev) && !is.na(candidateRTDev) && numLTE(abs(candidateRTDev), rtm)
                if (!levelOK)
                {
                    if (is.null(candidateRTDev) || is.na(candidateRTDev))
                        levelFailed <- "no retention time information available"
                    else
                        levelFailed <- sprintf("too high retention time deviation: %f/%f",
                                               abs(candidateRTDev), rtm)
                }
            }
            else if (type == "rank")
            {
                r <- if (getValType(val, type) == "formula") formRank else compRank
                maxR <- getVal(val, type, "max")
                levelOK <- !is.na(r) && r <= maxR
                if (!levelOK)
                    levelFailed <- if (is.na(r)) "candidate not ranked" else sprintf("ranked too low: %d/%d", r, maxR)
            }
            else if (type == "annMSMSSim")
            {
                sim <- if (getValType(val, type) == "formula") candidateAnnSimForm else candidateAnnSimComp
                minSim <- getVal(val, type, "min")
                levelOK <- !is.na(sim) && numGTE(sim, minSim)
                if (!levelOK)
                    levelFailed <- if (is.na(sim)) "no calculated similarity" else sprintf("similarity too low: %f/%f", sim, minSim)
            }
            else if (type == "annMSMSSimBoth")
            {
                minSim <- getVal(val, type, "min")
                levelOK <- !is.na(candidateAnnSimBoth) && numGTE(candidateAnnSimBoth, minSim)
                if (!levelOK)
                    levelFailed <- if (is.na(candidateAnnSimBoth)) "no calculated similarity" else
                        sprintf("similarity too low: %f/%f", candidateAnnSimBoth, minSim)
            }
            else if (type %in% allScores)
            {
                if (type %in% formCompScores)
                    isForm <- getValType(val, type) == "formula"
                else
                    isForm <- type %in% formScores
                
                if (isForm)
                    levelOK <- checkAnnotationScore(val, type, formRank, fRow, formTable, fRowNorm, formTableNorm,
                                                    getAllMergedConsCols(type, names(formTable), mFormNames))
                else
                    levelOK <- checkAnnotationScore(val, type, compRank, cRow, compTable, cRowNorm, compTableNorm,
                                                    getAllMergedConsCols(type, names(compTable), mCompNames))
                
                if (!isTRUE(levelOK))
                {
                    levelFailed <- levelOK
                    levelOK <- FALSE
                }
            }
            else
                stop(paste("Unknown ID level type:", type))
            
            if (!levelOK)
            {
                doLog(indent, "ID level failed: %s\n", levelFailed)
                return(FALSE)
            }
            doLog(indent, "ID level type passed!\n")
        }
        
        return(TRUE)
    }
    
    doLog(0, "Estimating identification level for '%s' to feature group '%s'\n---\n", candidateName, candidateFGroup)
    
    if (is.null(fRow))
    {
        if (is.null(candidateFormula))
            doLog(0, "NOTE: there is no formula data for available this candidate.\n")
        else
            doLog(0, "NOTE: the candidate formula could not be matched with the formula annotations.\n")
    }
    if (is.null(cRow))
    {
        if (is.null(candidateInChIKey1))
            doLog(0, "NOTE: there is no compound data (eg InChIKey) available for this candidate.\n")
        else
            doLog(0, "NOTE: the candidate compound could not be matched with the compound annotations.\n")
    }
    
    for (lvl in names(IDLevelRules))
    {
        doLog(0, "Checking level '%s'\n", lvl)
        if (checkLevelOK(IDLevelRules[[lvl]]))
        {
            doLog(0, "assigned level '%s'!\n", lvl)
            return(lvl)
        }
    }
    
    return(NA_character_)
}
