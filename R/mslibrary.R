#' @include main.R
#' @include workflow-step.R
NULL

#' @export
MSLibrary <- setClass("MSLibrary", slots = c(records = "data.table", spectra = "list"),
                      contains = "workflowStep")

setMethod("initialize", "MSLibrary", function(.Object, ...)
{
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})

#' @export
setMethod("records", "MSLibrary", function(obj) obj@records)

#' @export
setMethod("spectra", "MSLibrary", function(obj) obj@spectra)

#' @export
setMethod("length", "MSLibrary", function(x) nrow(records(x)))

#' @export
setMethod("names", "MSLibrary", function(x) names(spectra(x)))

#' @export
setMethod("show", "MSLibrary", function(object)
{
    callNextMethod()
    
    # UNDONE: more?
    printf("Total records: %d\n", length(object))
})

#' @export
setMethod("[", c("MSLibrary", "ANY", "missing", "missing"), function(x, i, ...)
{
    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, names(x))
        x <- delete(x, setdiff(names(x), i))
    }
    return(x)
})

#' @export
setMethod("[[", c("MSLibrary", "ANY", "missing"), function(x, i, j)
{
    assertExtractArg(i)
    return(x@spectra[[i]])
})

#' @export
setMethod("$", "MSLibrary", function(x, name)
{
    eval(substitute(x@spectra$NAME_ARG, list(NAME_ARG = name)))
})

#' @export
setMethod("as.data.table", "MSLibrary", function(x)
{
    allSpecs <- rbindlist(spectra(x), idcol = "DB_ID", fill = TRUE)
    return(merge(records(x), allSpecs, by = "DB_ID"))
})

#' @export
setMethod("delete", "MSLibrary", function(obj, i = NULL, j = NULL, ...)
{
    if (!is.function(i))
        i <- assertDeleteArgAndToChr(i, names(obj))
    checkmate::assert(
        checkmate::checkFunction(j, null.ok = TRUE),
        checkmate::checkIntegerish(j, any.missing = FALSE, null.ok = TRUE),
        .var.name = "j"
    )
    if (is.function(i) && is.function(j))
        stop("i and j cannot both be functions")

    if (length(i) == 0 || (!is.null(j) && length(j) == 0))
        return(obj) # nothing to remove...
    
    # i = vector remove specified records
    # i = function: remove returned records
    # j = vector/function: index of peaks to remove from spectra of i (all if i=NULL)
    
    if (is.null(j))
    {
        if (!is.function(i))
            obj@records <- obj@records[!DB_ID %chin% i]
        else
        {
            rm <- i(obj@records, ...)
            if (is.logical(rm))
                obj@records <- obj@records[!rm]
            else if (is.character(rm))
                obj@records <- obj@records[!rm %chin% DB_ID]
            else
                obj@records <- obj@records[setdiff(seq_len(nrow(obj@records)), rm)]
        }
        
        obj@spectra <- obj@spectra[obj@records$DB_ID]
    }
    else
    {
        for (rec in i)
        {
            if (is.function(j))
            {
                inds <- j(rec, obj@spectra[[rec]])
                if (is.logical(inds))
                    inds <- which(inds)
            }
            else # j = vector
                inds <- j[j <= nrow(obj@spectra)]
            if (length(inds) > 0)
                obj@spectra[[rec]] <- obj@spectra[[rec]][-inds, , drop = FALSE]
        }
        obj@spectra <- pruneList(obj@spectra, checkZeroRows = TRUE)
        obj@records <- obj@records[DB_ID %chin% names(obj@spectra)]
    }
    
    return(obj)
})

#' @export
setMethod("filter", "MSLibrary", function(obj, properties = NULL, massRange = NULL, mzRangeSpec = NULL,
                                          relMinIntensity = NULL, topMost = NULL, onlyAnnotated = FALSE,
                                          negate = FALSE)
{
    
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(properties, any.missing = FALSE, null.ok = TRUE, add = ac)
    if (!is.null(properties))
    {
        checkmate::assertNames(names(properties), type = "unique", subset.of = names(records(obj)), add = add)
        checkmate::qassertr(properties, "V")
    }
    aapply(assertRange, . ~ massRange + mzRangeSpec, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertNumber(relMinIntensity, null.ok = TRUE, add = ac)
    checkmate::assertCount(topMost, positive = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertFlag(onlyAnnotated, add = ac)
    checkmate::assertFlag(negate, add = ac)
    checkmate::reportAssertions(ac)
    
    mark <- if (negate) function(x) !x else function(x) x

    if (onlyAnnotated)
    {
        noAnnot <- sapply(spectra(obj), function(sp) is.null(sp[["annotation"]]))
        obj <- delete(obj, i = if (negate) !noAnnot else noAnnot)
    }
    
    if (!is.null(properties) || !is.null(massRange))
    {
        obj <- delete(obj, i = function(recs, ...)
        {
            recs <- copy(recs)
            recs[, keep := TRUE]
            
            if (!is.null(properties) && length(properties) > 0)
            {
                for (prop in names(properties))
                    recs[keep == TRUE, keep := mark(get(prop) %in% properties[[prop]])]
            }
            
            if (!is.null(massRange))
                recs[keep == TRUE, keep := !is.na(neutralMass) & mark(neutralMass %inrange% massRange)]
            
            return(!recs$keep)
        })
    }
    if (!is.null(mzRangeSpec) || !is.null(relMinIntensity))
    {
        obj <- delete(obj, j = function(rec, spec, ...)
        {
            keep <- rep(TRUE, nrow(spec))
            if (!is.null(mzRangeSpec))
                keep <- mark(spec[, "mz"] %inrange% mzRangeSpec)
            if (!is.null(relMinIntensity))
            {
                relInts <- spec[, "intensity"] / max(spec[, "intensity"])
                keep <- keep & mark(numGTE(relInts, relMinIntensity))
            }
            return(!keep)
        })
    }
    if (!is.null(topMost)) # NOTE: do after previous filters
    {
        if (!is.null(mzRangeSpec) || !is.null(relMinIntensity) || !is.null(topMost))
        obj <- delete(obj, j = function(rec, spec, ...)
        {
            err <- tryCatch(rep(TRUE, nrow(spec)), error=function(...) NULL)
            if (is.null(err)) browser()
            keep <- rep(TRUE, nrow(spec))
            if (nrow(spec) > topMost)
            {
                ord <- order(spec[, "intensity"], decreasing = !negate)
                keep <- seq_len(nrow(spec)) %in% ord[seq_len(topMost)] # NOTE: keep order
            }
            return(!keep)
        })
    }
    
    return(obj)
})

#' @export
setMethod("convertToSuspects", "MSLibrary", function(obj, adduct,
                                                     avgSpecParams = getDefAvgPListParams(minIntensityPre = 0,
                                                                                          minIntensityPost = 2,
                                                                                          topMost = 10),
                                                     collapse = TRUE, suspects = NULL)
{
    adduct <- checkAndToAdduct(adduct)
    
    ac <- checkmate::makeAssertCollection()
    assertAvgPListParams(avgSpecParams, add = ac)
    checkmate::assertFlag(collapse, add = ac)
    if (!is.null(suspects))
        assertSuspectList(suspects, FALSE, FALSE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
        stop("Cannot create suspect list: no data", call. = FALSE)

    libRecs <- records(obj)
    libRecs <- libRecs[Precursor_type == as.character(adduct)]
    if (nrow(libRecs) == 0)
        stop("No records found (for input adduct)", call. = FALSE)
    libSpecs <- spectra(obj)
    
    getAvgFrags <- function(specs, PrecursorMZ)
    {
        pls <- Map(specs, PrecursorMZ, f = function(sp, pmz)
        {
            sp <- copy(sp)[, c("mz", "intensity"), with = FALSE] # skip annotation, if present
            sp <- assignPrecursorToMSPeakList(sp, pmz)
            return(sp)
        })
        
        avgPL <- averageSpectra(pls, avgSpecParams$clusterMzWindow, avgSpecParams$topMost,
                                avgSpecParams$minIntensityPre, avgSpecParams$minIntensityPost,
                                avgSpecParams$avgFun, avgSpecParams$method, FALSE, avgSpecParams$retainPrecursorMSMS)
        doProgress()
        paste0(avgPL$mz, collapse = ";")
    }
        
    printf("Calculating MS/MS fragments...\n")
    
    if (!is.null(suspects))
    {
        ret <- if (is.data.table(suspects)) copy(suspects) else as.data.table(suspects)
        ret <- prepareSuspectList(ret, NULL, FALSE, calcMZs = FALSE)
        ret[, InChIKey1 := getIKBlock1(InChIKey)]
        
        if (any(is.na(ret$InChIKey1)))
            warning(paste("Ignored the following suspects because no InChIKey1 could be calculated:",
                           getStrListWithMax(ret$name, 10, ", ")))
        
        recs <- copy(libRecs)
        recs <- recs[!is.na(InChIKey) & !is.na(PrecursorMZ)]
        recs[, InChIKey1 := getIKBlock1(InChIKey)]
        
        withProg(uniqueN(ret$InChIKey1), FALSE, ret[, fragments_mz := sapply(InChIKey1, function(IK1)
        {
            recsSub <- recs[InChIKey1 == IK1]
            if (is.na(IK1) || nrow(recsSub) == 0)
                return("")
            return(getAvgFrags(libSpecs[recsSub$DB_ID], recsSub$PrecursorMZ))
        }), by = "InChIKey1"])
        
        ret[, InChIKey1 := NULL]
        
        frCount <- sum(nzchar(ret$fragments_mz))
        printf("Filled in fragments for %d/%d (%.2f%%) suspects\n", frCount, nrow(ret), frCount / nrow(ret) * 100)
    }
    else
    {
        ret <- copy(libRecs)
        
        if (collapse)
        {
            ret <- ret[!is.na(InChIKey) & !is.na(PrecursorMZ)]
            ret[, InChIKey1 := getIKBlock1(InChIKey)]

            frMZ <- withProg(uniqueN(ret$InChIKey1), FALSE, ret[, getAvgFrags(libSpecs[DB_ID], PrecursorMZ),
                                                                by = "InChIKey1"])
            setnames(frMZ, "V1", "fragments_mz")

            ret <- unique(ret, by = "InChIKey1")
            ret <- merge(ret, frMZ, by = "InChIKey1")
        }
        else
        {
            withProg(nrow(ret), FALSE, ret[, fragments_mz := sapply(libSpecs[ret$DB_ID], function(spec)
            {
                doProgress()
                paste0(spec[, "mz"], collapse = ";")
            })])
        }
        
        mapCols <- c(Name = "name",
                     SMILES = "SMILES",
                     InChI = "InChI",
                     InChIKey = "InChIKey",
                     formula = "formula",
                     Precursor_type = "adduct",
                     neutralMass = "neutralMass",
                     fragments_mz = "fragments_mz")
        mapCols <- mapCols[names(mapCols) %in% names(ret)]
        setnames(ret, names(mapCols), mapCols)
        ret <- ret[, mapCols, with = FALSE]
        ret <- prepareSuspectList(ret, NULL, FALSE, FALSE)
    }
    
    return(ret[])
})

setMethod("export", "MSLibrary", function(obj, type, out)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertChoice(type, "msp", add = ac)
    checkmate::assertPathForOutput(out, overwrite = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    # convert records to character matrix to simplfy Rcpp processing
    recs <- copy(records(obj))
    recs[, (names(recs)) := lapply(.SD, as.character)]
    recs <- as.matrix(recs)
    writeMSPLibrary(recs, spectra(obj), normalizePath(out))
})

setMethod("merge", c("MSLibrary", "MSLibrary"), function(x, y, ...)
{
    if (length(x) == 0)
        return(y)
    else if (length(y) == 0)
        return(x)
    
    if (any(is.na(records(x)$SPLASH)) || any(is.na(records(y)$SPLASH)))
        stop("x/y has missing SPLASH values. Please load the library with calcSPLASH=TRUE")
    
    unY <- records(y)[!SPLASH %chin% records(x)$SPLASH]$DB_ID
    y <- y[unY]
    
    recordsAll <- rbind(records(x), records(y), fill = TRUE)
    specsAll <- c(spectra(x), spectra(y))
    
    recordsAll <- makeDBIDsUnique(recordsAll)
    names(specsAll) <- recordsAll$DB_ID
    
    return(MSLibrary(records = recordsAll[], spectra = specsAll, algorithm = "merged"))
})


loadMSLibrary <- function(file, algorithm, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFileExists(file, "r", add = ac)
    checkmate::assertChoice(algorithm, c("msp", "json"), add = ac)
    checkmate::reportAssertions(ac)
    
    f <- switch(algorithm,
                msp = loadMSLibraryMSP,
                json = loadMSLibraryMoNAJSON,
                xcms = findFeaturesXCMS)
    
    f(file, ...)
}
