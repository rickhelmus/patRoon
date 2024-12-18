# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include workflow-step.R
NULL

#' Class to store data from a loaded MS library
#'
#' Stores the spectra and metadata from the records of an MS library.
#'
#' This class is used by \code{\link{loadMSLibrary}} to store the loaded MS library data.
#'
#' @param x,obj,object \code{MSLibrary} object to be accessed.
#'
#' @seealso \code{\link{loadMSLibrary}}
#'
#' @slot records A \code{\link{data.table}} with metadata for all records. Use the \code{records} method for access.
#' @slot spectra A \code{list} with all (annotated) spectra. Each spectrum is stored in a \code{\link{data.table}}. Use
#'   the \code{spectra} method for access.
#'
#' @templateVar seli records
#' @templateVar selOrderi names()
#' @templateVar dollarOpName record
#' @template sub_sel_del-args
#'
#' @templateVar class MSLibrary
#' @template class-hierarchy
#' 
#' @references \insertRef{Wohlgemuth2016}{patRoon} \cr\cr
#'   \addCitations{Rcpp}{1} \cr\cr \addCitations{Rcpp}{2} \cr\cr \addCitations{Rcpp}{3}
#'
#' @export
MSLibrary <- setClass("MSLibrary", slots = c(records = "data.table", spectra = "list"),
                      contains = "workflowStep")

setMethod("initialize", "MSLibrary", function(.Object, ...)
{
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})

#' @describeIn MSLibrary Accessor method for the \code{records} slot of an \code{MSLibrary} class.
#' @aliases records
#' @export
setMethod("records", "MSLibrary", function(obj) obj@records)

#' @describeIn MSLibrary Accessor method for the \code{spectra} slot of an \code{MSLibrary} class.
#' @aliases spectra
#' @export
setMethod("spectra", "MSLibrary", function(obj) obj@spectra)

#' @describeIn MSLibrary Obtains the total number of records stored.
#' @export
setMethod("length", "MSLibrary", function(x) nrow(records(x)))

#' @describeIn MSLibrary Obtains the names of the stored records (\code{DB_ID} field).
#' @export
setMethod("names", "MSLibrary", function(x) names(spectra(x)))

#' @describeIn MSLibrary Shows summary information for this object.
#' @export
setMethod("show", "MSLibrary", function(object)
{
    callNextMethod()
    
    printf("Total records: %d\n", length(object))
    
    totPeaks <- if (length(object) > 0) sum(sapply(spectra(object), nrow)) else 0
    totPeaksAnn <- if (length(object) > 0) sum(sapply(spectra(object), function(sp)
    {
        if (is.null(sp[["annotation"]]))
            return(0)
        return(sp[!is.na(annotation) & nzchar(annotation), .N])
    }))
    else
        0
    
    printf("Total peaks: %d\n", totPeaks)
    printf("Total annotated peaks: %d (%.2f%%)\n", totPeaksAnn, if (totPeaks > 0) totPeaksAnn * 100 / totPeaks else 0)
})

#' @describeIn MSLibrary Subset on records.
#' @param \dots Unused.
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

#' @describeIn MSLibrary Extracts a spectrum table for a record.
#' @export
setMethod("[[", c("MSLibrary", "ANY", "missing"), function(x, i, j)
{
    assertExtractArg(i)
    return(x@spectra[[i]])
})

#' @describeIn MSLibrary Extracts a spectrum table for a record.
#' @export
setMethod("$", "MSLibrary", function(x, name)
{
    eval(substitute(x@spectra$NAME_ARG, list(NAME_ARG = name)))
})

#' @describeIn MSLibrary Converts all the data (spectra and metadata) to a single \code{data.table}.
#' @export
setMethod("as.data.table", "MSLibrary", function(x)
{
    if (length(x) == 0)
        return(records(x))
    allSpecs <- rbindlist(spectra(x), idcol = "DB_ID", fill = TRUE)
    return(merge(records(x), allSpecs, by = "DB_ID"))
})

#' @templateVar where MSLibrary
#' @templateVar what full records or spectra
#' @template delete
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
        obj@spectra[i] <- Map(obj@spectra[i], i, f = function(sp, ind)
        {
            if (is.function(j))
            {
                rm <- j(ind, sp)
                if (is.logical(rm))
                    return(sp[!rm])
                return(sp[setdiff(seq_len(nrow(sp)), rm)])
            }
            
            # j = vector
            inds <- j[j <= nrow(sp)]
            return(if (length(inds) > 0) sp[-inds] else sp)
        })
        obj@spectra <- pruneList(obj@spectra, checkZeroRows = TRUE)
        obj@records <- obj@records[DB_ID %chin% names(obj@spectra)]
    }
    
    return(obj)
})

#' @describeIn MSLibrary Performs rule-based filtering of records and spectra. This may be especially to improve
#'   annotation with \code{\link{generateCompoundsLibrary}}.
#'
#' @param massRange Records with a neutral mass outside this range will be removed. Should be a two-sized \code{numeric}
#'   vector with the lower and upper mass range. Set to \code{NULL} to ignore.
#' @param mzRangeSpec Similar to the \code{massRange} argument, but removes any peaks from recorded mass spectra outside
#'   the given \emph{m/z} range.
#' @param relMinIntensity The minimum relative intensity (\samp{0-1}) of a mass peak to be kept. Set to \code{NULL} to
#'   ignore.
#' @param topMost Only keep \code{topMost} number of mass peaks for each spectrum. This filter is applied after others.
#'   Set to \code{NULL} to ignore.
#' @param onlyAnnotated If \code{TRUE} then only recorded spectra that are formula annotated are kept.
#' @param negate If \code{TRUE} then filters are performed in opposite manner.
#'
#' @templateVar getProps names(records)
#' @templateVar ex Instrument_type=c("LC-ESI-QTOF","LC-ESI-TOF")
#' @template filter-properties
#'
#' @return \code{filter} returns a filtered \code{MSLibrary} object.
#'
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
                keep <- mark(spec$mz %inrange% mzRangeSpec)
            if (!is.null(relMinIntensity))
            {
                relInts <- spec$intensity / max(spec$intensity)
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
            keep <- rep(TRUE, nrow(spec))
            if (nrow(spec) > topMost)
            {
                ord <- order(spec$intensity, decreasing = !negate)
                keep <- seq_len(nrow(spec)) %in% ord[seq_len(topMost)] # NOTE: keep order
            }
            return(!keep)
        })
    }
    
    return(obj)
})

#' @describeIn MSLibrary Converts the MS library data to a suspect list, which can be used with
#'   \code{\link{screenSuspects}}. See the \verb{Suspect conversion} section for details.
#'
#' @param adduct An \code{\link{adduct}} object (or something that can be converted to it with \code{\link{as.adduct}}).
#'   Any records with a different adduct (\code{Precursor_type}) are not considered. Alternatively, \code{adduct} can be
#'   set to \code{NULL} to not filter out any records. However, in this case \emph{no} MS/MS fragments will be added to
#'   the returned suspect list.
#' @param avgSpecParams A \code{list} with parameters used for averaging spectra. See \code{\link{getDefAvgPListParams}}
#'   for more details.
#' @param collapse Whether records with the same first-block \acronym{InChIKey} should be collapsed. See the
#'   \verb{Suspect conversion} section for details.
#' @param suspects If not \code{NULL} then this should be a suspect list (see \code{\link{screenSuspects}}) which will
#'   be amended with spectra data. See the \verb{Suspect conversion} section for details.
#'
#' @template spectrumType-arg
#'
#' @return \code{convertToSuspects} return a suspect list (\code{data.table}), which can be used with
#'   \code{\link{screenSuspects}}.
#'
#' @section Suspect conversion: The \code{convertToSuspects} method converts MS library data to a suspect list, which
#'   can be used with \emph{e.g.} \code{\link{screenSuspects}}. Furthermore, this function can also amend existing
#'   suspect lists with spectral data.
#'
#'   Conversion occurs in either of the following three methods: \enumerate{
#'
#'   \item \emph{Direct} (\code{collapse=FALSE} and \code{suspects=NULL}): each record is considered a suspect, and the
#'   resulting suspect list is generated directly by converting the records metadata. The \code{fragments_mz} column for
#'   each suspect is constructed from the mass peaks of the corresponding record.
#'
#'   \item \emph{Collapse} (\code{collapse=TRUE} and \code{suspects=NULL}): All records with the same first-block
#'   \acronym{InChIKey} are first merged, and their spectra are averaged using the parameters from the
#'   \code{avgSpecParams} argument (see \code{\link{getDefAvgPListParams}}). The suspect list is based on the merged
#'   records, where the \code{fragments_mz} column is constructed from the averaged spectra. This is generally a good
#'   default, especially with large MS libraries.
#'
#'   \item \emph{Amend} (\code{suspects} is not \code{NULL}): only those records are considered if their first-block
#'   \acronym{InChIKey} is present in the suspect list. The remaining records and their spectra are then collapsed as
#'   described for the \emph{Collapse} method, and the \code{fragments_mz} column for each suspect is set from the
#'   averaged spectra. If a suspect is not present in the library, its \code{fragments_mz} value will be empty. Note
#'   that any existing \code{fragments_mz} data will be overwritten.
#'
#'   }
#'
#' @templateVar whatCP input suspect list to \code{convertToSuspects}
#' @template chemPropCalc
#'
#' @export
setMethod("convertToSuspects", "MSLibrary", function(obj, adduct, spectrumType = "MS2",
                                                     avgSpecParams = getDefAvgPListParams(minIntensityPre = 0,
                                                                                          minIntensityPost = 2,
                                                                                          topMost = 10),
                                                     collapse = TRUE, suspects = NULL, prefCalcChemProps = TRUE,
                                                     neutralChemProps = FALSE)
{
    if (!is.null(adduct))
        adduct <- checkAndToAdduct(adduct)
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCharacter(spectrumType, min.len = 1, min.chars = 1, null.ok = TRUE, add = ac)
    assertAvgPListParams(avgSpecParams, add = ac)
    aapply(checkmate::assertFlag, . ~ prefCalcChemProps + neutralChemProps + collapse, fixed = list(add = ac))
    if (!is.null(suspects))
        assertSuspectList(suspects, FALSE, FALSE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
        stop("Cannot create suspect list: no data", call. = FALSE)
    
    calcMSMS <- !is.null(adduct)
    if (!is.null(suspects) && !calcMSMS)
    {
        warning("No adduct specified, nothing to do...", call. = FALSE)
        return(suspects)
    }

    hash <- makeHash(obj, adduct, spectrumType, avgSpecParams, collapse, suspects, prefCalcChemProps, neutralChemProps)
    cd <- loadCacheData("convertToSuspectsMSLibrary", hash)
    if (!is.null(cd))
        return(cd)
    
    libRecs <- records(obj)
    
    if (!is.null(adduct))
        libRecs <- libRecs[Precursor_type == as.character(adduct)]
    if (!is.null(spectrumType))
        libRecs <- libRecs[Spectrum_type %chin% spectrumType]
    if (nrow(libRecs) == 0)
        stop("No records found (for input adduct/spectrum type)", call. = FALSE)
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

    if (calcMSMS)   
        printf("Calculating MS/MS fragments...\n")
    
    if (!is.null(suspects))
    {
        # NOTE: we always calculate MS/MS fragments here, since an early exit is done above if calcMSMS==FALSE
        
        ret <- if (is.data.table(suspects)) copy(suspects) else as.data.table(suspects)
        # checkDesc = TRUE: we want to be able to calculate InChIKey1 values
        ret <- prepareSuspectList(ret, NULL, FALSE, checkDesc = TRUE, prefCalcChemProps = prefCalcChemProps,
                                  neutralChemProps = neutralChemProps, calcMZs = FALSE)
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
            
            if (calcMSMS)
            {
                frMZ <- withProg(uniqueN(ret$InChIKey1), FALSE, ret[, getAvgFrags(libSpecs[DB_ID], PrecursorMZ),
                                                                    by = "InChIKey1"])
                setnames(frMZ, "V1", "fragments_mz")
                
                ret <- unique(ret, by = "InChIKey1")
                ret <- merge(ret, frMZ, by = "InChIKey1")
            }
            else
                ret <- unique(ret, by = "InChIKey1")
        }
        else if (calcMSMS)
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
                     neutralMass = "neutralMass",
                     molNeutralized = "molNeutralized",
                     fragments_mz = "fragments_mz")
        if (!is.null(adduct))
            mapCols <- c(mapCols, Precursor_type = "adduct")
        mapCols <- mapCols[names(mapCols) %in% names(ret)]
        setnames(ret, names(mapCols), mapCols)
        ret <- ret[, mapCols, with = FALSE]
        ret <- prepareSuspectList(ret, NULL, FALSE, FALSE, FALSE, FALSE)
    }
    
    saveCacheData("convertToSuspectsMSLibrary", ret, hash)
    
    return(ret[])
})

#' @describeIn MSLibrary Exports the library data to a \file{.msp} file. The export is accelerated by an \code{C++}
#'   interface with \CRANpkg{Rcpp}.
#'
#' @param type The export type. Currently just \code{"msp"}.
#' @param out The file path to the output library file.
#'
#' @note \code{export} does not split any \code{Synon} data that was merged when the library was loaded.
#' 
#' @export
setMethod("export", "MSLibrary", function(obj, type = "msp", out)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertChoice(type, "msp", add = ac)
    checkmate::assertPathForOutput(out, overwrite = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    # convert records to character matrix to simplify Rcpp processing
    recs <- copy(records(obj))
    recs[, (names(recs)) := lapply(.SD, as.character)]
    recs <- as.matrix(recs)
    writeMSPLibrary(recs, spectra(obj), normalizePath(out, mustWork = FALSE))
})

#' @describeIn MSLibrary Merges two \code{MSLibrary} objects (\code{x} and \code{y}). The records from \code{y} that are
#'   unique are added to \code{x}. Records that were already in \code{x} are simply ignored. The
#'   \href{https://splash.fiehnlab.ucdavis.edu/}{SPLASH} values are used to test equality between records, hence, the
#'   \code{calcSPLASH} argument to \code{\link{loadMSLibrary}} should be \code{TRUE}.
#'
#' @param y The \code{MSLibrary} to be merged with \code{x}.
#'
#' @return \code{merge} returns a merged \code{MSLibrary} object.
#'
#' @export
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

#' Loading of MS library data
#'
#' Loads, parses, verifies and curates MS library data, \emph{e.g.} obtained from MassBank.
#'
#' @templateVar func loadMSLibrary
#' @templateVar what loads MS library data
#' @templateVar ex1 loadMSLibraryMSP
#' @templateVar ex2 loadMSLibraryMoNAJSON
#' @templateVar algos msp,json
#' @templateVar algosSuffix MSP,MoNAJSON
#' @templateVar ret MSLibrary
#' @template generic-algo
#'
#' @param file A \code{character} string that specifies the the file path to the library.
#' @param \dots Any parameters to be passed to the selected MS library loading algorithm.
#'
#' @return A \code{\link{MSLibrary}} object containing the loaded library data.
#'
#' @export
loadMSLibrary <- function(file, algorithm, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFileExists(file, "r", add = ac)
    checkmate::assertChoice(algorithm, c("msp", "json"), add = ac)
    checkmate::reportAssertions(ac)
    
    f <- switch(algorithm,
                msp = loadMSLibraryMSP,
                json = loadMSLibraryMoNAJSON)
    
    f(file, ...)
}
