# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include workflow-step.R
NULL

#' Class containing MS Peak Lists
#'
#' Contains all MS (and MS/MS where available) peak lists for a \code{\link{featureGroups}} object.
#'
#' Objects for this class are returned by \code{\link{generateMSPeakLists}}.
#'
#' @param reAverage Set to \code{TRUE} to regenerate group averaged MS peak lists. \strong{NOTE} it is very important
#'   that any annotation data relying on MS peak lists (formulae/compounds) are regenerated afterwards! Otherwise it is
#'   likely that \emph{e.g.} plotting methods will use wrong MS/MS data.
#' @param MSLevel The MS level: \samp{1} for regular MS, \samp{2} for MSMS.
#' @param \dots Further arguments passed to \code{\link[graphics]{plot}}.
#' 
#'   \setsPassedArgs1{MSPeakLists}
#'
#' @template plot-lim
#' @template specSimParams-arg
#'
#' @slot peakLists Contains a list of all MS (and MS/MS) peak lists. Use the \code{peakLists} method for access.
#' @slot metadata Metadata for all spectra used to generate peak lists. Follows the format of the \code{peakLists} slot.
#' @slot averagedPeakLists A \code{list} with averaged MS (and MS/MS) peak lists for each feature group.
#' @slot avgPeakListArgs A \code{list} with arguments used to generate feature group averaged MS(/MS) peak lists.
#' @slot origFGNames A \code{character} with the original input feature group names.
#'
#' @templateVar seli analyses
#' @templateVar selOrderi analyses()
#' @templateVar selj feature groups
#' @templateVar selOrderj groupNames()
#' @templateVar optionalji TRUE
#' @templateVar del TRUE
#' @templateVar deli feature groups
#' @templateVar delj MS peaks
#' @templateVar deljtype numeric indices (rows)
#' @templateVar delfwhat feature group
#' @templateVar delfa the peak list table (a \code{data.table}), feature group name, analysis (\code{NULL} if averaged peak list), type (\code{"MS"} or \code{"MSMS"})
#' @templateVar delfr the peak list indices (rows) to be removed (specified as an \code{integer} or \code{logical} vector)
#' @templateVar dollarOpName feature group
#' @template sub_sel_del-args
#'
#' @section Isolating precursor data: Formula calculation typically relies on evaluating the measured isotopic pattern
#'   from the precursor to score candidates. Some algorithms (currently only \command{GenForm}) penalize candidates if
#'   mass peaks are present in MS1 spectra that do not contribute to the isotopic pattern. Since these spectra are
#'   typically very 'noisy' due to background and co-eluting ions, an additional filtering step may be recommended prior
#'   to formula calculation. During this precursor isolation step all mass peaks are removed that are (1) not the
#'   precursor and (2) not likely to be an isotopologue of the precursor. To determine potential isotopic peaks the
#'   following parameters are used:
#'
#'   \itemize{
#'
#'   \item \code{maxIsotopes} The maximum number of isotopes to consider. For instance, a value of \samp{5} means that
#'   \code{M+0} (\emph{i.e.} the monoisotopic peak) till \code{M+5} is considered. All mass peaks outside this range are
#'   removed.
#'
#'   \item \code{mzDefectRange} A two-sized \code{vector} specifying the minimum (can be negative) and maximum
#'   \emph{m/z} defect deviation compared to the precursor \emph{m/z} defect. When chlorinated, brominated or other
#'   compounds with strong \emph{m/z} defect in their isotopologues are to be considered a higher range may be desired.
#'   On the other hand, for natural compounds this range may be tightened. Note that the search range is propegated with
#'   increasing distance from the precursor, \emph{e.g.} the search range is doubled for \code{M+2}, tripled for
#'   \code{M+3} etc.
#'
#'   \item \code{intRange} A two-sized \code{vector} specifying the minimum and maximum relative intensity range
#'   compared to the precursor. For instance, \code{c(0.001, 2)} removes all peaks that have an intensity below 0.1\% or
#'   above 200\% of that of the precursor.
#'
#'   \item \code{z} The \code{z} value (\emph{i.e.} absolute charge) to be considerd. For instance, a value of \code{2}
#'   would look for \code{M+0.5}, \code{M+1} etc. Note that the \code{mzDefectRange} is adjusted accordingly
#'   (\emph{e.g.} halved if \code{z=2}).
#'
#'   \item \code{maxGap} The maximum number of missing adjacent isotopic peaks ('gaps'). If the (rounded) \emph{m/z}
#'   difference to the previous peak exceeds this value then this and all next peaks will be removed. Similar to
#'   \code{z}, the maximum gap is automatically adjusted for \code{charge}.
#'
#'   }
#'
#'   These parameters should be in a \code{list} that is passed to the \code{isolatePrec} argument to \code{filter}. The
#'   default values can be obtained with the \code{getDefIsolatePrecParams} function:
#'
#' @eval paste0("@@section Isolating precursor data:", getDefIsolatePrecParamsRD())
#'
#' @templateVar class MSPeakLists
#' @template class-hierarchy
#'
#' @param obj,x,object The \code{\link{MSPeakLists}} object to access.
#' @export
MSPeakLists <- setClass("MSPeakLists",
                        slots = c(peakLists = "list", metadata = "list", averagedPeakLists = "list",
                                  avgPeakListArgs = "list", origFGNames = "character"),
                        contains = "workflowStep")


setMethod("initialize", "MSPeakLists", function(.Object, setIDs = TRUE, doAverage = TRUE, ...)
{
    .Object <- callNextMethod(.Object, ...)

    # average if not already done (eg unset objects may do themselves)
    if (doAverage)
        .Object@averagedPeakLists <- averageMSPeakLists(.Object)
    .Object@peakLists <- makeEmptyListNamed(.Object@peakLists)
    .Object@averagedPeakLists <- makeEmptyListNamed(.Object@averagedPeakLists)
    
    if (setIDs)
        .Object <- assignMSPLIDs(.Object)

    return(.Object)
})

setMethod("averageMSPeakLists", "MSPeakLists", function(obj)
{
    # UNDONE: use cache sets?
    
    cat("Generating averaged peak lists for all feature groups...\n")
    
    pLists <- peakLists(obj)
    hash <- makeHash(pLists, obj@avgPeakListArgs)
    avgPLists <- loadCacheData("MSPeakListsAvg", hash)
    
    # figure out feature groups from (non-averaged) peak lists
    gNames <- unique(unlist(sapply(pLists, names, simplify = FALSE), use.names = FALSE))
    gNames <- intersect(obj@origFGNames, gNames) # sort to original order
    gCount <- length(gNames)
    
    if (gCount == 0 || length(obj@avgPeakListArgs) == 0)
        avgPLists <- makeEmptyListNamed(list())
    else if (is.null(avgPLists))
    {
        prog <- openProgBar(0, gCount)
        
        avgPLists <- lapply(seq_len(gCount), function(grpi)
        {
            plistsMS <- lapply(pLists, function(pl) pl[[gNames[grpi]]][["MS"]])
            plistsMS <- pruneList(plistsMS, checkZeroRows = TRUE)
            avgPLMS <- averageSpectra(plistsMS, obj@avgPeakListArgs$clusterMzWindow, obj@avgPeakListArgs$topMost,
                                      obj@avgPeakListArgs$minIntensityPre,
                                      obj@avgPeakListArgs$minIntensityPost,
                                      obj@avgPeakListArgs$avgFun, obj@avgPeakListArgs$method,
                                      obj@avgPeakListArgs$pruneMissingPrecursorMS, TRUE)
            
            plistsMSMS <- lapply(pLists, function(pl) pl[[gNames[grpi]]][["MSMS"]])
            plistsMSMS <- pruneList(plistsMSMS, checkZeroRows = TRUE)
            avgPLMSMS <- averageSpectra(plistsMSMS, obj@avgPeakListArgs$clusterMzWindow, obj@avgPeakListArgs$topMost,
                                        obj@avgPeakListArgs$minIntensityPre,
                                        obj@avgPeakListArgs$minIntensityPost,
                                        obj@avgPeakListArgs$avgFun, obj@avgPeakListArgs$method,
                                        FALSE, obj@avgPeakListArgs$retainPrecursorMSMS)
            
            results <- pruneList(list(MS = if (nrow(avgPLMS) > 0) avgPLMS else NULL,
                                      MSMS = if (nrow(avgPLMSMS) > 0) avgPLMSMS else NULL))
            
            if (!any(avgPLMS$precursor))
                warning(sprintf("Couldn't find back any precursor m/z from (averaged) MS peak list of group %s!", gNames[grpi]))
            
            setTxtProgressBar(prog, grpi)
            return(results)
        })
        names(avgPLists) <- gNames
        
		avgPLists <- pruneList(avgPLists, checkEmptyElements = TRUE)
		
        setTxtProgressBar(prog, gCount)
        close(prog)
        
        saveCacheData("MSPeakListsAvg", avgPLists, hash)
    }
    else
        cat("Done!\n")
    
    return(avgPLists)
})

#' @describeIn MSPeakLists Accessor method to obtain the MS peak lists.
#' @return \code{peakLists} returns a nested list containing MS (and MS/MS where
#'   available) peak lists per feature group and per analysis. The format is:
#'   \code{[[analysis]][[featureGroupName]][[MSType]][[PeakLists]]} where
#'   \code{MSType} is either \code{"MS"} or \code{"MSMS"} and \code{PeakLists} a
#'   \code{\link{data.table}} containing all \emph{m/z} values (\code{mz}
#'   column) and their intensities (\code{intensity} column). In addition, the
#'   peak list tables may contain a \code{cmp} column which contains an unique
#'   alphabetical identifier to which isotopic cluster (or "compound") a mass
#'   belongs (only supported by MS peak lists generated by Bruker tools at the
#'   moment).
#'
#' @aliases peakLists
#' @export
setMethod("peakLists", "MSPeakLists", function(obj) obj@peakLists)

#' @describeIn MSPeakLists Accessor method to obtain the feature group averaged
#'   MS peak lists.
#' @return \code{averagedPeakLists} returns a nested list of feature group
#'   averaged peak lists in a similar format as \code{peakLists}.
#' @aliases averagedPeakLists
#' @export
setMethod("averagedPeakLists", "MSPeakLists", function(obj) obj@averagedPeakLists)

#' @templateVar class MSPeakLists
#' @templateVar what analyses
#' @template strmethod
#' @export
setMethod("analyses", "MSPeakLists", function(obj) names(obj@peakLists))

#' @templateVar class MSPeakLists
#' @templateVar what feature groups
#' @template strmethod
#' @export
setMethod("groupNames", "MSPeakLists", function(obj) names(obj@averagedPeakLists))

#' @describeIn MSPeakLists Obtain total number of \emph{m/z} values.
#' @export
setMethod("length", "MSPeakLists", function(x) sum(unlist(recursiveApplyDT(x@peakLists, nrow))))

#' @describeIn MSPeakLists Shows summary information for this object.
#' @export
setMethod("show", "MSPeakLists", function(object)
{
    callNextMethod()

    if (length(object) == 0)
    {
        mspcount <- msmspcount <- totpcount <- totpcount <- mslcount <- msmslcount <- totlcount <- 0
        amspcount <- amsmspcount <- atotpcount <- atotpcount <- amslcount <- amsmslcount <- atotlcount <- 0
    }
    else
    {
        mspcount <- sum(sapply(object@peakLists, function(pa) sum(sapply(pa, function(pf) if (!is.null(pf[["MS"]])) nrow(pf[["MS"]]) else 0))))
        msmspcount <- sum(sapply(object@peakLists, function(pa) sum(sapply(pa, function(pf) if (!is.null(pf[["MSMS"]])) nrow(pf[["MSMS"]]) else 0))))
        totpcount <- sum(mspcount, msmspcount)
        mslcount <- sum(sapply(object@peakLists, function(pa) sum(sapply(pa, function(pf) !is.null(pf[["MS"]])))))
        msmslcount <- sum(sapply(object@peakLists, function(pa) sum(sapply(pa, function(pf) !is.null(pf[["MSMS"]])))))
        totlcount <- sum(mslcount, msmslcount)

        anacount <- length(object@peakLists)
        atotpcount <- totpcount / anacount
        amspcount <- mspcount / anacount
        amsmspcount <- msmspcount / anacount
        atotlcount <- totlcount / anacount
        amslcount <- mslcount / anacount
        amsmslcount <- msmslcount / anacount
    }

    printf("Total peak count: %d (MS: %d - MS/MS: %d)\n", totpcount, mspcount, msmspcount)
    printf("Average peak count/analysis: %.0f (MS: %.0f - MS/MS: %.0f)\n", atotpcount,
           amspcount, amsmspcount)
    printf("Total peak lists: %d (MS: %d - MS/MS: %d)\n", totlcount, mslcount, msmslcount)
    printf("Average peak lists/analysis: %.0f (MS: %.0f - MS/MS: %.0f)\n", atotlcount,
           amslcount, amsmslcount)
})

#' @describeIn MSPeakLists Subset on analyses/feature groups.
#' @export
setMethod("[", c("MSPeakLists", "ANY", "ANY", "missing"), function(x, i, j, ..., reAverage = FALSE, drop = TRUE)
{
    checkmate::assertFlag(reAverage)

    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, analyses(x))
        x <- delete(x, k = setdiff(analyses(x), i), reAverage = reAverage)
    }

    if (!missing(j))
    {
        j <- assertSubsetArgAndToChr(j, groupNames(x))
        x <- delete(x, i = setdiff(groupNames(x), j), reAverage = reAverage)
    }

    return(x)
})

#' @describeIn MSPeakLists Extract a list with MS and MS/MS (if available) peak
#'   lists. If the second argument (\code{j}) is not specified the averaged peak
#'   lists for the group specified by the first argument (\code{i}) will be
#'   returned.
#' @export
setMethod("[[", c("MSPeakLists", "ANY", "ANY"), function(x, i, j)
{
    assertExtractArg(i)
    if (!missing(j))
        assertExtractArg(j)

    if (!missing(j))
    {
        # both arguments specified, return regular peak lists

        if (!is.character(i))
            i <- analyses(x)[i]

        if (!is.character(j))
            j <- groupNames(x)[j]

        return(x@peakLists[[c(i, j)]])
    }

    # else return averaged peak lists for specified feature group

    if (!is.character(i))
        i <- groupNames(x)[i]

    return(x@averagedPeakLists[[i]])
})

#' @describeIn MSPeakLists Extract group averaged MS peaklists for a feature group.
#' @export
setMethod("$", "MSPeakLists", function(x, name)
{
    eval(substitute(x@averagedPeakLists[[NAME_ARG]], list(NAME_ARG = name)))
})

#' @describeIn MSPeakLists Returns all MS peak list data in a table.
#'
#' @param averaged If \code{TRUE} then feature group averaged peak list data is
#'   used.
#'
#' @template as_data_table-args
#'
#' @export
setMethod("as.data.table", "MSPeakLists", function(x, fGroups = NULL, averaged = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", null.ok = TRUE, add = ac)
    checkmate::assertFlag(averaged, add = ac)
    checkmate::reportAssertions(ac)

    if (averaged)
        ret <- rbindlist(lapply(averagedPeakLists(x), rbindlist, idcol = "type"), idcol = "group")
    else
    {
        ret <- rbindlist(lapply(peakLists(x), function(pa)
        {
            rbindlist(lapply(pa, rbindlist, idcol = "type"), idcol = "group")
        }), idcol = "analysis")
    }

    if (!is.null(fGroups))
    {
        ret[, c("ret", "group_mz") := groupInfo(fGroups)[group, c("rts", "mzs")]][]
        setcolorder(ret, c("group", "ret", "group_mz"))
        if (!averaged)
            setcolorder(ret, "analysis")
    }

    return(ret)
})

#' @templateVar where MSPeakLists
#' @templateVar what peaks from MS peak lists
#' @template delete
#' @param k A vector with analyses (\code{character} with names or \code{integer} with indices) for which the data
#'   should be deleted. If \code{k!=NULL} then deletions will \emph{not} occur on group averaged peak lists. Otherwise,
#'   if \code{k=NULL} then deltion occurs on \emph{both} group averaged and analysis specific peak lists.
#' @export
setMethod("delete", "MSPeakLists", function(obj, i = NULL, j = NULL, k = NULL, reAverage = FALSE, ...)
{
    doAnaOnly <- !is.null(k)
    
    ac <- checkmate::makeAssertCollection()
    i <- assertDeleteArgAndToChr(i, groupNames(obj), add = ac)
    checkmate::assert(
        checkmate::checkFunction(j, null.ok = TRUE),
        checkmate::checkIntegerish(j, null.ok = TRUE),
        .var.name = "j"
    )
    k <- assertDeleteArgAndToChr(k, analyses(obj), add = ac)
    checkmate::assertFlag(reAverage, add = ac)
    checkmate::reportAssertions(ac)
    
    if ((length(i) == 0 && length(k) == 0) || (!is.null(j) && length(j) == 0))
        return(obj) # nothing to remove...

    # i/j = NULL; k = vector: remove analyses (subset) (averagedPeakLists ignored)
    # i != NULL; j = NULL; k = vector: remove groups from specified analyses (averagedPeakLists ignored)
    # i = NULL; j != NULL; k = vector: remove from all groups from specified analyses (averagedPeakLists ignored)
    # i != NULL; j != NULL; k = vector: remove from specified groups from specified analyses (averagedPeakLists ignored)
    # i = vector; j = NULL; k = NULL: remove specified groups (peakLists and averagedPeakLists)
    # i = vector; j != NULL; k = NULL: remove from specified groups (peakLists and averagedPeakLists)
    # i = NULL; j != NULL; k = NULL: remove from all groups (peakLists and averagedPeakLists)

    if (!is.null(j) && !is.function(j))
    {
        origj <- j
        j <- function(pl, ...) origj
    }
    
    doDelPLs <- function(PL, grp, ana, ...)
    {
        for (t in c("MS", "MSMS"))
        {
            if (!is.null(PL[[t]]))
            {
                rm <- j(PL[[t]], grp, ana, t, ...)
                PL[[t]] <- if (is.logical(rm))
                    PL[[t]][!rm]
                else
                    PL[[t]][setdiff(seq_len(nrow(PL[[t]])), rm)]
            }
        }
        return(pruneList(PL, checkZeroRows = TRUE))
    }
    
    if (doAnaOnly && setequal(i, groupNames(obj)) && is.null(j))
    {
        # analyses subset
        obj@peakLists <- pruneList(obj@peakLists[setdiff(analyses(obj), k)])
    }
    else
    {
        if (is.null(j))
        {
            # remove results from groups completely
            keepFG <- setdiff(groupNames(obj), i)
            obj@peakLists <- lapply(obj@peakLists, function(x) return(pruneList(x[keepFG])))
            obj@peakLists <- pruneList(obj@peakLists, TRUE)
            
            if (!doAnaOnly)
                obj@averagedPeakLists <- pruneList(obj@averagedPeakLists[keepFG], TRUE)
        }
        else
        {
            obj@peakLists[k] <- Map(obj@peakLists[k], k, f = function(x, ana)
            {
                x[i] <- Map(x[i], i, f = doDelPLs, MoreArgs = list(ana, ...))
                return(pruneList(x, checkEmptyElements = TRUE))
            })

            if (!doAnaOnly)
            {
                obj@averagedPeakLists[i] <- Map(obj@averagedPeakLists[i], i, f = doDelPLs, MoreArgs = list(ana = NULL, ...))
                obj@averagedPeakLists <- pruneList(obj@averagedPeakLists, checkEmptyElements = TRUE)
            }
        }
        
    }

    if (reAverage)
    {
        obj@averagedPeakLists <- averageMSPeakLists(obj)
        obj <- assignMSPLIDs(obj) # re-generate as IDs are cleared by re-averaging
    }
    
    return(obj)
})

#' @describeIn MSPeakLists provides post filtering of generated MS peak lists, which may further enhance quality of
#'   subsequent workflow steps (\emph{e.g.} formulae calculation and compounds identification) and/or speed up these
#'   processes. The filters are applied to all peak lists for each analysis. These peak lists are subsequently averaged
#'   to update group averaged peak lists. However, since version \samp{1.1}, the resulting feature group lists are
#'   \emph{not} filtered afterwards.
#'
#' @param absMSIntThr,absMSMSIntThr,relMSIntThr,relMSMSIntThr Absolute/relative intensity threshold for MS or MS/MS peak
#'   lists. \code{NULL} for none.
#' @param topMSPeaks,topMSMSPeaks Only consider this amount of MS or MS/MS peaks with highest intensity. \code{NULL} to
#'   consider all.
#' @param minMSMSPeaks If the number of peaks in an MS/MS peak list (\strong{excluding} the precursor peak) is lower
#'   than this it will be completely removed. Set to \code{NULL} to ignore.
#' @param isolatePrec If not \code{NULL} then value should be a \code{list} with parameters used for isolating the
#'   precursor and its isotopes in MS peak lists (see \verb{Isolating precursor data}). Alternatively, \code{TRUE} to
#'   apply the filter with default settings (as given with \code{getDefIsolatePrecParams}).
#' @param deIsotopeMS,deIsotopeMSMS Remove any isotopic peaks in MS or MS/MS peak lists. This may improve data
#'   processing steps which do not assume the presence of isotopic peaks (e.g. MetFrag for MS/MS). Note that
#'   \code{getMzRPeakLists} does not (yet) support flagging of isotopes.
#' @param withMSMS If set to \code{TRUE} then only results will be retained for which MS/MS data is available. if
#'   \code{negate=TRUE} then only results \emph{without} MS/MS data will be retained.
#' @param annotatedBy Either a \code{\link{formulas}} or \code{\link{compounds}} object, or a \code{list} with both. Any
#'   MS/MS peaks that are \emph{not} annotated by any of the candidates in the specified objects are removed.
#' @param retainPrecursorMSMS If \code{TRUE} then precursor peaks will never be filtered out from MS/MS peak lists (note
#'   that precursors are never removed from MS peak lists). The \code{negate} argument does not affect this setting.
#' @param negate If \code{TRUE} then filters are applied in opposite manner.
#'
#' @export
setMethod("filter", "MSPeakLists", function(obj, absMSIntThr = NULL, absMSMSIntThr = NULL, relMSIntThr = NULL,
                                            relMSMSIntThr = NULL, topMSPeaks = NULL, topMSMSPeaks = NULL,
                                            minMSMSPeaks = NULL, isolatePrec = NULL, deIsotopeMS = FALSE,
                                            deIsotopeMSMS = FALSE, withMSMS = FALSE, annotatedBy = NULL,
                                            retainPrecursorMSMS = TRUE, reAverage = FALSE,
                                            negate = FALSE)
{
    if (is.logical(isolatePrec) && isolatePrec == TRUE)
        isolatePrec <- getDefIsolatePrecParams()

    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertNumber, . ~ absMSIntThr + absMSMSIntThr + relMSIntThr + relMSMSIntThr,
           lower = 0, finite = TRUE, null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertCount, . ~ topMSPeaks + topMSMSPeaks + minMSMSPeaks, positive = TRUE,
           null.ok = TRUE, fixed = list(add = ac))
    assertPListIsolatePrecParams(isolatePrec, add = ac)
    aapply(checkmate::assertFlag, . ~ deIsotopeMS + deIsotopeMSMS + withMSMS + retainPrecursorMSMS + reAverage + negate,
           fixed = list(add = ac))
    checkmate::assert(
        checkmate::checkNull(annotatedBy),
        checkmate::checkClass(annotatedBy, "formulas"),
        checkmate::checkClass(annotatedBy, "compounds"),
        checkmate::checkList(annotatedBy, c("formulas", "compounds"), any.missing = FALSE, min.len = 1, unique = TRUE),
        .var.name = "annotatedBy"
    )
    checkmate::reportAssertions(ac)

    if (length(obj) == 0)
        return(obj)

    hash <- makeHash(obj, absMSIntThr, absMSMSIntThr, relMSIntThr, relMSMSIntThr, topMSPeaks, topMSMSPeaks, minMSMSPeaks,
                     isolatePrec, deIsotopeMS, deIsotopeMSMS, withMSMS, annotatedBy, retainPrecursorMSMS, reAverage,
                     negate)
    cache <- loadCacheData("filterMSPeakLists", hash)
    if (!is.null(cache))
        return(cache)
    
    if (!is.null(annotatedBy) && !is.list(annotatedBy))
        annotatedBy <- list(annotatedBy)
    
    oldn <- length(obj)
    
    obj <- delete(obj, j = function(pl, grp, ana, type)
    {
        plF <- copy(pl)
        
        isReAveraged <- is.null(ana) && reAverage
        
        if (type == "MS")
        {
            plF <- doMSPeakListFilter(plF, absMSIntThr, relMSIntThr, topMSPeaks, NULL, deIsotopeMS, TRUE, negate)
            if (!is.null(isolatePrec))
                plF <- isolatePrecInMSPeakList(plF, isolatePrec, negate)
        }
        else
        {
            if (!isReAveraged && !is.null(annotatedBy))
            {
                allAnnPLIDs <- unique(unlist(lapply(annotatedBy, function(ab)
                {
                    if (is.null(ab[[grp]]))
                        return(integer())
                    return(unlist(lapply(ab[[grp]]$fragInfo, "[[", "PLID")))
                })))
                
                if (length(allAnnPLIDs) == 0)
                {
                    if (!negate)
                        plF <- plF[retainPrecursorMSMS & precursor]
                }
                else
                {
                    if (negate)
                        plF[, keep := !ID %in% allAnnPLIDs]
                    else
                        plF[, keep := ID %in% allAnnPLIDs]
                    
                    if (retainPrecursorMSMS)
                        plF[precursor == TRUE, keep := TRUE]
                    plF <- plF[keep == TRUE][, keep := NULL]
                }
            }

            plF <- doMSPeakListFilter(plF, absMSMSIntThr, relMSMSIntThr, topMSMSPeaks, minMSMSPeaks, deIsotopeMSMS,
                                      retainPrecursorMSMS, negate)
        }
        
        return(!pl$ID %in% plF$ID)
    }, reAverage = reAverage)

    # apply after other filters as these may have removed all MSMS peaks
    if (withMSMS)
    {
        obj <- delete(obj, j = function(pl, grp, ana, type)
        {
            keep <- if (type == "MSMS")
                TRUE
            else if (!is.null(ana))
                !is.null(obj[[ana, grp]][["MSMS"]])
            else
                !is.null(obj[[grp]][["MSMS"]])
            
            if (negate)
                keep <- !keep
            
            return(if (keep) integer() else seq_len(nrow(pl)))
        }, reAverage = reAverage)
    }
    
    saveCacheData("filterMSPeakLists", obj, hash)

    newn <- length(obj)
    printf("Done! Filtered %d (%.2f%%) MS peaks. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)

    return(obj)
})

#' @describeIn MSPeakLists Plots a spectrum using MS or MS/MS peak lists for a given feature group. Two spectra can be
#'   compared when two feature groups are specified.
#'
#' @param groupName The name of the feature group for which a plot should be made. To compare spectra, two group names
#'   can be specified.
#' @param analysis The name of the analysis for which a plot should be made. If \code{NULL} then data from the feature
#'   group averaged peak list is used. When comparing spectra, either \code{NULL} or the analyses for both spectra
#'   should be specified.
#' @param title The title of the plot. If \code{NULL} a title will be automatically made.
#'
#' @template plot-lim
#'
#' @export
setMethod("plotSpectrum", "MSPeakLists", function(obj, groupName, analysis = NULL, MSLevel = 1, title = NULL,
                                                  specSimParams = getDefSpecSimParams(),
                                                  xlim = NULL, ylim = NULL, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCharacter(groupName, min.len = 1, max.len = 2, min.chars = 1, add = ac)
    checkmate::assertCharacter(analysis, min.len = 1, max.len = 2, min.chars = 1, null.ok = TRUE, add = ac)
    if (!is.null(analysis) && length(analysis) != length(groupName))
        stop("Lengths of analysis and groupName should be equal.")
    assertSpecSimParams(specSimParams, add = ac)
    checkmate::assertChoice(MSLevel, 1:2, add = ac)
    assertXYLim(xlim, ylim, add = ac)
    checkmate::reportAssertions(ac)

    if (length(obj) == 0)
        return(NULL)

    setTitle <- is.null(title)
    if (setTitle)
        title <- getMSPeakListPlotTitle(MSLevel, analysis, groupName)
    
    if (length(groupName) == 1)
    {
        spec <- getSpec(obj, groupName, MSLevel, analysis)
        if (is.null(spec))
            return(NULL)
        
        makeMSPlot(getMSPlotData(spec, 2), 1, xlim, ylim, main = title, ...)
    }
    else
    {
        if (setTitle)
        {
            sim <- spectrumSimilarity(obj, groupName[1], groupName[2], analysis[1], analysis[2], MSLevel, specSimParams,
                                      NAToZero = TRUE, drop = TRUE)
            title <- c(title, sprintf("Similarity: %.2f", sim))
        }
        
        binnedPLs <- getBinnedPLPair(obj, groupName, analysis, MSLevel, specSimParams, "unique", mustExist = TRUE)
        plotData <- getMSPlotDataOverlay(binnedPLs, TRUE, FALSE, 2, "overlap")
        makeMSPlotOverlay(plotData, title, 1, xlim, ylim, ...)
    }
})

setMethod("plotSpectrumHash", "MSPeakLists", function(obj, groupName, analysis = NULL, MSLevel = 1, title = NULL,
                                                      specSimParams = getDefSpecSimParams(), ...)
{
    if (length(groupName) == 1)
        sp <- getSpec(obj, groupName, MSLevel, analysis)
    else
        sp <- c(getSpec(obj, groupName[1], MSLevel, analysis[1]),
                getSpec(obj, groupName[2], MSLevel, analysis[2]))
    return(makeHash(sp, title, specSimParams, ...))
})

#' @describeIn MSPeakLists Calculates the spectral similarity between two or more spectra.
#'
#' @param groupName1,groupName2 The names of the feature groups for which the comparison should be made. If both
#'   arguments are specified then a comparison is made with the spectra specified by \code{groupName1} \emph{vs} those
#'   specified by \code{groupName2}. The length of either can be \samp{>1} to generate a comparison matrix.
#'   Alternatively, if \code{groupName2} is \code{NULL} then all the spectra specified in \code{groupName1} will be
#'   compared with eachother, \emph{i.e.} resulting in a square similarity matrix.
#' @param analysis1,analysis2 The name of the analysis (analyses) for the comparison. If \code{NULL} then data from the
#'   feature group averaged peak list is used. Otherwise, should be the same length as
#'   \code{groupName1}/\code{groupName2}.
#' @param NAToZero Set to \code{TRUE} to convert \code{NA} similarities (\emph{i.e.} when no similarity could be
#'   calculated) to zero values.
#' @param drop If set to \code{TRUE} and if the comparison is made between two spectra then \code{\link{drop}} is used
#'   to reduce the \code{matrix} return value to a \code{numeric} vector.
#'
#' @author For \code{spectrumSimilarity}: major contributions by Bas van de Velde for spectral binning and similarity
#'   calculation.
#'   
#' @section Source: \code{spectrumSimilarity}: The principles of spectral binning and cosine similarity calculations
#'   were loosely was based on the code from \code{SpectrumSimilarity()} function of \pkg{OrgMassSpecR}.
#'
#' @aliases spectrumSimilarity
#' @export
setMethod("spectrumSimilarity", "MSPeakLists", function(obj, groupName1, groupName2 = NULL, analysis1 = NULL,
                                                        analysis2 = NULL, MSLevel = 1,
                                                        specSimParams = getDefSpecSimParams(),
                                                        NAToZero = FALSE, drop = TRUE)
{
    # NOTE: keep args in sync with sets method

    if (length(obj) == 0)
        return(NULL)
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertSubset, . ~ groupName1 + groupName2, empty.ok = c(FALSE, TRUE),
           fixed = list(choices = groupNames(obj), add = ac))
    aapply(checkmate::assertSubset, . ~ analysis1 + analysis2, empty.ok = TRUE,
           fixed = list(choices = analyses(obj), add = ac))
    checkmate::assertChoice(MSLevel, 1:2, add = ac)
    assertSpecSimParams(specSimParams, add = ac)
    aapply(checkmate::assertFlag, . ~ NAToZero + drop, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (is.null(groupName2))
    {
        # calculate dist matrix
        PLP <- getSimPLAndPrec(obj, groupName1, analysis1, MSLevel, specSimParams, 1)
        if (is.null(PLP))
        {
            sims <- matrix(NA_real_, length(groupName1), length(groupName1))
            rownames(sims) <- colnames(sims) <- names(groupName1)
        }
        else
        {
            sims <- specDistMatrix(PLP$specs, specSimParams$method, specSimParams$shift, PLP$precs,
                                   specSimParams$mzWeight, specSimParams$intWeight, specSimParams$absMzDev)
            rownames(sims) <- colnames(sims) <- names(PLP$specs)
            sims <- expandFillSpecSimilarities(sims, groupName1, groupName1)
        }
    }
    else
    {
        PLP1 <- getSimPLAndPrec(obj, groupName1, analysis1, MSLevel, specSimParams, 1)
        PLP2 <- getSimPLAndPrec(obj, groupName2, analysis2, MSLevel, specSimParams, 2)
        if (is.null(PLP1) || is.null(PLP2))
        {
            sims <- matrix(NA_real_, length(groupName1), length(groupName2))
            rownames(sims) <- groupName1; colnames(sims) <- groupName2
        }
        else
        {
            sims <- specDistRect(PLP1$specs, PLP2$specs, specSimParams$method, specSimParams$shift, PLP1$precs,
                                 PLP2$precs, specSimParams$mzWeight, specSimParams$intWeight, specSimParams$absMzDev)
            rownames(sims) <- names(PLP1$specs); colnames(sims) <- names(PLP2$specs)
            sims <- expandFillSpecSimilarities(sims, groupName1, groupName2)
        }
    }
    
    if (NAToZero)
        sims[is.na(sims)] <- 0
    
    return(if (drop && length(sims) == 1) drop(sims) else sims)
})

#' Generation of MS Peak Lists
#'
#' Functionality to convert MS and MS/MS data into MS peak lists.
#'
#' Formula calculation and identification tools rely on mass spectra that belong to features of interest. For
#' processing, MS (and MS/MS) spectra are typically reduced to a table with a column containing measured \emph{m/z}
#' values and a column containing their intensities. These 'MS peak lists' can then be used for
#' \link[=generateFormulas]{formula generation} and \link[=generateCompounds]{compound generation}.
#'
#' MS and MS/MS peak lists are first generated for all features (or a subset, if the \code{topMost} argument is set).
#' During this step multiple spectra over the feature elution profile are averaged. Subsequently, peak lists will be
#' generated for each feature group by averaging peak lists of the features within the group. Functionality that uses
#' peak lists will either use data from individual features or from group averaged peak lists. For instance, the former
#' may be used by formulae calculation, while compound identification and plotting functionality typically uses group
#' averaged peak lists.
#'
#' @templateVar func generateMSPeakLists
#' @templateVar what generate MS peak lists
#' @templateVar ex1 generateMSPeakListsMzR
#' @templateVar ex2 generateMSPeakListsDA
#' @templateVar algos bruker,brukerfmf,mzr
#' @templateVar algosSuffix DA,DAFMF,MzR
#' @templateVar ret MSPeakLists
#' @template generic-algo
#'
#' @param fGroups The \code{\link{featureGroups}} object from which MS peak lists should be extracted.
#' @param \dots Any parameters to be passed to the selected MS peak lists generation algorithm.
#'
#' @template centroid_note
#'
#' @return A \code{\link{MSPeakLists}} object.
#'
#' @section Sets workflows: With a \link[=sets-workflow]{sets workflow}, the feature group averaged peak lists are made
#'   per set. This is important, because for averaging peak lists cannot be mixed, for instance, when different
#'   ionization modes were used to generate the sets. The group averaged peaklists are then simply combined and labelled
#'   in the final peak lists. However, please note that annotation and other functionality typically uses only the set
#'   specific peak lists, as this functionality cannot work with mixed peak lists.
#'
#' @templateVar what generateMSPeakLists
#' @template main-rd-method
#' @export
setMethod("generateMSPeakLists", "featureGroups", function(fGroups, algorithm, ...)
{
    checkmate::assertChoice(algorithm, c("bruker", "brukerfmf", "mzr"))
    
    f <- switch(algorithm,
                bruker = generateMSPeakListsDA,
                brukerfmf = generateMSPeakListsDAFMF,
                mzr = generateMSPeakListsMzR)

    f(fGroups, ...)
})

