#' @include main.R
#' @include workflow-step.R
NULL

#' Class containing MS Peak Lists
#'
#' Contains all MS (and MS/MS where available) peak lists for a
#' \code{\link{featureGroups}} object.
#'
#' Objects for this class are returned by \link[=MSPeakLists-generation]{MS peak
#' lists generators}.
#'
#' @slot peakLists Contains a list of all MS (and MS/MS) peak lists. Use the
#'   \code{peakLists} method for access.
#' @slot metadata Metadata for all spectra used to generate peak lists. Follows
#'   the format of the \code{peakLists} slot.
#' @slot averagedPeakLists A \code{list} with averaged MS (and MS/MS) peak lists
#'   for each feature group.
#' @slot avgPeakListArgs A \code{list} with arguments used to generate feature
#'   group averaged MS(/MS) peak lists.
#' @slot origFGNames A \code{character} with the original input feature group
#'   names.
#'
#' @templateVar seli analyses
#' @templateVar selOrderi analyses()
#' @templateVar selj feature groups
#' @templateVar selOrderj groupNames()
#' @templateVar optionalji TRUE
#' @templateVar dollarOpName feature group
#' @template sub_op-args
#'
#' @section Isolating precursor data: Formula calculation typically relies on
#'   evaluating the measured isotopic pattern from the precursor to score
#'   candidates. Some algorithms (currently only \command{GenForm}) penalize
#'   candidates if mass peaks are present in MS1 spectra that do not contribute
#'   to the isotopic pattern. Since these spectra are typically very 'noisy' due
#'   to background and co-eluting ions, an additional filtering step may be
#'   recommended prior to formula calculation. During this precursor isolation
#'   step all mass peaks are removed that are (1) not the precursor and (2) not
#'   likely to be an isotopologue of the precursor. To determine potential
#'   isotopic peaks the following parameters are used:
#'
#'   \itemize{
#'
#'   \item \code{maxIsotopes} The maximum number of isotopes to consider. For
#'   instance, a value of \samp{5} means that \code{M+0} (\emph{i.e.} the
#'   monoisotopic peak) till \code{M+5} is considered. All mass peaks outside
#'   this range are removed.
#'
#'   \item \code{mzDefectRange} A two-sized \code{vector} specifying the minimum
#'   (can be negative) and maximum \emph{m/z} defect deviation compared to the
#'   precursor \emph{m/z} defect. When chlorinated, brominated or other
#'   compounds with strong \emph{m/z} defect in their isotopologues are to be
#'   considered a higher range may be desired. On the other hand, for natural
#'   compounds this range may be tightened. Note that the search range is
#'   propegated with increasing distance from the precursor, \emph{e.g.} the
#'   search range is doubled for \code{M+2}, tripled for \code{M+3} etc.
#'
#'   \item \code{intRange} A two-sized \code{vector} specifying the minimum and
#'   maximum relative intensity range compared to the precursor. For instance,
#'   \code{c(0.001, 2)} removes all peaks that have an intensity below 0.1\% or
#'   above 200\% of that of the precursor.
#'
#'   \item \code{z} The \code{z} value (\emph{i.e.} absolute charge) to be
#'   considerd. For instance, a value of \code{2} would look for \code{M+0.5},
#'   \code{M+1} etc. Note that the \code{mzDefectRange} is adjusted accordingly
#'   (\emph{e.g.} halved if \code{z=2}).
#'
#'   \item \code{maxGap} The maximum number of missing adjacent isotopic peaks
#'   ('gaps'). If the (rounded) \emph{m/z} difference to the previous peak
#'   exceeds this value then this and all next peaks will be removed. Similar to
#'   \code{z}, the maximum gap is automatically adjusted for \code{charge}.
#'
#'   }
#'
#'   These parameters should be in a \code{list} that is passed to the
#'   \code{isolatePrec} argument to \code{filter}. The default values can be
#'   obtained with the \code{getDefIsolatePrecParams} function:
#'
#' @eval paste0("@@section Isolating precursor data:",
#'   getDefIsolatePrecParamsRD())
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


setMethod("initialize", "MSPeakLists", function(.Object, ...)
{
    .Object <- callNextMethod(.Object, ...)

    .Object@averagedPeakLists <- averageMSPeakLists(.Object)
    .Object@peakLists <- makeEmptyListNamed(.Object@peakLists)
    .Object@averagedPeakLists <- makeEmptyListNamed(.Object@averagedPeakLists)
    
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
        avgPLists <- list()
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
#' @param reAverage Set to \code{TRUE} to regenerate averaged MS peak lists
#'   after subsetting analyses.
#' @export
setMethod("[", c("MSPeakLists", "ANY", "ANY", "missing"), function(x, i, j, ..., reAverage = TRUE, drop = TRUE)
{
    checkmate::assertFlag(reAverage)

    # non-existing indices result in NULL values --> prune

    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, analyses(x))
        x@peakLists <- pruneList(x@peakLists[i])

        # update group averaged peak lists
        if (reAverage)
            x@averagedPeakLists <- averageMSPeakLists(x)
    }

    if (!missing(j))
    {
        j <- assertSubsetArgAndToChr(j, groupNames(x))
        x@peakLists <- sapply(x@peakLists, function(a) return(pruneList(a[j])),
                              simplify = FALSE)
        x@peakLists <- pruneList(x@peakLists, TRUE)

        x@averagedPeakLists <- pruneList(x@averagedPeakLists[j], TRUE)
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

#' @describeIn MSPeakLists provides post filtering of generated MS peak lists,
#'   which may further enhance quality of subsequent workflow steps (\emph{e.g.}
#'   formulae calculation and compounds identification) and/or speed up these
#'   processes. The filters are applied to all peak lists for each analysis.
#'   These peak lists are subsequently averaged to update group averaged peak
#'   lists. However, since version \samp{1.1}, the resulting feature group lists are
#'   \emph{not} filtered afterwards.
#'   
#' @param absMSIntThr,absMSMSIntThr,relMSIntThr,relMSMSIntThr Absolute/relative
#'   intensity threshold for MS or MS/MS peak lists. \code{NULL} for none.
#' @param topMSPeaks,topMSMSPeaks Only consider this amount of MS or MS/MS peaks
#'   with highest intensity. \code{NULL} to consider all.
#' @param isolatePrec If not \code{NULL} then value should be a \code{list} with
#'   parameters used for isolating the precursor and its isotopes in MS peak
#'   lists (see \verb{Isolating precursor data}). Alternatively, \code{TRUE} to
#'   apply the filter with default settings (as given with
#'   \code{getDefIsolatePrecParams}).
#' @param deIsotopeMS,deIsotopeMSMS Remove any isotopic peaks in MS or MS/MS
#'   peak lists. This may improve data processing steps which do not assume the
#'   presence of isotopic peaks (e.g. MetFrag for MS/MS). Note that
#'   \code{getMzRPeakLists} does not (yet) support flagging of isotopes.
#' @param withMSMS If set to \code{TRUE} then only results will be retained for
#'   which MS/MS data is available. if \code{negate=TRUE} then only results
#'   \emph{without} MS/MS data will be retained.
#' @param retainPrecursorMSMS If \code{TRUE} then precursor peaks will never be
#'   filtered out from MS/MS peak lists (note that precursors are never removed
#'   from MS peak lists). The \code{negate} argument does not affect this
#'   setting.
#' @param negate If \code{TRUE} then filters are applied in opposite manner.
#'
#' @export
setMethod("filter", "MSPeakLists", function(obj, absMSIntThr = NULL, absMSMSIntThr = NULL, relMSIntThr = NULL,
                                            relMSMSIntThr = NULL, topMSPeaks = NULL, topMSMSPeaks = NULL,
                                            isolatePrec = NULL, deIsotopeMS = FALSE, deIsotopeMSMS = FALSE,
                                            withMSMS = FALSE, retainPrecursorMSMS = TRUE, negate = FALSE)
{
    if (is.logical(isolatePrec) && isolatePrec == TRUE)
        isolatePrec <- getDefIsolatePrecParams()

    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertNumber, . ~ absMSIntThr + absMSMSIntThr + relMSIntThr + relMSMSIntThr,
           lower = 0, finite = TRUE, null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertCount, . ~ topMSPeaks + topMSMSPeaks, positive = TRUE,
           null.ok = TRUE, fixed = list(add = ac))
    assertPListIsolatePrecParams(isolatePrec, add = ac)
    aapply(checkmate::assertFlag, . ~ deIsotopeMS + deIsotopeMSMS + withMSMS + retainPrecursorMSMS + negate, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    if (length(obj) == 0)
        return(obj)

    hash <- makeHash(obj, absMSIntThr, absMSMSIntThr, relMSIntThr, relMSMSIntThr,
                     topMSPeaks, topMSMSPeaks, isolatePrec, deIsotopeMS, deIsotopeMSMS,
                     withMSMS, retainPrecursorMSMS, negate)
    cache <- loadCacheData("filterMSPeakLists", hash)
    if (!is.null(cache))
        return(cache)

    doFilterGroups <- function(pl)
    {
        gCount <- length(pl)
        if (gCount == 0)
            return(pl)

        prog <- openProgBar(0, gCount)

        pln <- names(pl)
        pl <- lapply(seq_along(pl), function(grpi)
        {
            if (withMSMS)
            {
                if (!negate && is.null(pl[[grpi]][["MSMS"]]))
                    return(list())
                if (negate && !is.null(pl[[grpi]][["MSMS"]]))
                    return(list())
            }

            ret <- list()
            if (!is.null(pl[[grpi]][["MS"]]))
                ret$MS <- doMSPeakListFilter(pl[[grpi]]$MS, absMSIntThr, relMSIntThr, topMSPeaks, deIsotopeMS, TRUE, negate)
            if (!is.null(pl[[grpi]][["MSMS"]]))
                ret$MSMS <- doMSPeakListFilter(pl[[grpi]]$MSMS, absMSMSIntThr, relMSMSIntThr, topMSMSPeaks,
                                               deIsotopeMSMS, retainPrecursorMSMS, negate)

            if (!is.null(isolatePrec))
                ret$MS <- isolatePrecInMSPeakList(ret$MS, isolatePrec, negate)

            setTxtProgressBar(prog, grpi)

            return(pruneList(ret, checkZeroRows = TRUE))
        })
        names(pl) <- pln

        pl <- pl[sapply(pl, function(p) !is.null(p[["MS"]]) || !is.null(p[["MSMS"]]))]

        setTxtProgressBar(prog, gCount)
        close(prog)

        return(pl)
    }

    pLists <- peakLists(obj)
    oldn <- length(obj)

    for (anai in seq_along(pLists))
    {
        printf("Filtering MS peak lists for %d feature groups in analysis '%s'...\n", length(pLists[[anai]]),
               names(pLists)[anai])
        pLists[[anai]] <- doFilterGroups(pLists[[anai]])
    }

    obj@peakLists <- pLists
    
    # update group averaged peak lists
    obj@averagedPeakLists <- averageMSPeakLists(obj)

    saveCacheData("filterMSPeakLists", obj, hash)

    newn <- length(obj)
    printf("Done! Filtered %d (%.2f%%) MS peaks. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)

    return(obj)
})

#' @describeIn MSPeakLists Plots a spectrum using MS or MS/MS peak lists for a
#'   given feature group.
#'
#' @param groupName The name of the feature group for which a plot should be
#'   made.
#' @param analysis The name of the analysis for which a plot should be made. If
#'   \code{NULL} then data from the feature group averaged peak list is used.
#' @param MSLevel The MS level: \samp{1} for regular MS, \samp{2} for MSMS.
#' @param title The title of the plot. If \code{NULL} a title will be
#'   automatically made.
#' @param \dots Further arguments passed to \code{\link[graphics]{plot}}.
#'
#' @template useGGplot2
#'
#' @template plot-lim
#'
#' @return \code{plotSpectrum} will return a \code{\link[=ggplot2]{ggplot
#'   object}} if \code{useGGPlot2} is \code{TRUE}.
#'
#' @export
setMethod("plotSpectrum", "MSPeakLists", function(obj, groupName, analysis = NULL, MSLevel = 1, title = NULL,
                                                  useGGPlot2 = FALSE, xlim = NULL, ylim = NULL, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertString(groupName, min.chars = 1, add = ac)
    checkmate::assertString(analysis, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertChoice(MSLevel, 1:2, add = ac)
    checkmate::assertFlag(useGGPlot2, add = ac)
    assertXYLim(xlim, ylim, add = ac)
    checkmate::reportAssertions(ac)

    if (length(obj) == 0)
        return(NULL)

    spec <- getSpec(obj, groupName, MSLevel, analysis)
    if (is.null(spec))
        return(NULL)

    if (is.null(title))
        title <- getMSPeakListPlotTitle(MSLevel, analysis, groupName)
    
    if (useGGPlot2)
        return(makeMSPlotGG(getMSPlotData(spec, 2)) + ggtitle(title))

    makeMSPlot(getMSPlotData(spec, 2), xlim, ylim, main = title, ...)
})


#' @templateVar func generateMSPeakLists
#' @templateVar what generate MS peak lists
#' @templateVar ex1 generateMSPeakListsMzR
#' @templateVar ex2 generateMSPeakListsDA
#' @templateVar algos bruker,brukerfmf,mzr
#' @template generic-algo
#'
#' @rdname MSPeakLists-generation
#' @aliases generateMSPeakLists
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

