# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

emptyMSPeakList <- function(abundanceColumn, avgCols)
{
    ret <- data.table(mz = numeric(), intensity = numeric())
    if (!is.null(abundanceColumn))
        ret[, (paste0(abundanceColumn, c("_rel", "_abs"))) := numeric()]
    if (!is.null(avgCols))
        ret[, (avgCols) := numeric()]
    ret[, precursor := logical()]
    return(ret)
}
    

#' Parameters for averaging MS peak list data
#'
#' Create parameter lists for averaging MS peak list data.
#'
#' The parameters set used for averaging peak lists are set by the \code{avgFeatParams} and \code{avgFGroupParams}
#' arguments to \code{\link{generateMSPeakLists}} and its related algorithm specific functions. The parameters are
#' specified as a named \code{list} with the following values: \itemize{
#'
#' \item \code{method},\code{clusterMzWindow} The cluster method and window (see \link[=cluster-params]{clustering
#' parameters}) used to average mass spectra. \code{clusterMzWindow} is defaulted as \code{defaultLim("mz", "medium")} (see
#' \link{limits}).
#'
#' \item \code{topMost} Only retain this maximum number of MS peaks when generating averaged spectra. Lowering this
#' number may exclude more irrelevant (noisy) MS peaks and decrease processing time, whereas higher values may avoid
#' excluding lower intense MS peaks that may still be of interest.
#'
#' \item \code{minIntensityPre} MS peaks with intensities below this value will be removed (applied prior to selection
#' by \code{topMost} and averaging).
#'
#' \item \code{minIntensityPost} MS peaks with intensities below this value will be removed (after averaging).
#'
#' \item \code{minIntensityIMS} MS peaks in spectra of raw IMS frames with intensities below this value will be removed
#' (applied prior to any other treatment steps).
#'
#' \item \code{minAbundanceAbs},\code{minAbundanceRel} Minimum absolute/relative abundance of an MS peak across the
#' spectra that are averaged. If \code{minAbundanceAbs} exceeds the number of spectra then the threshold is
#' automatically lowered to the number of spectra.
#'
#' \item \code{minAbundanceIMSAbs},\code{minAbundanceIMSRel} Minimum absolute/relative abundance of an MS peak across
#' the spectra that are summed within an IMS frame. If \code{minAbundanceIMSAbs} exceeds the number of spectra in the
#' IMS frame then the threshold is automatically lowered to the number of spectra.
#'
#' \item \code{withPrecursorMS} For MS data only: ignore any spectra that do not contain the precursor peak.
#'
#'   For IMS data this excludes MS spectra within an IMS frame that do not contain the precursor peak, typically due to
#'   mobility separation. Hence, setting this option performs some crude cleanup of MS spectra, even for features for
#'   which no mobilities were assigned (\emph{e.g.} non-IMS workflows).
#'
#' \item \code{pruneMissingPrecursorMS} For MS data only: if \code{TRUE} then peak lists without a precursor peak are
#' removed. Note that even when this is set to \code{FALSE}, functionality that relies on MS (not MS/MS) peak lists
#' (\emph{e.g.} formulae calculation) will still skip calculation if a precursor is not found.
#'
#' \item \code{retainPrecursorMSMS} For MS/MS data only: if \code{TRUE} then always retain the precursor mass peak even
#' if is not amongst the \code{topMost} peaks. Note that MS precursor mass peaks are always kept. Furthermore, note that
#' precursor peaks in both MS and MS/MS data may still be removed by intensity thresholds (this is unlike the
#' \code{\link[=filter,MSPeakLists-method]{filter}} method function).
#'
#' }
#'
#' The \code{getDefAvgPListParams} function can be used to generate a default parameter list. The current defaults are:
#'
#' @eval paste("@@details", getDefAvgPListParamsRD())
#'
#' @param \dots Optional named arguments that override defaults.
#'
#' @return \code{getDefAvgPListParams} returns a \code{list} with the peak list averaging parameters.
#'
#' @export
getDefAvgPListParams <- function(...)
{
    def <- list(clusterMzWindow = defaultLim("mz", "medium"),
                topMost = 50,
                minIntensityPre = 500,
                minIntensityPost = 500,
                minAbundanceAbs = 0,
                minAbundanceRel = 0,
                smoothWindowIMS = 0,
                halfWindowIMS = 2,
                minIntensityIMS = 25,
                maxGapIMS = defaultLim("mz", "medium"),
                method = "distance_mean",
                withPrecursorMS = TRUE,
                pruneMissingPrecursorMS = TRUE,
                retainPrecursorMSMS = TRUE)
    return(modifyList(def, list(...)))
}

# For docs
# nocov start
getDefAvgPListParamsRD <- function()
{
    def <- getDefAvgPListParams()
    def <- sapply(def, function(v) if (is.character(v)) paste0("\"", v, "\"") else v)
    return(paste0("\\code{", names(def), "=", def, "}", collapse = "; "))
}
# nocov end

#' @details The \code{getDefIsolatePrecParams} is used to create a parameter
#'   list for isolating the precursor and its isotopes (see \verb{Isolating precursor data}).
#' @rdname MSPeakLists-class
#' @export
getDefIsolatePrecParams <- function(...)
{
    def <- list(maxIsotopes = 5,
                mzDefectRange = c(-0.01, 0.01),
                intRange = c(0.001, 2),
                z = 1,
                maxGap = 2)
    return(modifyList(def, list(...)))
}

# For docs
# nocov start
getDefIsolatePrecParamsRD <- function()
{
    def <- getDefIsolatePrecParams()
    def <- sapply(def, function(v) if (length(v) == 2) sprintf("c(%s, %s)", v[1], v[2]) else as.character(v))
    return(paste0("\\code{", names(def), "=", def, "}", collapse = "; "))
}
# nocov end

#' MS spectral similarity calculation parameters
#'
#' Parameters relevant for calculation of similarities between mass spectra.
#'
#' For the calculation of spectral similarities the following parameters exist:
#'
#' \itemize{
#'
#' \item \code{method} The similarity method: either \code{"cosine"} or \code{"jaccard"}.
#'
#' \item \code{removePrecursor} If \code{TRUE} then precursor peaks (\emph{i.e.} the mass peak corresponding to the
#' feature) are removed prior to similarity calculation.
#'
#' \item \code{mzWeight},\code{intWeight} Mass and intensity weights used for cosine calculation.
#'
#' \item \code{absMzDev} Maximum absolute \emph{m/z} deviation between mass peaks, used for binning spectra. Defaults to
#' \code{defaultLim("mz", "medium")} (see \link{limits}).
#'
#' \item \code{relMinIntensity} The minimum relative intensity for mass peaks (\samp{0-1}). Peaks with lower intensities
#' are not considered for similarity calculation. The relative intensities are called after the precursor peak is
#' removed when \code{removePrecursor=TRUE}.
#'
#' \item \code{minPeaks} Only consider spectra that have at least this amount of peaks (\emph{after} the spectrum is
#' filtered).
#'
#' \item \code{shift} If and how shifting is applied prior to similarity calculation. Valid options are: \code{"none"}
#' (no shifting), \code{"precursor"} (all mass peaks of the second spectrum are shifted by the mass difference between
#' the precursors of both spectra) or \code{"both"} (the spectra are first binned without shifting, and peaks still
#' unaligned are then shifted as is done when \code{shift="precursor"}).
#'
#' \item \code{setCombinedMethod} \setsWF Determines how spectral similarities from different sets are combined.
#' Possible values are \code{"mean"}, \code{"min"} or \code{"max"}, which calculates the combined value as the mean,
#' minimum or maximum value, respectively. \code{NA} values (\emph{e.g.} if a set does not have peak list data to
#' combine) are removed in advance.
#'
#' }
#'
#' These parameters are typically passed as a named \code{list} as the \code{specSimParams} argument to functions that
#' do spectral similarity calculations. The \code{getDefSpecSimParams} function can be used to generate such parameter
#' list with defaults.
#'
#' @param \dots optional named arguments that override defaults.
#'
#' @name specSimParams
#' @export
getDefSpecSimParams <- function(...)
{
    def <- list(method = "cosine",
                removePrecursor = FALSE,
                mzWeight = 0,
                intWeight = 1,
                absMzDev = defaultLim("mz", "medium"),
                relMinIntensity = 0.05,
                minPeaks = 1,
                shift = "none",
                setCombineMethod = "mean")

    return(modifyList(def, list(...)))
}

averageSpectraList <- function(spectraList, clusterMzWindow, topMost, minIntensityPre, minIntensityPost,
                               minAbundanceRel, minAbundanceAbs, method, assignPrecursor, withPrecursor,
                               pruneMissingPrecursor, retainPrecursor)
{
    # pre-treat
    spectraList <- lapply(spectraList, function(spectra)
    {
        lapply(spectra, function(s)
        {
            if (withPrecursor && !any(s$precursor))
                return(s[0])
            
            s <- s[intensity >= minIntensityPre]
            if (nrow(s) > topMost)
            {
                ord <- order(s$intensity, decreasing = TRUE)
                keep <- ord[seq_len(topMost)]
                if (retainPrecursor)
                    keep <- union(keep, which(s$precursor))
                s <- s[keep]
            }
            return(s)
        })
    })

    averagedSpectra <- doAverageSpectraList(spectraList, method, clusterMzWindow, minIntensityPost, minAbundanceRel,
                                            minAbundanceAbs)

    # post-treat
    averagedSpectra <- lapply(seq_along(averagedSpectra), function(i)
    {
        avgSpec <- averagedSpectra[[i]]
        setDT(avgSpec)
        
        if (nrow(avgSpec) == 0)
            return(avgSpec)
        
        if (assignPrecursor)
        {
            spcomb <- rbindlist(spectraList[[i]])
            precMZs <- spcomb[precursor == TRUE]$mz
            if (length(precMZs) > 0)
                avgSpec <- assignPrecursorToMSPeakList(avgSpec, mean(precMZs))
            else
                avgSpec[, precursor := FALSE]
        }
        else
            avgSpec[, precursor := FALSE]
        
        if (pruneMissingPrecursor && !any(avgSpec$precursor))
            avgSpec <- avgSpec[0]
        
        return(avgSpec)
    })
        
    return(averagedSpectra)
}

# get corresponding mz of feature from MS peaklist
getMZIndexFromMSPeakList <- function(featMZ, plist)
{
    ret <- which.min(abs(plist$mz - featMZ))
    if (abs(plist$mz[ret] - featMZ) <= 0.01) # UNDONE: make range configurable?
        return(ret)
    return(NA)
}

assignPrecursorToMSPeakList <- function(MSPeakList, precursorMZ)
{
    if (nrow(MSPeakList) == 0)
        MSPeakList[, precursor := logical()]
    else
    {
        MSPeakList[, precursor := FALSE]
        ind <- getMZIndexFromMSPeakList(precursorMZ, MSPeakList)
        if (!is.na(ind))
            MSPeakList[ind, precursor := TRUE]
    }
    return(MSPeakList)
}

assignMSPLIDs <- function(MSPeakLists)
{
    addIDs <- function(pl)
    {
        add <- function(p)
        {
            p <- copy(p)
            p[, ID := seq_len(.N)]
            setcolorder(p, "ID")
            return(p)
        }
        if (!is.null(pl[["MS"]]))
            pl$MS <- add(pl$MS)
        if (!is.null(pl[["MSMS"]]))
            pl$MSMS <- add(pl$MSMS)
        
        return(pl)
    }
    MSPeakLists@peakLists <- lapply(MSPeakLists@peakLists, function(pla) lapply(pla, addIDs))
    MSPeakLists@averagedPeakLists <- lapply(MSPeakLists@averagedPeakLists, addIDs)
    return(MSPeakLists)
}

deIsotopeMSPeakList <- function(MSPeakList, negate)
{
    # UNDONE: this is DA specific and not tested for a while, maybe update some day...
    
    if (nrow(MSPeakList) == 0)
        return(MSPeakList)

    if (is.null(MSPeakList[["cmp"]]))
        stop("No isotope information available. Note that this is currently only implemented for DataAnalysis peak lists (if configured properly, see ?generateMSPeakLists.")

    # make sure most intense ions top the table
    MSPeakList <- MSPeakList[order(mz, -intensity)]

    unique_iso <- sapply(seq_along(MSPeakList$cmp), function(i)
    {
        # first and unassigned compounds are always unique
        if (i == 1 || !nzchar(MSPeakList$cmp[i]))
            return(TRUE)

        # peak may belong to multiple isotope compounds (separated by whitespace)
        cmp <- strsplit(MSPeakList$cmp[i], "\\s+")

        # isotope compound present before this entry?
        othercmp <- MSPeakList[seq_len(i - 1)][nzchar(cmp)]$cmp
        for (ocmp in othercmp)
        {
            if (any(cmp %in% strsplit(ocmp, "\\s+")))
                return(FALSE)
        }

        return(TRUE)
    }, USE.NAMES = FALSE)

    if (negate)
        unique_iso <- !unique_iso

    return(MSPeakList[unique_iso])
}

doMSPeakListFilter <- function(pList, absIntThr, relIntThr, topMost, minPeaks, maxMZOverPrec, minAbundanceFeatAbs,
                               minAbundanceFeatRel, minAbundanceFGroupAbs, minAbundanceFGroupRel, deIsotope, removeMZs,
                               retainPrecursor, precursorMZ, mzWindow, negate)
{
    ret <- pList

    intPred <- if (negate) function(i, t) i < t else function(i, t) i >= t

    if (!is.null(absIntThr))
        ret <- ret[intPred(intensity, absIntThr)]

    if (!is.null(relIntThr) && nrow(ret) > 0)
    {
        thr <- max(ret$intensity) * relIntThr
        ret <- ret[intPred(intensity, thr)]
    }

    if (!is.null(topMost) && nrow(ret) > topMost)
    {
        ord <- order(ret$intensity, decreasing = !negate)
        ret <- ret[seq_len(.N) %in% ord[seq_len(topMost)]] # NOTE: keep order
    }

    if (!is.null(maxMZOverPrec) && length(precursorMZ) != 0)
    {
        thr <- precursorMZ + maxMZOverPrec
        pred <- if (negate) function(x) x >= thr else function(x) x < thr
        ret <- ret[pred(mz)]
    }
    
    if (!is.null(minAbundanceFeatAbs))
    {
        if (!is.null(ret[["feat_abundance_abs"]]))
            ret <- ret[intPred(feat_abundance_abs, minAbundanceFeatAbs)] # averaged PL
        else
            ret <- ret[intPred(abundance_abs, minAbundanceFeatAbs)] # feat PL
    }
    if (!is.null(minAbundanceFeatRel))
    {
        if (!is.null(ret[["feat_abundance_rel"]]))
            ret <- ret[intPred(feat_abundance_rel, minAbundanceFeatRel)] # averaged PL
        else
            ret <- ret[intPred(abundance_rel, minAbundanceFeatRel)] # feat PL
    }
    if (!is.null(minAbundanceFGroupAbs) && !is.null(ret[["fgroup_abundance_abs"]]))
        ret <- ret[intPred(fgroup_abundance_abs, minAbundanceFGroupAbs)]
    if (!is.null(minAbundanceFGroupRel) && !is.null(ret[["fgroup_abundance_rel"]]))
        ret <- ret[intPred(fgroup_abundance_rel, minAbundanceFGroupRel)]

    if (deIsotope)
        ret <- deIsotopeMSPeakList(ret, negate)

    if (!is.null(removeMZs) && nrow(ret) > 0 && length(removeMZs) > 0)
    {
        pred <- function(mz) sapply(mz, function(x) !any(numLTE(abs(x - removeMZs), mzWindow)))
        if (negate)
            pred <- Negate(pred)
        ret <- ret[pred(mz)]
    }
    
    # re-add precursor if necessary
    if (retainPrecursor && !any(ret$precursor))
    {
        prec <- pList[precursor == TRUE]
        if (nrow(prec) > 0)
        {
            ret <- rbind(ret, prec)
            setorderv(ret, "mz")
        }
    }

    if (!is.null(minPeaks))
    {
        notEnough <- (nrow(ret)-1) < minPeaks # -1: don't count precursor peaks
        if (negate != notEnough)
            ret <- ret[0]
    }
    
    return(ret)
}

isolatePrecInMSPeakList <- function(plist, isolatePrec, negate)
{
    prec <- plist[precursor == TRUE]
    if (nrow(prec) == 1) # 0 if no precursor
    {
        s <- seq_len(isolatePrec$maxIsotopes) / isolatePrec$z
        mzranges <- matrix(c(s + isolatePrec$mzDefectRange[1] * s,
                             s + isolatePrec$mzDefectRange[2] * s), ncol = 2) + prec$mz
        keep <- plist[, precursor | (
            # only keep peaks with reasonable intensities
            intensity %between% (isolatePrec$intRange * prec$intensity) &

            # only keep with reasonably close m/z to precursor, taking
            # in to account larger windows for larger distances
            inrange(mz, mzranges[, 1], mzranges[, 2]))]
        plist <- plist[if (negate) !keep else keep]

        # remove gaps
        gaps <- plist[round(mz - shift(mz)) > (isolatePrec$maxGap / isolatePrec$z), which = TRUE]
        if (length(gaps) > 0)
        {
            sq <- seq_len(gaps[1] - 1)
            plist <- plist[if (negate) -sq else sq]
        }
    }

    return(plist)
}

getSpec <- function(MSPeakLists, groupName, MSLevel, analysis)
{
    MSInd <- if (MSLevel == 1) "MS" else "MSMS"
    if (!is.null(analysis))
        spec <- MSPeakLists[[analysis, groupName]][[MSInd]]
    else
        spec <- MSPeakLists[[groupName]][[MSInd]]
    return(spec)
}

getMSPeakListPlotTitle <- function(MSLevel, analysis, groupName)
{
    MSInd <- if (MSLevel == 1) "MS" else "MSMS"
    if (!is.null(analysis))
        return(sprintf("%s (%s) %s", groupName, analysis, MSInd))
    return(paste(groupName, MSInd))
}

# NOTE: This function is mainly for debugging the C++ version
# nocov start
binPeakLists <- function(pl1, pl2, shift, absMzDev)
{
    prep <- function(pl)
    {
        pl <- copy(pl)
        pl[, c("index", "low", "high") := .(seq_len(nrow(pl)), mz - absMzDev, mz + absMzDev)]
        setkeyv(pl, c("low", "high"))
        return(pl)
    }
    pl1 <- prep(pl1); pl2 <- prep(pl2)
    
    if (shift != "none")
    {
        if (!any(pl1$precursor) || !any(pl2$precursor))
            stop("Cannot shift spectra: one or both lack precursor ion!")
        precDiff <- pl2[precursor == TRUE]$mz - pl1[precursor == TRUE]$mz
        if (shift == "precursor")
            pl2[, mz := mz - precDiff]
        else # both
        {
            # first bin as normal (recursive call)
            binNone <- binPeakLists(pl1, pl2, "none", absMzDev)
            pl1Unique <- setnames(binNone[intensity_2 == 0, -"intensity_2"], "intensity_1", "intensity")
            pl2Unique <- setnames(binNone[intensity_1 == 0, -"intensity_1"], "intensity_2", "intensity")
            
            # bin missing with shift
            pl2Unique[, mz := mz - precDiff]
            binShift <- binPeakLists(pl1Unique, pl2Unique, "none", absMzDev)
            
            # merge both: add missing from binNone
            ret <- rbind(binNone[intensity_1 != 0 & intensity_2 != 0], binShift)
            setorderv(ret, "mz")
            return(ret)
        }
    }
    
    ov <- foverlaps(pl1, pl2)
    # NOTE: i.xx cols are from left
    
    ret <- ov[, c("i.mz", "intensity", "i.intensity"), with = FALSE]
    setnames(ret, c("i.mz", "intensity", "i.intensity"), c("mz", "intensity_2", "intensity_1"))

    # all left entries should be in there, right can miss and should be added manually
    stopifnot(all(pl1$index %in% ov$i.index))
    ret <- rbind(ret, setnames(pl2[!index %in% ov$index, c("mz", "intensity")], "intensity", "intensity_2"), fill = TRUE)
    
    setnafill(ret, fill = 0)
    setorderv(ret, "mz")
    setcolorder(ret, c("mz", "intensity_1", "intensity_2"))
    
    return(ret)
}

# NOTE: This function is mainly for debugging the C++ version
specSimilarityR <- function(pl1, pl2, method, shift = "none", removePrecursor = FALSE, mzWeight = 0, intWeight = 1,
                           absMzDev = defaultLim("mz", "medium"), relMinIntensity = 0.1)
{
    # UNDONE: refs
    
    # code contributed by Bas van de Velde
    # cosine similarity from OrgMassSpecR
    
    binnedPL <- as.data.table(binSpecCPP(pl1, pl2, shift, absMzDev))
    # binnedPL <- binPeakLists(pl1, pl2, shift, absMzDev)
    
    # remove precursor
    if (removePrecursor)
        binnedPL <- binnedPL[precursor == FALSE]
    
    # normalize
    binnedPL[, c("intensity_1", "intensity_2") :=
                 .(normalize(intensity_1, FALSE), normalize(intensity_2, FALSE))]
    
    # intensity filter
    binnedPL[intensity_1 < relMinIntensity, intensity_1 := 0]
    binnedPL[intensity_2 < relMinIntensity, intensity_2 := 0]
    binnedPL <- binnedPL[intensity_1 != 0 | intensity_2 != 0]
    
    u <- binnedPL$mz^mzWeight * binnedPL$intensity_1^intWeight
    v <- binnedPL$mz^mzWeight * binnedPL$intensity_2^intWeight

    return(switch(method,
                  cosine = as.vector((u %*% v) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))),
                  pearson = cor(u, v, method = "pearson"),
                  spearman = cor(u, v, method = "spearman"),
                  jaccard = binnedPL[intensity_1 != 0 & intensity_2 != 0, .N] / nrow(binnedPL)))
}
# nocov end

getSimPLAndPrec <- function(MSPeakLists, group, analysis, MSLevel, specSimParams, nr)
{
    if (is.null(analysis))
        analysis <- rep(list(NULL), length(group))
    else if (length(analysis) != length(group))
        stop(sprintf("Length of analysis%d must equal the length of groupName%d", nr, nr))
    names(analysis) <- group
    
    specs <- pruneList(Map(getSpec, group, analysis, MoreArgs = list(MSPeakLists = MSPeakLists, MSLevel = MSLevel)))
    if (length(specs) == 0)
        return(NULL)
    
    analysis <- analysis[names(specs)]
    
    precs <- NULL
    if (specSimParams$shift != "none")
    {
        precs <- mapply(specs, names(specs), analysis, FUN = function(spec, gn, ana)
        {
            ret <- spec[precursor == TRUE]$mz
            if (length(ret) == 0)
            {
                if (MSLevel == 2)
                {
                    precSpec <- getSpec(MSPeakLists, gn, 1, ana)
                    ret <- if (is.null(precSpec))
                        numeric()
                    else
                        precSpec[precursor == TRUE]$mz
                }
            }
            return(if (length(ret) == 0) NA_real_ else ret)            
        })
        precsNA <- names(which(is.na(precs)))
        if (length(precsNA) != 0)
            warning("Some pecursor ions are unknown, returning NAs")
        wh <- setdiff(names(precs), precsNA)
        precs <- precs[wh]; specs <- specs[wh]
    }
    else
        precs <- setNames(rep(0, length(specs)), names(specs))
    
    specs <- pruneList(lapply(specs, prepSpecSimilarityPL, removePrecursor = specSimParams$removePrecursor,
                              relMinIntensity = specSimParams$relMinIntensity, minPeaks = specSimParams$minPeaks),
                       checkZeroRows = TRUE)
    
    return(list(specs = specs, precs = precs))
}

prepSpecSimilarityPL <- function(pl, removePrecursor, relMinIntensity, minPeaks)
{
    pl <- copy(pl)
    
    if (removePrecursor)
        pl <- pl[precursor == FALSE]
    
    if (relMinIntensity > 0 && nrow(pl) > 0)
    {
        minInt <- relMinIntensity * max(pl$intensity)
        pl <- pl[numGTE(intensity, minInt)]
    }
    
    if (nrow(pl) < minPeaks)
        return(pl[FALSE])
    
    return(pl)
}

getBinnedPLPair <- function(MSPeakLists, groupNames, analyses, MSLevel, specSimParams, uniqueName, mustExist,
                            normalizedIntensities = TRUE)
{
    PLP1 <- getSimPLAndPrec(MSPeakLists, groupNames[1], analyses[1], MSLevel, specSimParams, 1)
    PLP2 <- getSimPLAndPrec(MSPeakLists, groupNames[2], analyses[2], MSLevel, specSimParams, 2)
    
    if (is.null(PLP1) || is.null(PLP2))
    {
        if (!mustExist)
        {
            # create dummy output
            dummySpec <- function(PLP)
            {
                if (is.null(PLP))
                {
                    emptyPL <- emptyMSPeakList(NULL, NULL)
                    emptyPL[, c("ID", "mergedBy") := character()]
                    setcolorder(emptyPL, "ID")
                    return(emptyPL)
                }
                sp <- PLP$specs[[1]]
                sp[, mergedBy := uniqueName]
                sp[, intensity := intensity / max(intensity)]
                sp <- removeDTColumnsIfPresent(sp, c("fgroup_abundance_rel", "fgroup_abundance_abs", "feat_abundance_rel",
                                                     "feat_abundance_abs"))
                return(sp)
            }
            return(list(spec1 = dummySpec(PLP1), spec2 = dummySpec(PLP2)))
        }
        
        if (is.null(PLP1))
            stop("Could not obtain first spectrum")
        if (is.null(PLP2))
            stop("Could not obtain second spectrum")
    }
    
    precDiff <- 0
    if (specSimParams$shift != "none")
    {
        if (is.na(PLP1$precs) || is.na(PLP2$precs))
            stop("One or both precursor ions are unknown, can't calculate shift")
        precDiff <- PLP2$precs - PLP1$precs
    }
    
    bin <- as.data.table(binSpectra(PLP1$specs[[1]], PLP2$specs[[1]], specSimParams$shift, precDiff, specSimParams$absMzDev))
    bin[, mergedBy := fifelse(intensity_1 != 0 & intensity_2 != 0, "overlap", uniqueName)]
    
    getSpecFromBin <- function(nr)
    {
        othernr <- if (nr == 1) 2 else 1
        # remove other columns and rename intensity
        rmCols <- paste0(c("intensity_", "ID_"), othernr)
        ret <- setnames(bin[get(paste0("intensity_", nr)) != 0, setdiff(names(bin), rmCols), with = FALSE],
                        paste0(c("intensity_", "ID_"), nr), c("intensity", "ID"))
        setorderv(ret, "ID") # restore order
        
        if (!normalizedIntensities)
        {
            PLP <- if (nr == 1) PLP1 else PLP2
            ret[, intensity := PLP$specs[[1]]$intensity[match(ID, PLP$specs[[1]]$ID)]]
        }
        
        # re-add precursor
        ret[, precursor := if (nr == 1) PLP1$specs[[1]][match(ret$ID, ID)]$precursor else PLP2$specs[[1]][match(ret$ID, ID)]$precursor]
        
        setcolorder(ret, c("ID", "mz", "intensity", "precursor", "mergedBy"))
        
        return(ret)
    }
    
    return(list(spec1 = getSpecFromBin(1), spec2 = getSpecFromBin(2)))
}

expandFillSpecSimilarities <- function(sims, groupName1, groupName2)
{
    if (length(sims) == 0)
    {
        sims <- matrix(NA_real_, nrow = length(groupName1), ncol = length(groupName2))
        rownames(sims) <- groupName1; colnames(sims) <- groupName2
    }
    else
    {
        nMissing <- length(groupName1) - nrow(sims)
        
        if (nMissing > 0)
        {
            # add missing rows
            sims <- do.call(rbind, c(list(sims), setNames(as.list(rep(NA_real_, nMissing)),
                                                          setdiff(groupName1, rownames(sims)))))
            sims <- sims[groupName1, , drop = FALSE]
        }
        
        nMissing <- length(groupName2) - ncol(sims)
        if (nMissing > 0)
        {
            # add missing columns
            sims <- do.call(cbind, c(list(sims), setNames(as.list(rep(NA_real_, nMissing)),
                                                          setdiff(groupName2, colnames(sims)))))
            sims <- sims[, groupName2, drop = FALSE]
        }
    }
    
    return(sims)
}

mergeBinnedAndAnnPL <- function(binPL, annPL, which)
{
    # get rid of duplicate columns (except ID)
    annPL <- annPL[, c("ID", setdiff(names(annPL), names(binPL))), with = FALSE]
    
    annPL <- merge(binPL, annPL, by = "ID")
    
    annPL[, which := which]
    return(annPL)
}
