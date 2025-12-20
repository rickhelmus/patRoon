# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include features.R
NULL

getPiekEICsInfo <- function(params, IMS, suspects, MS2Info, verbose)
{
    # UNDONE: also support other binning approaches?
    binsMZ <- seq(params$mzRange[1], params$mzRange[2], by = params$mzStep * 0.5)
    names(binsMZ) <- paste0("bin_M", binsMZ)
    
    EICInfo <- data.table(mzmin = binsMZ, mzmax = binsMZ + params$mzStep, EIC_ID_MZ = names(binsMZ))
    
    if (nrow(EICInfo) == 0)
        stop("Cannot form EICs with the given m/z parameters (there are no bins)! Please check the m/z range and step size.", call. = FALSE)
    
    if (IMS)
    {
        binsIMS <- seq(params$mobRange[1], params$mobRange[2], by = params$mobStep * 0.5)
        names(binsIMS) <- paste0("bin_I", binsIMS)
        tab <- CJ(EIC_ID_MZ = names(binsMZ), EIC_ID_IMS = names(binsIMS), sorted = FALSE)
        tab[, c("mobmin", "mobmax") := .(binsIMS[EIC_ID_IMS], binsIMS[EIC_ID_IMS] + params$mobStep)]
        EICInfo <- merge(EICInfo, tab, by = "EIC_ID_MZ", sort = FALSE)
    }
    else
        EICInfo[, c("mobmin", "mobmax") := 0]
    
    if (params$filter == "suspects")
    {
        if (IMS && params$filterIMS == "suspects")
            EICInfo[, keep := filterEICBins(mzmin, params$mzStep, mobmin, params$mobStep,
                                            suspects$mz - params$mzWindow, suspects$mz + params$mzWindow,
                                            suspects$mobility - params$IMSWindow, suspects$mobility + params$IMSWindow)]
        else
            EICInfo[, keep := filterEICBins(mzmin, params$mzStep, numeric(), 0,
                                            suspects$mz - params$mzWindow, suspects$mz + params$mzWindow, 0, 0)]
        if (verbose)
        {
            printf("Removed %d (%.2f%%) EICs after suspect filtering. Remaining: %d\n", sum(!EICInfo$keep),
                   sum(!EICInfo$keep) / nrow(EICInfo) * 100, sum(EICInfo$keep))
        }
        EICInfo <- EICInfo[keep == TRUE][, keep := NULL]
    }
    else if (params$filter == "ms2")
    {
        if (IMS && params$filterIMS ==  "ms2")
            EICInfo[, keep := filterEICBins(mzmin, params$mzStep, mobmin, params$mobStep, MS2Info$mzmin, MS2Info$mzmax,
                                            MS2Info$mobmin, MS2Info$mobmax)]
        else
            EICInfo[, keep := filterEICBins(mzmin, params$mzStep, numeric(), 0, MS2Info$mzmin, MS2Info$mzmax, 0, 0)]
        if (verbose)
        {
            printf("Removed %d (%.2f%%) EICs after MS2 filtering. Remaining: %d\n", sum(!EICInfo$keep),
                   sum(!EICInfo$keep) / nrow(EICInfo) * 100, sum(EICInfo$keep))
        }
        EICInfo <- EICInfo[keep == TRUE][, keep := NULL]
    }
    
    EICInfo[, EIC_ID := paste0("EIC_", .I)]
    
    return(EICInfo)
}

#' @rdname features-class
#' @export
featuresPiek <- setClass("featuresPiek", slots = c(mzProfiles = "list", EIMs = "list"), contains = "features")

setMethod("initialize", "featuresPiek",
          function(.Object, ...) callNextMethod(.Object, algorithm = "piek", ...))

#' @rdname features-class
#' @export
setMethod("delete", "featuresPiek", function(obj, i = NULL, j = NULL, ...)
{
    old <- obj
    obj <- callNextMethod()

    if (length(obj@mzProfiles) > 0)
    {
        obj@mzProfiles <- obj@mzProfiles[intersect(names(obj@EIMs), analyses(obj))]
        obj@mzProfiles <- Map(obj@mzProfiles, featureTable(obj), f = \(eims, ft) eims[intersect(names(eims), ft$ID)])
    }
    
    if (length(obj@EIMs) > 0)
    {
        obj@EIMs <- obj@EIMs[intersect(names(obj@EIMs), analyses(obj))]
        obj@EIMs <- Map(obj@EIMs, featureTable(obj), f = \(eims, ft) eims[intersect(names(eims), ft$ID)])
    }
    
    return(obj)
})


#' Find features using piek
#'
#' Uses the \code{piek} algorithm to find features.
#'
#' @templateVar algo piek
#' @templateVar do automatically find features
#' @templateVar generic findFeatures
#' @templateVar algoParam piek
#' @template algo_generator
#'
#' @details The \code{piek} algorithm extends and improves on the simple and fast feature detection algorithm introduced
#'   by \insertRef{Dietrich2021}{patRoon}. This algorithm first forms extracted ion chromatograms (EICs) and
#'   subsequently performs automatic peak detection to generate features. The piek algorithm introduces the following
#'   improvements and changes: \itemize{
#'
#'     \item Support for IMS workflows.
#'
#'     \item The \link{msdata} interface is used to efficiently form EICs from the raw data. All the file formats and
#'     types that are supported by \link{msdata}. This includes IMS data, even if not used for feature detection, which
#'     allows the use of IMS data directly in non-IMS or \link[=assignMobilities_feat]{post mobility assignment}
#'     workflows.
#'
#'     \item The EICs may be formed by different approaches: \enumerate{
#'
#'       \item the use of pre-defined m/z bins with fixed widths. This is the approach of the original algorithm by
#'       \insertRef{Dietrich2021}{patRoon} and applicable for most workflows.
#'
#'       \item the extension of m/z bins with ion mobility data to form two-dimensional bins (m/z + mobility).
#'
#'       \item generate EICs from a suspect list (m/z or m/z + mobility)
#'
#'       \item generate EICs from the precursor ions that were found in an data-dependent (DDA) MS/MS experiment (m/z)
#'      or Bruker DDA-PASEF experiment (m/z + mobility)
#'
#'     }
#'
#'     \item The original peak detection algorithm was further optimized or can be be exchanged with others: see
#'     \code{\link{getDefPeakParams}} for details. If retention times are available, \emph{i.e.} when EICs are formed
#'     from the third approach (and retention times are available in the suspect list) or fourth approach, then these
#'     will be used to filter the detected peaks.
#'
#'     \item Several filters are applied to improve the data and reduce redundancy. These are discussed in the next
#'     sections.
#'
#'   }
#'
#'   The inclusion of ion mobility data in the EIC formation initiates a \link[=assignMobilities_feat]{direct mobility
#'   assignment} workflow. Combinations of these approaches are also possible in IMS workflows
#'   (see the \verb{EIC formation parameters} section below).
#'
#'   The m/z and mobility values (in IMS workflows) assigned to the feature are derived from the weighted mean of the
#'   base peaks from the EIC data in the retention time range of the feature.
#'
#' @param genEICParams A \code{list} of parameters for the EIC generation. See the \verb{EIC formation parameters}
#'   section below. The \code{getPiekEICParams} is used to generate the parameter list.
#' @param peakParams A \code{list} of parameters for the peak detection. See \code{\link{getDefPeakParams}} for details.
#' @param suspects A suspect list that needs to be specified if EICs are formed from suspect data. See
#'   \link[=suspect-screening]{suspect screening} for details on the suspect list format.
#' @param adduct An \code{\link{adduct}} object (or something that can be converted to it with \code{\link{as.adduct}}).
#'   Examples: \code{"[M-H]-"}, \code{"[M+Na]+"}. Only needs to be specified when EICs are formed from suspect data.
#' @param assignMethod Should be \code{"basepeak"} or \code{"weighted.mean"}. This parameter sets how measured
#'   \emph{m/z} or mobilities across the EIC datapoints are handled for feature assignment. If
#'   \code{assignMethod="basepeak"}, then the value of the base peak (=highest intensity peak) from each EIC datapoint
#'   is taken. If \code{assignMethod="weighted.mean"} then the intensity weighted mean is calculated of the values that
#'   fall within the EIC bin.
#' @param assignAggr Should be \code{"max"} or \code{"weighted.mean"}. This parameter sets how the \emph{m/z} or
#'   mobility values determined for each datapoint (see \code{assignMethod}) are aggregated to determine the final
#'   feature value. If \code{assignAggr="weighted.mean"} then the intensity weighted mean is used. With
#'   \code{assignAggr="max"} the data at the the highest intensity of the chromatographic peak is taken.
#'
#' @inheritParams findFeatures
#' @template minIntensityIMS-arg
#'
#' @section EIC formation parameters: The \code{genEICParams} argument to \code{findFeaturesPiek} configures the
#'   formation of EICs. The \code{getPiekEICParams} function should be used to generate the parameter list.
#'
#'   The following general parameters exist: \itemize{
#'
#'     \item \code{filter} Sets how m/z data is used for the formation of EICs. Possible values are \code{"bins"} (use
#'     equally sized m/z bins), \code{"suspects"} (use the unique m/z values from a suspect list) and \code{"ms2"} (use
#'     the m/z values from the precursors detected in a data-dependent MS/MS experiment).
#'
#'     \item \code{filterIMS} Equivalent as \code{filter}, but for ion mobility data. If \code{NULL}, no IMS data will
#'     be used and no \link[=assignMobilities_feat]{direct mobility assignment} IMS workflow is initiated. Different
#'     values for \code{filter} and \code{filterIMS} can be specified to combine approaches. The following
#'     combinations are supported: \itemize{
#'
#'       \item \code{filter="bins"} and \code{filterIMS="bins"}
#'
#'       \item \code{filter="suspects"} and \code{filterIMS="suspects"}
#'
#'       \item \code{filter="suspects"} and \code{filterIMS="bins"}
#'
#'       \item \code{filter="ms2"} and \code{filterIMS="ms2"}
#'
#'       \item \code{filter="ms2"} and \code{filterIMS="bins"}
#'
#'     }
#'
#'     Currently only Bruker DDA-PASEF experiments provide the data needed for the \code{"ms2"} approach.
#'
#'     \item \code{retRange} A \code{numeric} vector of length two that specifies the retention time range for the EICs.
#'     Data outside this range is excluded. Set to \code{NULL} to use the full range.
#'
#'     \item \code{gapFactor} A \code{numeric} that configures gap filling for EICs. See \code{\link{getDefEICParams}}
#'     for further details.
#'
#'     \item \code{minEICIntensity} The minimum intensity of the highest data point in the EIC. Used to filter EICs.
#'
#'     \item \code{minEICAdjTime},\code{minEICAdjPoints},\code{minEICAdjIntensity} The EIC should have at least a
#'     continuous signal of \code{minEICAdjTime} seconds and \code{minEICAdjPoints} data points, where the continuity is
#'     defined by data points with an intensity of at least \code{minEICAdjIntensity} high. Set \code{minEICAdjTime} or
#'     \code{minEICAdjPoints} to zero to disable continuity checks for time or data points, respectively. Set
#'     \code{minEICAdjIntensity} to zero to completely disable continuity checks.
#'
#'     \item \code{topMostEIC} Only keep this number of top-most intense EICs. The intensity is derived from the data
#'     point with the highest intensity in the EIC. Set to zero to always select all EICs.
#'
#'     \item \code{topMostEICPre} Equivalent to \code{topMostEIC}, but used for pre-checking EICs in IMS workflows
#'     (discussed in the next section).
#'
#'   }
#'
#'   The following parameters are specific for EIC binning: \itemize{
#'
#'     \item \code{mzRange},\code{mzStep} Configures the formation of m/z bins. \code{mzRange} is a numeric vector of
#'     length two that specifies the min/max m/z range. \code{mzStep} specifies the bin widths.
#'
#'     \item \code{mobRange},\code{mobStep} Equivalent to above, but for ion mobility binning (\emph{i.e.} if
#'     \code{filterIMS="bins"}).
#'
#'   }
#'
#'   The following parameters are specifically for EICs from suspect data: \itemize{
#'
#'     \item \code{rtWindow},\code{mzwindow},\code{IMSWindow}: If retention times are present in the suspect list:
#'     specify the retention time, m/z and mobility (if \code{filterIMS="suspects"}) tolerance window to match features
#'     with suspects. This is done when eliminating features with deviating retention times. Set \code{rtWindow=Inf} to
#'     disable this step.
#'
#'     The \code{mzWindow} and \code{IMSWindow} parameters are also used to set the data windows for the EICs.
#'
#'     Defaults to \code{defaultLim("retention", "medium")}, \code{defaultLim("mz", "medium")} and
#'     \code{defaultLim("mobility", "medium")}, see \link{limits}.
#'
#'     \item \code{skipInvalid},\code{prefCalcChemProps},\code{neutralChemProps} Controls preparing the suspect list
#'     data. See \code{\link{screenSuspects}}.
#'
#'
#'   }
#'
#'   The following parameters are specifically for EICs from MS/MS data: \itemize{
#'
#'     \item \code{rtWindow} Eliminates any features without an MS/MS spectrum within this retention time window. Set
#'     \code{rtWindow=Inf} to disable this filter. Defaults to \code{defaultLim("retention", "very_narrow")} (see
#'     \link{limits}).
#'
#'     \item \code{mzWindow},\code{IMSWindow} The m/z and mobility (if \code{filterIMS="ms2"}) tolerance windows, used to \enumerate{
#'
#'       \item match features to MS/MS spectrum retention times (see \code{rtWindow}).
#'
#'       \item set the data windows for the EICs.
#'
#'       \item used as clustering width to average the m/z and mobility data of MS/MS precursor (see the
#'       \verb{Elimination and averaging of redundant data} section).
#'
#'     }
#'
#'     Defaults to \code{defaultLim("mz", "narrow")} and \code{defaultLim("mobility", "medium")} (see \link{limits}).
#'
#'     \item \code{minTIC} The minimum total ion current (TIC) signal for an MS/MS spectrum to be considered. Can be
#'     increased to eliminate features with low intensity MS/MS data.
#'
#'     \item \code{clusterMethod} The clustering method to average MS/MS precursor data (discussed further in the
#'     \verb{Elimination and averaging of redundant data} section).
#'
#'   }
#'
#' @section Elimination and averaging of redundant data: If EICs are formed from suspect data, then any duplicates in
#'   the suspect list are eliminated to avoid duplicate EIC formation. Duplicate suspects are defined by very close m/z
#'   and mobility (in IMS workflows), with tolerance windows defined by \code{defaultLim("mz", "very_narrow")} and
#'   \code{defaultLim("mobility", "very_narrow")}, respectively (see \link{limits}).
#'
#'   If EICs are formed from MS/MS data, then the m/z and mobility data (in IMS workflows) of the MS/MS precursors are
#'   clustered and averaged to avoid duplicate EIC formation. The clustering width and method are defined by the
#'   \code{mzWindow}, \code{IMSWindow} and \code{clusterMethod} parameters, respectively (see previous section and
#'   \link[=cluster-params]{clustering parameters}).
#'
#'   In IMS workflows a 'pre-check' is formed to reduce the number of EICs that should be subjected to peak detection.
#'   This is especially important when both \code{filter} and \code{filterIMS} are set to \code{"bins"}, as this may
#'   lead to millions of EIC bins in total. During the pre-check EICs are formed by only considering the m/z dimension,
#'   and these are subsequently filtered by parameters described in the previous section. The final EICs for feature
#'   detection are then only formed if they have m/z data that was not removed during the pre-check.
#'
#'   After feature detection, duplicate features are detected by very close retention time, m/z, and mobility (in IMS
#'   workflows), and only the feature with highest intensity is kept. The equivalence is determined from the tolerance
#'   windows defined by \code{defaultLim("retention", "very_narrow")} \code{defaultLim("mz", "very_narrow")} and
#'   \code{defaultLim("mobility", "very_narrow")}, respectively (see \link{limits}).
#'
#'
#' @templateVar what \code{findFeaturesPiek}
#' @template uses-msdata
#'
#' @references \insertAllCited{}
#'
#' @inherit findFeatures return
#'
#' @export
findFeaturesPiek <- function(analysisInfo, genEICParams, peakParams, suspects = NULL, adduct = NULL, IMS = FALSE,
                             assignMethod = "basepeak", assignAggr = "weighted.mean", minIntensityIMS = 25,
                             assignRTWindow = defaultLim("retention", "very_narrow"), EICBatchSize = Inf, verbose = TRUE)
{
    # UNDONE: add refs to docs, and highlight changes
    # UNDONE: test empties, eg no EICs
    
    if (!is.null(adduct))
        adduct <- checkAndToAdduct(adduct)
    
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    assertPiekGenEICParams(genEICParams, add = ac)
    assertFindPeakParams(peakParams, add = ac)
    if (genEICParams$filter == "suspects")
        assertSuspectList(suspects, needsAdduct = is.null(adduct), skipInvalid = genEICParams$skipInvalid, null.ok = FALSE,
                          add = ac)
    checkmate::assertFlag(IMS, add = ac)
    checkmate::assertChoice(assignMethod, c("basepeak", "weighted.mean"), add = ac)
    checkmate::assertChoice(assignAggr, c("max", "weighted.mean"), add = ac)
    aapply(checkmate::assertNumber, . ~ minIntensityIMS + assignRTWindow, lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::assertNumber(EICBatchSize, lower = 1, finite = FALSE, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    maybePrintf <- \(...) if (verbose) printf(...)
    
    if (genEICParams$filter == "suspects")
    {
        if (nrow(suspects) == 0)
            stop("The suspect list empty.", call. = FALSE)
        
        suspects <- prepareSuspectList(suspects, adduct = adduct, skipInvalid = genEICParams$skipInvalid,
                                       checkDesc = TRUE, prefCalcChemProps = genEICParams$prefCalcChemProps,
                                       neutralChemProps = genEICParams$neutralChemProps)
        
        whNotInRange <- suspects[!mz %between% genEICParams$mzRange, which = TRUE]
        if (length(whNotInRange) > 0)
            warning(sprintf("The following %d suspect rows have m/z values outside the binning range and will be ignored: %s",
                            length(whNotInRange), paste0(whNotInRange, collapse = ", ")),
                    call. = FALSE)
        
        if (genEICParams$filterIMS == "suspects")
        {
            suspects[, mobility_susp := selectFromSuspAdductCol(suspects, "mobility", if (!is.null(adduct)) as.character(adduct))]
            whMissing <- which(is.na(suspects$mobility_susp))
            if (length(whMissing) > 0)
            {
                if (length(whMissing) == nrow(suspects))
                    stop("The suspect list contains no mobility data! ",
                         "Please note that CCS data first must be converted to mobilities, e.g. with assignMobilities().",
                         call. = FALSE)
                
                warning(sprintf("The suspect list contains missing mobility data (rows %s), these will be ignored. ",
                                paste0(whMissing, collapse = ", ")),
                        "Please note that CCS data first must be converted to mobilities, e.g. with assignMobilities().",
                        call. = FALSE)
                suspects <- suspects[-whMissing]
            }
            suspectsOrig <- suspects # to report original row numbers below
            suspects <- expandSuspMobilities(suspects)
            
            whNotInRange <- suspects[!mobility %between% genEICParams$mobRange, which = TRUE]
            if (length(whNotInRange) > 0)
                warning(sprintf("The following %d suspect rows have mobilities values outside the binning range and will be ignored: %s",
                                length(whNotInRange), paste0(match(suspects$name[whNotInRange], suspectsOrig$name), collapse = ", ")),
                        call. = FALSE)
            
        }
    }
    
    if (is.null(genEICParams[["retRange"]]))
        genEICParams$retRange <- c(0, 0)
    
    cacheDB <- openCacheDBScope()
    baseHash <- makeHash(genEICParams, peakParams, suspects, adduct, assignMethod, assignAggr, minIntensityIMS,
                         assignRTWindow)
    anaHashes <- getMSFileHashesFromAvailBackend(analysisInfo, needIMS = IMS)
    anaHashes <- sapply(anaHashes, makeHash, baseHash)
    cachedData <- pruneList(loadCacheData("featuresPiek", anaHashes, simplify = FALSE, dbArg = cacheDB))
    if (length(cachedData) > 0)
    {
        names(cachedData) <- names(anaHashes)[match(names(cachedData), anaHashes)]
        anaInfoTBD <- analysisInfo[!analysis %in% names(cachedData)]
    }
    else
        anaInfoTBD <- analysisInfo
    
    getEICsAna <- function(backend, EICInfo, mode, topMost)
    {
        ret <- getEICList(backend, EICInfo$mzmin, EICInfo$mzmax, genEICParams$retRange[1], genEICParams$retRange[2],
                          EICInfo$mobmin, EICInfo$mobmax, gapFactor = genEICParams$gapFactor, mzExpIMSWindow = 0,
                          minIntensityIMS = minIntensityIMS, mode = mode, sumWindowMZ = genEICParams$sumWindowMZ,
                          sumWindowMob = genEICParams$sumWindowMob, smoothWindowMZ = genEICParams$smoothWindowMZ,
                          smoothExtMZ = genEICParams$smoothExtMZ, smoothWindowMob = genEICParams$smoothWindowMob,
                          smoothExtMob = genEICParams$smoothExtMob, saveMZProfiles = genEICParams$saveMZProfiles,
                          saveEIMs = genEICParams$saveEIMs, pad = FALSE, minEICIntensity = genEICParams$minEICIntensity,
                          minEICAdjTime = genEICParams$minEICAdjTime, minEICAdjPoints = genEICParams$minEICAdjPoints,
                          minEICAdjIntensity = genEICParams$minEICAdjIntensity, topMost = topMost)
        names(ret) <- EICInfo$EIC_ID
        if (mode != "test")
            ret <- pruneList(ret, checkZeroRows = TRUE, keepAttr = TRUE)
        return(ret)
    }
    
    assignMZOrMobsToPeaks <- function(peaks, EICInfo, what)
    {
        chkCol <- paste0(what, if (assignAggr == "max") "BPMax" else "BP")
        minCol <- if (what == "mz") "mzmin" else "mobmin"
        centCol <- paste0(what, "Centered")
        step <- if (what == "mz") genEICParams$mzStep else genEICParams$mobStep
        peaks[, binStart := EICInfo[match(peaks$EIC_ID, EIC_ID)][[minCol]]]
        peaks[, (centCol) := between(get(chkCol), binStart + step/4, binStart + step/4*3)]
        if (assignAggr == "max")
        {
            if (assignMethod == "basepeak")
                peaks[, (what) := get(paste0(what, "BPMax"))]
            else
                peaks[, (what) := get(paste0(what, "Max"))]
        }
        else # assignAggr == "weighted.mean"
        {
            if (assignMethod == "basepeak")
                peaks[, (what) := get(paste0(what, "BP"))]
            # else already weighted mean
        }
        # NOTE: centered is kept and removed later to allow filtering redundant peaks
        peaks[, c(paste0(what, c("BP", "BPMax", "Max")), "binStart") := NULL]
        return(peaks)
    }
    
    getIMSProfiles <- function(profAttr, peaks, EICs, EICInfo)
    {
        ret <- list()
        profAttr <- pruneList(profAttr) # NOTE: EICs/EICInfo is already pruned
        if (length(profAttr) > 0)
        {
            ret <- Map(peaks$EIC_ID, peaks$ret, f = function(EIC_ID, ret)
            {
                profAttr[[match(EIC_ID, EICInfo$EIC_ID)]][[which.min(abs(EICs[[EIC_ID]][, "time"] - ret))]]
            })
            names(ret) <- peaks$ID
        }
        return(ret)
    }
    
    fList <- list()
    mzProfiles <- list()
    EIMs <- list()
    if (nrow(anaInfoTBD) > 0)
    {
        fList <- applyMSData(anaInfoTBD, needIMS = IMS, showProgress = FALSE, func = function(ana, path, backend)
        {
            maybePrintf("\n-------\nFinding features for '%s'...\n", ana)
            
            openMSReadBackend(backend, path)
         
            MS2Info <- NULL
            if (genEICParams$filter == "ms2")
            {
                MS2Info <- if (IMS)
                    setDT(getIMSIsolationInfo(backend))
                else
                    setDT(getMSMetadata(backend, 2))
                setnames(MS2Info, c("mobStart", "mobEnd"), c("mobmin", "mobmax"), skip_absent = TRUE)
                if (!is.finite(genEICParams$mzIsoWindow))
                {
                    MS2Info[, c("mzmin", "mzmax") := .(precursorMZ - isolationRangeMin,
                                                       precursorMZ + isolationRangeMax)]
                }
                else
                {
                    MS2Info[, c("mzmin", "mzmax") := .(precursorMZ - min(genEICParams$mzIsoWindow, isolationRangeMin),
                                                       precursorMZ + min(genEICParams$mzIsoWindow, isolationRangeMax))]
                }
                MS2Info <- MS2Info[TIC >= genEICParams$minTIC]
            }
            
            EICs <- EICInfo <- NULL
            EICInfoMZ <- getPiekEICsInfo(genEICParams, FALSE, suspects, MS2Info, verbose)
            EICInfo <- if (IMS)
            {
                EICInfoMob <- getPiekEICsInfo(genEICParams, TRUE, suspects, MS2Info, verbose)
                
                if (nrow(EICInfoMob) > genEICParams$minEICsIMSPreCheck)
                {
                    maybePrintf("Pre-checking %d m/z EICs... ", nrow(EICInfoMZ))
                    testEICs <- getEICsAna(backend, EICInfoMZ, "test", genEICParams$topMostEICMZ)
                    testEICs <- unlist(testEICs)
                    EICInfoMZ <- EICInfoMZ[EIC_ID %chin% names(testEICs)[testEICs]]
                    maybePrintf("Done! Eliminated %d (%.2f%%) EICs\n", sum(!testEICs),
                                (sum(!testEICs) / length(testEICs)) * 100)
                    
                    # remove complete m/z bins that were filtered out before
                    temp <- EICInfoMZ[, c("mzmin", "mzmax"), with = FALSE]
                    setkeyv(temp, c("mzmin", "mzmax"))
                    ov <- foverlaps(EICInfoMob, temp, type = "within", nomatch = NULL, which = TRUE)
                    EICInfoMob[ov$xid]
                }
                else
                    EICInfoMob
            }
            else
                EICInfoMZ
            
            EICInfoSplit <- if (nrow(EICInfo) > EICBatchSize)
                split(EICInfo, ceiling(seq_len(nrow(EICInfo)) / EICBatchSize))
            else
                list(EICInfo)
            
            peaksRes <- Map(EICInfoSplit, seq_along(EICInfoSplit), f = function(EICInfoBatch, batch)
            {
                if (length(EICInfoSplit) > 1)
                    maybePrintf("Processing batch %d/%d...\n", batch, length(EICInfoSplit))
                
                EICs <- if (IMS)
                {
                    maybePrintf("Loading %d m/z+mobility EICs... ", nrow(EICInfoBatch))
                    getEICsAna(backend, EICInfoBatch, "full", genEICParams$topMostEICMZMob)
                }
                else
                {
                    maybePrintf("Loading %d m/z EICs... ", nrow(EICInfoBatch))
                    getEICsAna(backend, EICInfoBatch, "full_mz", genEICParams$topMostEICMZ)
                }
                maybePrintf("Done!\n")
                
                EICInfoBatch <- EICInfoBatch[EIC_ID %chin% names(EICs)] # omit missing
                
                maybePrintf("Finding peaks in remaining %d EICs... ", length(EICs))
                peaks <- findPeaksInEICs(EICs, peakParams, withMobility = IMS, calcStats = TRUE,
                                         assignRTWindow = assignRTWindow,
                                         sumWindowMZ = if (backend$getHaveIMS()) genEICParams$sumWindowMZ else 0,
                                         sumWindowMob = if (backend$getHaveIMS()) genEICParams$sumWindowMob else 0,
                                         logPath = file.path("log", "featEICs", paste0(ana, ".txt")), cacheDB = cacheDB)
                maybePrintf("Done! Found %d peaks.\n", nrow(peaks))
                
                peaks <- assignMZOrMobsToPeaks(peaks, EICInfoBatch, "mz")
                if (IMS)
                    peaks <- assignMZOrMobsToPeaks(peaks, EICInfoBatch, "mobility")

                peaks <- peaks[mzCentered == TRUE]
                if (IMS)
                    peaks <- peaks[mobilityCentered == TRUE]
                
                mzProfs <- if (genEICParams$saveMZProfiles)
                    getIMSProfiles(attr(EICs, "mzProfiles"), peaks, EICs, EICInfoBatch)
                EIMs <- if (genEICParams$saveEIMs)
                    getIMSProfiles(attr(EICs, "EIMs"), peaks, EICs, EICInfoBatch)
                
                return(list(peaks = peaks, mzProfs = mzProfs, EIMs = EIMs))
            })
            
            peaks <- rbindlist(lapply(peaksRes, `[[`, "peaks"), fill = TRUE) # NOTE: need to fill for empties
            dups <- findFeatTableDups(peaks$ret, peaks$retmin, peaks$retmax,  peaks$mz,
                                      if (IMS) peaks$mobility else numeric(),
                                      peaks$intensity, defaultLim("retention", "narrow"),
                                      defaultLim("mz", "medium"), defaultLim("mobility", "medium"))
            peaks <- peaks[!dups]
            
            if (genEICParams$filter == "suspects")
            {
                # only keep peaks that match with a suspect
                checkRT <- is.finite(genEICParams$rtWindow) && !is.null(suspects[["rt"]])
                checkMob <- IMS && genEICParams$filterIMS == "suspects" && !is.null(suspects[["mobility"]])
                keep <- filterPiekResults(peaks$ret, peaks$mz, if (checkMob) peaks$mobility else numeric(),
                                          if (checkRT) suspects$rt else numeric(), suspects$mz - genEICParams$mzWindow,
                                          suspects$mz + genEICParams$mzWindow,
                                          if (checkMob) suspects$mobility - genEICParams$IMSWindow else numeric(),
                                          if (checkMob) suspects$mobility + genEICParams$IMSWindow else numeric(),
                                          genEICParams$rtWindow)
                peaks <- peaks[keep == TRUE]
            }
            else if (genEICParams$filter == "ms2" && is.finite(genEICParams$rtWindow))
            {
                checkMob <- IMS && genEICParams$filterIMS == "ms2"
                mobCheckMin <- mobCheckMax <- numeric()
                if (checkMob)
                {
                    mobMeans <- mapply(MS2Info$mobmin, MS2Info$mobmax, FUN = mean)
                    mobWidths <- (MS2Info$mobmax - MS2Info$mobmin) / 2
                    mobWidths <- fifelse(mobWidths < genEICParams$IMSWindow, genEICParams$IMSWindow, mobWidths)
                    mobCheckMin <- mobMeans - mobWidths
                    mobCheckMax <- mobMeans + mobWidths
                }
                keep <- filterPiekResults(peaks$ret, peaks$mz, if (checkMob) peaks$mobility else numeric(),
                                          MS2Info$time, MS2Info$mzmin, MS2Info$mzmax, mobCheckMin, mobCheckMax,
                                          genEICParams$rtWindow)
                peaks <- peaks[keep == TRUE]
            }
            
            peaks <- removeDTColumnsIfPresent(peaks, c("mzCentered", "mobilityCentered", "keep"))

            maybePrintf("%d peaks remain after filtering.\n", nrow(peaks))            
            
            if (genEICParams$saveMZProfiles)
                mzProfiles[[ana]] <<- Reduce(c, lapply(peaksRes, `[[`, "mzProfs"))[peaks$ID]
            if (genEICParams$saveEIMs)
                EIMs[[ana]] <<- Reduce(c, lapply(peaksRes, `[[`, "EIMs"))[peaks$ID]
            
            peaks[, EIC_ID := NULL][]
            if (IMS)
                peaks[, ims_parent_ID := NA_character_]
            
            return(peaks)
        })
        
        saveCacheDataList("featuresPiek", fList[anaInfoTBD$analysis], anaHashes[anaInfoTBD$analysis], dbArg = cacheDB)
        if (genEICParams$saveMZProfiles && length(mzProfiles) > 0)
            saveCacheDataList("featuresPiekMZProfiles", mzProfiles[anaInfoTBD$analysis], anaHashes[anaInfoTBD$analysis],
                              dbArg = cacheDB)
        if (genEICParams$saveEIMs && length(EIMs) > 0)
            saveCacheDataList("featuresPiekEIMs", EIMs[anaInfoTBD$analysis], anaHashes[anaInfoTBD$analysis],
                              dbArg = cacheDB)
    }
    
    if (length(cachedData) > 0)
    {
        fList <- c(fList, cachedData)
        fList <- fList[analysisInfo$analysis] # put original order
        if (genEICParams$saveMZProfiles)
        {
            cachedMZProfiles <- pruneList(loadCacheData("featuresPiekMZProfiles", names(cachedData),
                                                        simplify = FALSE, dbArg = cacheDB))
            mzProfiles <- c(mzProfiles, cachedMZProfiles)
            mzProfiles <- mzProfiles[analysisInfo$analysis]
        }
        if (genEICParams$saveEIMs)
        {
            cachedEIMs <- pruneList(loadCacheData("featuresPiekEIMs", names(cachedData),
                                                  simplify = FALSE, dbArg = cacheDB))
            EIMs <- c(EIMs, cachedEIMs)
            EIMs <- EIMs[analysisInfo$analysis]
        }
    }
    
    if (verbose)
    {
        printf("\n-------\nDone!\n")
        printFeatStats(fList)
    }
    
    return(featuresPiek(analysisInfo = analysisInfo, features = fList, hasMobilities = IMS, mzProfiles = mzProfiles,
                        EIMs = EIMs))
}

#' @param \dots Any additional parameters to be set in the returned parameter list. These will override the defaults.
#'   See the \verb{EIC formation parameters} section for details.
#'
#' @return \code{getPiekEICParams} returns a \code{list} of parameters for the EIC generation, which is used to set
#'   the \code{genEICParams} argument to \code{findFeaturesPiek}.
#'
#' @rdname findFeaturesPiek
#' @export
getPiekEICParams <- function(...)
{
    ret <- list(filter = "none", filterIMS = "none", mzRange = c(80, 600), mzStep = 0.02, mobRange = c(0.5, 1.3),
                mobStep = 0.04, retRange = NULL, gapFactor = 3, sumWindowMZ = defaultLim("retention", "very_narrow"),
                sumWindowMob = defaultLim("retention", "very_narrow"), smoothWindowMZ = 0, smoothWindowMob = 0,
                saveMZProfiles = FALSE, saveEIMs = FALSE, minEICIntensity = 5000, minEICAdjTime = 5,
                minEICAdjPoints = 5, minEICAdjIntensity = 250, topMostEICMZ = 10000, topMostEICMZMob = 10000,
                minEICsIMSPreCheck = 50000)
    
    ret <- modifyList(ret, list(...), keep.null = TRUE)
    
    maybeSetDefault <- function(param, default)
    {
        if (is.null(ret[[param]]))
            ret[[param]] <<- default
    }
    
    maybeSetDefault("smoothExtMZ", ret$mzStep)
    maybeSetDefault("smoothExtMob", ret$mobStep)
    
    if (ret$filter != "suspects" && ret$filterIMS == "suspects")
        stop("filterIMS can only be 'suspects' if filter is also set to 'suspects'", call. = FALSE)
    if (ret$filter != "ms2" && ret$filterIMS == "ms2")
        stop("filterIMS can only be 'ms2' if filter is also set to 'ms2'", call. = FALSE)
    
    if (ret$filter == "suspects")
    {
        maybeSetDefault("rtWindow", defaultLim("retention", "medium"))
        maybeSetDefault("mzWindow", defaultLim("mz", "medium"))
        maybeSetDefault("skipInvalid", TRUE)
        maybeSetDefault("prefCalcChemProps", TRUE)
        maybeSetDefault("neutralChemProps", FALSE)
    }
    else if (ret$filter == "ms2")
    {
        maybeSetDefault("rtWindow", defaultLim("retention", "very_narrow"))
        maybeSetDefault("mzIsoWindow", 1)
        maybeSetDefault("minTIC", 10000)
    }
    
    if (ret$filterIMS != "none")
        maybeSetDefault("IMSWindow", defaultLim("mobility", "medium"))
    
    return(ret)
}
