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
        # NOTE: MS2Info is already adjusted for mz and mob tolerances
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
#'   subsequently performs automatic peak detection to generate features. The \code{piek} algorithm introduces the
#'   following improvements and changes: \itemize{
#'
#'     \item Support for IMS-HRMS workflows.
#'
#'     \item The \link{msdata} interface is used to efficiently form EICs from the raw data. All the file formats and
#'     types can be used that are supported by \link{msdata}. This includes IMS data, even if not used for feature
#'     detection, which allows the use of IMS data directly in non-IMS or \link[=assignMobilities_feat]{post mobility
#'     assignment} workflows.
#'     
#'     \item The EIC binning approach can be extended with the mobility dimension to support \link[=assignMobilities_feat]{direct mobility
#'     assignment} workflows.
#'     
#'     \item The EIC bins can be filtered with suspect or MS2 data to speed up feature detection.
#'     
#'     \item Several filters are available to eliminate EICs with are likely devoid of any signal of interest.
#'
#'     \item The original peak detection algorithm was further optimized or can be be exchanged with others: see
#'     \code{\link{getDefPeakParams}} for details.
#'
#'     \item Several filters are available to improve the data and reduce redundancy: \itemize{
#'     
#'     \item The original redundancy detection, which performs a second feature detection with EIC bins that are shifted
#'     by 50\% width and eliminates features with \code{m/z} values outside the center of any bin, was extended for IMS
#'     support.
#'
#'     \item Redundant features across bins are eliminated if with close retention time, \code{m/z}, mobility and
#'     chromatographic overlap. The most intense feature is kept.
#'     
#'     \item Data from suspects or MS2 precursors that was used to pre-filter EICs, can also be used to filter the final
#'     feature list.
#'     
#'     }
#'     
#'     \item Various small bug fixes and improvements for the original code.
#'
#'   }
#'   
#' @param genEICParams A \code{list} of parameters for the EIC generation. See the \verb{EIC generation parameters}
#'   section below. The \code{getPiekEICParams} function is used to generate the parameter list.
#' @param peakParams A \code{list} of parameters for the peak detection. See \code{\link{getDefPeakParams}} for details.
#' @param IMS Set to \code{TRUE} to use IMS data to resolve features. This initiates a
#'   \link[=assignMobilities_feat]{direct mobility assignment} workflow. If \code{IMS=FALSE} then IMS data can still be
#'   used for feature detection.
#' @param suspects The suspect list to be used for suspect pre-filtering of EIC bins. See
#'   \link[=suspect-screening]{suspect screening} for details on the suspect list format and \verb{EIC generation
#'   parameters} to enable suspect filtering.
#'   
#'   \strong{NOTE}: Suspect matching can only be performed by mobilities and not \acronym{CCS} values. The
#'   \code{\link[=assignMobilities_susp]{assignMobilities}} method should be used to convert any \acronym{CCS} data in
#'   advance.
#' @param adduct An \code{\link{adduct}} object (or something that can be converted to it with \code{\link{as.adduct}}).
#'   Examples: \code{"[M-H]-"}, \code{"[M+Na]+"}. Only needs to be specified if \code{suspects} is set.
#' @param assignMethod Should be \code{"basepeak"} or \code{"weighted.mean"}. This parameter sets how measured
#'   \emph{m/z} or mobilities across the EIC datapoints are handled for feature assignment. If
#'   \code{assignMethod="basepeak"}, then the value of the base peak (=highest intensity peak) from each EIC datapoint
#'   is taken. If \code{assignMethod="weighted.mean"} then the intensity weighted mean is calculated of the values that
#'   fall within the EIC bin.
#' @param assignRTWindow The retention time window (+/- seconds) used for aggregating EIC datapoints to assign feature
#'   \emph{m/z} and mobility data, using an intensity weighted mean. The maximum window is always bound by the feature
#'   retention time range. Increasing this number may improve accuracy by averaging more points. However, decreasing the
#'   window may reduce inaccuracies due to inclusion of data from closely eluting features (with similar \emph{m/z} and
#'   mobility) or noisy data from the chromatographic peak extremes. If \code{assignRTWindow=0} then only the EIC
#'   datapoint at the feature retention time is used.
#'
#'   The assignment window is automatically adjusted for the values set for \code{sumWindowMZ} and \code{sumWindowMob}
#'   (see \verb{EIC generation parameters}).
#' @param rtWindowDup,mzWindowDup,mobWindowDup The retention time (seconds), \emph{m/z} and mobility windows used to
#'   identify duplicate (redundant) features detected in multiple EIC bins. These values default to
#'   \code{defaultLim("retention", "very_narrow")}, \code{defaultLim("mz", "medium")} and \code{defaultLim("mobility",
#'   "medium")}, respectively (see \link{limits}).
#' @param minPeakOverlapDup The minimum overlap (fraction between 0 and 1) in retention time between two features to be
#'   considered a duplicate.
#' @param EICBatchSize The number of EICs to be processed in a single batch. Decreasing this number will reduce memory
#'   usage, at the cost of speed. Set to \code{Inf} to process all EICs in a single batch.
#' @param keepDups Set to \code{TRUE} to keep duplicate features and features with non-centered \emph{m/z} or mobility
#'   values. This is primarily intended for debugging, but can be useful to investigate why features are missing or
#'   optimize tolerance windows for duplicate feature detection.
#'
#' @inheritParams findFeatures
#' @template minIntensityIMS-arg
#'
#' @section IMS workflows: In IMS workflows (\emph{IMS=TRUE}), a 'pre-check' is performed to avoid excessive numbers of
#'   two-dimensional bins for EIC formation and peak detection. These EICs are formed by only considering the m/z
#'   dimension, and subsequently filtered by the parameters described in the \verb{EIC generation parameters} section.
#'   The final EICs for feature detection are then only formed if they have m/z data that was not removed during the
#'   pre-check.
#'
#'   The \code{m/z} and mobility data from IMS-HRMS data is typically not or partially centroided. The feature
#'   \code{m/z} and mobility values are derived from \code{m/z} or mobility \emph{versus} intensity profiles. The
#'   profiles are generated for each EIC timepoint, and the value at the maximum intensity or intensity weighted mean of
#'   the profile is used to derive the intermediate values (configured by \code{assignMethod}). Several parameters exist
#'   to improve the profile data (see next section).
#'
#' @section EIC generation parameters: The \code{genEICParams} argument to \code{findFeaturesPiek} configures the
#'   generation of EICs. The \code{getPiekEICParams} function should be used to generate the parameter list.
#'
#'   The following general parameters exist: \itemize{
#'
#'     \item \code{filter} Controls the pre-filtering of EIC bins with \emph{m/z} data. Should be \code{"none"} (no
#'     filtering), \code{"suspects"} (filter with suspect data) or \code{"ms2"} (filter with data from precursors
#'     detected in a data-dependent MS/MS experiment).
#'     
#'     \item \code{mzRange},\code{mzStep} Configures the formation of \emph{m/z} bins. \code{mzRange} is a numeric
#'     vector of length two that specifies the min/max \emph{m/z} range. \code{mzStep} specifies the bin widths.
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
#'     \item \code{topMostEICMZ} Only keep this number of top-most intense EICs. The intensity is derived from the data
#'     point with the highest intensity in the EIC. Set to zero to always select all EICs.
#'     
#'     For IMS workflows, this parameter is \emph{only} used to limit the number of EICs resulting from the 'pre-check'
#'     in the \emph{m/z} dimension.
#'
#'   }
#'   
#'   The following parameters are specifically used for IMS workflows: \itemize{
#'   
#'   \item \code{filterIMS} Similar to the \code{filter} parameter, but controls how mobility data is used for pre-filtering of EIC bins.
#'     
#'     Different values for \code{filter} and \code{filterIMS} can be specified: \itemize{
#'
#'       \item \code{filter="none"} and \code{filterIMS="none"}
#'
#'       \item \code{filter="suspects"} and \code{filterIMS="suspects"}
#'
#'       \item \code{filter="suspects"} and \code{filterIMS="none"} (only use \emph{m/z} filtering)
#'
#'       \item \code{filter="ms2"} and \code{filterIMS="ms2"}
#'
#'       \item \code{filter="ms2"} and \code{filterIMS="none"}
#'
#'     }
#'     
#'     Currently only Bruker DDA-PASEF experiments provide the data needed for \code{"ms2"} filtering.
#'   
#'     \item \code{mobRange},\code{mobStep} Equivalent to \code{mzRange} and \code{mzStep}, but for ion mobility binning.
#'   
#'     \item \code{sumWindowMZ},\code{sumWindowMob} The retention time window (+/- s) used to sum adjacent datapoints
#'     for the determination of intermediate EIC \emph{m/z} and mobility values. This data is aggregated to determine
#'     the final feature values (see also the \code{assignRTWindow} argument). Set to \samp{0} to not sum any adjacent
#'     timepoints. Larger values can generally improve accuracy for noisy data (\emph{e.g.} from TIMS), but care must be
#'     taken to stay below the expected minimum chromatographic peak width to avoid inclusion of data from other
#'     features. Defaults to \code{defaultLim("retention", "very_narrow")} (see \link{limits}).
#'     
#'     \item \code{smoothWindowMZ},\code{smoothWindowMob} The window size used to perform centered moving average
#'     smoothing on intensity data of the \emph{m/z} and mobility profiles used to determine intermediate EIC values.
#'     Smoothing of noisy data (\emph{e.g.} TIMS) is highly recommended to improve accuracy and consistency. Set to
#'     \code{0} to disable smoothing.
#'     
#'     \item \code{smoothExtMZ},\code{smoothExtMob} The \code{m/z} or mobility window to extend the smoothing at the
#'     edges of the EIC bin. This is recommended to improve smoothing, \emph{e.g.} when the peak profile is only
#'     partially captured in the bin. Defaults to the bin width, \emph{i.e.} data from an adjacent bin on each side is
#'     additionally included for smoothing. The final smoothed data is only taken from the actual EIC bin. Set to
#'     \code{0} to disable extension.
#'     
#'     \item \code{saveMZProfiles},\code{saveEIMs} Set to \code{TRUE} to save the \emph{m/z} and mobility profiles for
#'     each feature. Only the profiles at the feature retention time is saved. This can be useful for debugging or
#'     parameter optimization, but will increase memory usage and processing times.
#'     
#'     \item \code{topMostEICMob} Equivalent to \code{topMostEICMZ}, used to reduce the final two-dimensional EIC bins
#'     with \code{m/z} and mobility information.
#'     
#'     \item \code{minEICsIMSPreCheck} Only perform the \code{m/z} pre-check if the number of two-dimensional EIC bins
#'     is at least \code{minEICsIMSPreCheck}.
#'   
#'   }
#'
#'   The following parameters are specifically for when suspect data is used to pre-filter EIC bins: \itemize{
#'
#'     \item \code{rtWindow},\code{mzwindow},\code{IMSWindow}: The retention time, \emph{m/z} and mobility tolerance
#'     windows for suspect data. These are used for: \enumerate{
#'     
#'       \item Pre-filtering of EIC bins with suspect data, \emph{i.e.} larger tolerances will lead to more EIC bins
#'       being kept. (only applicable for \code{mzWindow} and \code{IMSWindow}).
#'       
#'       \item Matching the final features to suspect data. \code{rtWindow=Inf} can be used to disable retention time
#'       matching.
#'     
#'     }
#'     
#'     Defaults to \code{defaultLim("retention", "medium")}, \code{defaultLim("mz", "medium")} and
#'     \code{defaultLim("mobility", "medium")}, see \link{limits}.
#'
#'     \item \code{skipInvalid},\code{prefCalcChemProps},\code{neutralChemProps} Controls preparing the suspect list
#'     data. See \code{\link{screenSuspects}}.
#'
#'   }
#'
#'   The following parameters are specifically for when MS2 data is used to pre-filter EICs: \itemize{
#'
#'     \item \code{rtWindow} Eliminates any features without an MS/MS spectrum within this retention time window. Set
#'     \code{rtWindow=Inf} to disable this filter. Defaults to \code{defaultLim("retention", "very_narrow")} (see
#'     \link{limits}).
#'     
#'     \item \code{mzIsoWindow} The maximum \emph{m/z} window considered for MS/MS precursors that were isolated by DDA.
#'     These \emph{m/z} isolation windows are used to pre-filter EICs and match the final features. Setting
#'     \code{mzIsoWindow} to a value lower than typical instrument isolation windows will make feature detection more
#'     specific, as features need to be more close to the triggered DDA precursor \code{m/z} values. In contrast, larger
#'     values for \code{mzIsoWindow} allows to include features that were not specifically targeted by DDA, but may
#'     still have MS/MS data as their \emph{m/z} could still fall within the MS/MS isolation window. The effective
#'     window used will never exceed the instrumental isolation window. Setting \code{mzIsoWindow=Inf} will always use
#'     instrumental windows.
#'     
#'     \item \code{IMSWindow} The mobility tolerance window to match DDA MS/MS precursors in IMS workflows. Used for
#'     pre-filtering EICs and the final features. To match DDA precursor data, the measured mobility range of the
#'     corresponding MS/MS data is used as the mobility window. This window is then adjusted to be at least +/-
#'     \code{IMSWindow}. Defaults to \code{defaultLim("mobility", "medium")} (see \link{limits})
#'
#'     \item \code{minTIC} The minimum total ion current (TIC) signal for an MS/MS spectrum to be considered. Can be
#'     increased to eliminate features with low intensity MS/MS data.
#'
#'   }
#'
#' @templateVar what \code{findFeaturesPiek}
#' @templateVar noProfile TRUE
#' @template uses-msdata
#'
#' @references \insertAllCited{}
#'
#' @inherit findFeatures return
#'
#' @export
findFeaturesPiek <- function(analysisInfo, genEICParams, peakParams, IMS = FALSE, suspects = NULL, adduct = NULL,
                             assignMethod = "basepeak", assignRTWindow = defaultLim("retention", "very_narrow"),
                             rtWindowDup = defaultLim("retention", "narrow"),
                             mzWindowDup = defaultLim("mz", "medium"),
                             mobWindowDup = defaultLim("mobility", "medium"), minPeakOverlapDup = 0.25,
                             minIntensityIMS = 25, EICBatchSize = Inf, keepDups = FALSE, verbose = TRUE)
{
    # UNDONE: add refs to docs, and highlight changes
    # UNDONE: test empties, eg no EICs
    
    if (!is.null(adduct))
        adduct <- checkAndToAdduct(adduct)
    
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    assertPiekGenEICParams(genEICParams, add = ac)
    assertFindPeakParams(peakParams, add = ac)
    checkmate::assertFlag(IMS, add = ac)
    if (genEICParams$filter == "suspects")
        assertSuspectList(suspects, needsAdduct = is.null(adduct), skipInvalid = genEICParams$skipInvalid, null.ok = FALSE,
                          add = ac)
    aapply(checkmate::assertNumber, . ~ minIntensityIMS + assignRTWindow + rtWindowDup + mzWindowDup, mobWindowDup,
           lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::assertChoice(assignMethod, c("basepeak", "weighted.mean"), add = ac)
    checkmate::assertNumber(minPeakOverlapDup, lower = 0, upper = 1, finite = TRUE, add = ac)
    if (!is.infinite(EICBatchSize))
        checkmate::assertCount(EICBatchSize, positive = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ keepDups + verbose, fixed = list(add = ac))
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
                         "Please note that CCS data first must be converted to mobilities, e.g. with assignMobilities(), ",
                         "and check if the adduct argument is correcttly set for the suspect list.",
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
    baseHash <- makeHash(genEICParams, peakParams, IMS, suspects, adduct, assignMethod, assignRTWindow, rtWindowDup,
                         mzWindowDup, mobWindowDup, minPeakOverlapDup, minIntensityIMS, keepDups)
    anaHashes <- getMSFileHashesFromAvailBackend(analysisInfo, needTypes = if (IMS) "ims" else c("ims", "centroid"))
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
                          EICInfo$mobmin, EICInfo$mobmax, gapFactor = genEICParams$gapFactor,
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
        chkCol <- paste0(what, "BP")
        minCol <- if (what == "mz") "mzmin" else "mobmin"
        centCol <- paste0(what, "Centered")
        binMinCol <- paste0("binMin", if (what == "mz") "MZ" else "Mob")
        binMaxCol <- paste0("binMax", if (what == "mz") "MZ" else "Mob")
        step <- if (what == "mz") genEICParams$mzStep else genEICParams$mobStep
        peaks[, (binMinCol) := EICInfo[match(peaks$EIC_ID, EIC_ID)][[minCol]]]
        peaks[, (binMaxCol) := get(binMinCol) + step]
        peaks[, (centCol) := between(get(chkCol), get(binMinCol) + step/4, get(binMinCol) + step/4*3)]
        if (assignMethod == "basepeak")
            peaks[, (what) := get(paste0(what, "BP"))]
        # else already weighted mean
        
        # NOTE: centered is kept and removed later to allow filtering redundant peaks
        peaks[, c(paste0(what, "BP")) := NULL]
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
        fList <- applyMSData(anaInfoTBD, needTypes = if (IMS) "ims" else c("ims", "centroid"), showProgress = FALSE, func = function(ana, path, backend)
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
                
                # reduce m/z tolerance if needed
                if (!is.finite(genEICParams$mzIsoWindow))
                {
                    MS2Info[, c("mzmin", "mzmax") := .(precursorMZ - isolationRangeMin,
                                                       precursorMZ + isolationRangeMax)]
                }
                else
                {
                    MS2Info[, c("mzmin", "mzmax") := .(precursorMZ - pmin(genEICParams$mzIsoWindow, isolationRangeMin),
                                                       precursorMZ + pmin(genEICParams$mzIsoWindow, isolationRangeMax))]
                }
                MS2Info <- MS2Info[TIC >= genEICParams$minTIC]
                if (IMS)
                {
                    setnames(MS2Info, c("mobStart", "mobEnd"), c("mobmin", "mobmax"))
                    # grow mob tolerance if needed
                    MS2Info[, mobCenter := (mobmin + mobmax) / 2]
                    MS2Info[, mobmin := pmin(mobmin, mobCenter - genEICParams$IMSWindow)]
                    MS2Info[, mobmax := pmax(mobmax, mobCenter + genEICParams$IMSWindow)]
                }
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

                if (!keepDups)
                {
                    peaks <- peaks[mzCentered == TRUE]
                    if (IMS)
                        peaks <- peaks[mobilityCentered == TRUE]
                }
                
                mzProfs <- if (genEICParams$saveMZProfiles)
                    getIMSProfiles(attr(EICs, "mzProfiles"), peaks, EICs, EICInfoBatch)
                EIMs <- if (genEICParams$saveEIMs)
                    getIMSProfiles(attr(EICs, "EIMs"), peaks, EICs, EICInfoBatch)
                
                return(list(peaks = peaks, mzProfs = mzProfs, EIMs = EIMs))
            })
            
            peaks <- rbindlist(lapply(peaksRes, `[[`, "peaks"), fill = TRUE) # NOTE: need to fill for empties
            # NOTE: only do centered here in case keepDups==T
            peaksCentered <- peaks[mzCentered == TRUE]
            if (IMS)
                peaksCentered <- peaksCentered[mobilityCentered == TRUE]
            dups <- findFeatTableDups(peaksCentered$ret, peaksCentered$retmin, peaksCentered$retmax,
                                      peaksCentered$mz,
                                      if (IMS) peaksCentered$mobility else numeric(),
                                      peaksCentered$intensity, rtWindowDup, mzWindowDup, mobWindowDup,
                                      minPeakOverlapDup)
            if (!keepDups)
                peaks <- peaks[dups == 0] # peaks == peaksCentered
            else
            {
                dupsID <- rep(NA_character_, length(dups))
                dupsID[dups != 0] <- peaksCentered$ID[dups[dups != 0]]
                peaks[peaksCentered, dup := dupsID, on = "ID"]
            }
            
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
                keep <- filterPiekResults(peaks$ret, peaks$mz, if (checkMob) peaks$mobility else numeric(),
                                          MS2Info$time, MS2Info$mzmin, MS2Info$mzmax, MS2Info$mobmin, MS2Info$mobmax,
                                          genEICParams$rtWindow)
                peaks <- peaks[keep == TRUE]
            }
            
            peaks <- removeDTColumnsIfPresent(peaks, "keep")
            if (!keepDups)
                peaks <- removeDTColumnsIfPresent(peaks, c("mzCentered", "mobilityCentered"))

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
    
    return(featuresPiek(analysisInfo = analysisInfo, features = fList, hasMobilities = IMS, fromIMS = IMS,
                        mzProfiles = mzProfiles, EIMs = EIMs))
}

#' @param \dots Any additional parameters to be set in the returned parameter list. These will override the defaults.
#'   See the \verb{EIC generation parameters} section for details.
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
