#' @include features.R
NULL

removeDuplicateFeatsSusps <- function(tab, checkRet, selectTopIntens, nameCol)
{
    # marks duplicates in feature or suspect table, and keeps the feature with the highest intensity

    if (nrow(tab) <= 1)
        return(tab)

    isSame <- function(x, y)
    {
        yes <- numLTE(abs(x$mz - y$mz), defaultLim("mz", "very_narrow"))
        if (checkRet)
            yes <- yes & ((is.na(x$ret) & is.na(y$ret)) | numLTE(abs(x$ret - y$ret), defaultLim("retention", "very_narrow")))
        if (!is.null(x[["mobility"]]))
            yes <- yes & numLTE(abs(x$mobility - y$mobility), defaultLim("mobility", "very_narrow"))
        return(yes)
    }
    
    tab <- copy(tab)
    tab[, duplicate := FALSE]
    for (r in seq_len(nrow(tab) - 1))
    {
        if (tab$duplicate[r])
            next
        rTab <- tab[r]
        nextTab <- tab[seq(r + 1, nrow(tab))]
        nextTab <- nextTab[duplicate == FALSE & isSame(rTab, nextTab)]
        if (nrow(nextTab) > 0)
        {
            if (selectTopIntens)
            {
                maxInt <- max(rTab$intensity, nextTab$intensity)
                tab[get(nameCol) %chin% c(rTab[[nameCol]], nextTab[[nameCol]]), duplicate := !numEQ(intensity, maxInt)]
            }
            else
            {
                # otherwise just keep current and mark all from next rows
                tab[get(nameCol) %chin% nextTab[[nameCol]], duplicate := TRUE]
            }
        }
    }
    return(tab[duplicate == FALSE][, duplicate := NULL])
}

getPiekEICsInfo <- function(params, suspects, withIMS, MS2Info)
{
    EICInfo <- NULL
    if (params$methodMZ == "bins")
    {
        # UNDONE: also support other binning approaches?
        binsMZ <- seq(params$mzRange[1], params$mzRange[2], by = params$mzStep * 0.5)
        names(binsMZ) <- paste0("bin_M", binsMZ)
        
        EICInfo <- data.table(mzmin = binsMZ, mzmax = binsMZ + params$mzStep, EIC_ID_MZ = names(binsMZ))
    }
    else if (params$methodMZ == "suspects")
        EICInfo <- data.table(mzmin = suspects$mz - params$mzWindow, mzmax = suspects$mz + params$mzWindow)
    else if (params$methodMZ == "ms2")
        EICInfo <- data.table(mzmin = MS2Info$mz - params$mzWindow, mzmax = MS2Info$mz + params$mzWindow)
    
    if (withIMS)
    {
        if (params$methodIMS == "bins")
        {
            binsIMS <- NULL
            {
                binsIMS <- seq(params$mobRange[1], params$mobRange[2], by = params$mobStep * 0.5)
                names(binsIMS) <- paste0("bin_I", binsIMS)
            }
            tab <- CJ(EIC_ID_MZ = names(binsMZ), EIC_ID_IMS = names(binsIMS), sorted = FALSE)
            tab[, c("mobmin", "mobmax") := .(binsIMS[EIC_ID_IMS], binsIMS[EIC_ID_IMS] + params$mobStep)]
            EICInfo <- merge(EICInfo, tab, by = "EIC_ID_MZ", sort = FALSE)
        }
        else if (params$methodIMS == "suspects")
        {
            EICInfo <- EICInfo[, c("mobmin", "mobmax") := .(suspects$mobility - params$IMSWindow,
                                                            suspects$mobility + params$IMSWindow)]
        }
        else if (params$methodIMS == "ms2")
        {
            EICInfo <- EICInfo[, c("mobmin", "mobmax") := .(MS2Info$mobility - params$IMSWindow,
                                                            MS2Info$mobility + params$IMSWindow)]
        }
    }
    else
        EICInfo[, c("mobmin", "mobmax") := 0]
    
    EICInfo[, EIC_ID := paste0("EIC_", .I)]
    
    return(EICInfo)
}

#' @rdname features-class
#' @export
featuresPiek <- setClass("featuresPiek", contains = "features")

setMethod("initialize", "featuresPiek",
          function(.Object, ...) callNextMethod(.Object, algorithm = "piek", ...))

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
#'   section below. The \code{getPiekGenEICParams} is used to generate the parameter list.
#' @param peakParams A \code{list} of parameters for the peak detection. See \code{\link{getDefPeakParams}} for details.
#' @param suspects A suspect list that needs to be specified if EICs are formed from suspect data. See
#'   \link[=suspect-screening]{suspect screening} for details on the suspect list format.
#' @param adduct An \code{\link{adduct}} object (or something that can be converted to it with \code{\link{as.adduct}}).
#'   Examples: \code{"[M-H]-"}, \code{"[M+Na]+"}. Only needs to be specified when EICs are formed from suspect data.
#'
#' @inheritParams findFeatures
#' @template minIntensityIMS-arg
#'
#' @section EIC formation parameters: The \code{genEICParams} argument to \code{findFeaturesPiek} configures the
#'   formation of EICs. The \code{getPiekGenEICParams} function should be used to generate the parameter list.
#'
#'   The following general parameters exist: \itemize{
#'
#'     \item \code{methodMZ} Sets how m/z data is used for the formation of EICs. Possible values are \code{"bins"} (use
#'     equally sized m/z bins), \code{"suspects"} (use the unique m/z values from a suspect list) and \code{"ms2"} (use
#'     the m/z values from the precursors detected in a data-dependent MS/MS experiment).
#'
#'     \item \code{methodIMS} Equivalent as \code{methodMZ}, but for ion mobility data. If \code{NULL}, no IMS data will
#'     be used and no \link[=assignMobilities_feat]{direct mobility assignment} IMS workflow is initiated. Different
#'     values for \code{methodMZ} and \code{methodIMS} can be specified to combine approaches. The following
#'     combinations are supported: \itemize{
#'
#'       \item \code{methodMZ="bins"} and \code{methodIMS="bins"}
#'
#'       \item \code{methodMZ="suspects"} and \code{methodIMS="suspects"}
#'
#'       \item \code{methodMZ="suspects"} and \code{methodIMS="bins"}
#'
#'       \item \code{methodMZ="ms2"} and \code{methodIMS="ms2"}
#'
#'       \item \code{methodMZ="ms2"} and \code{methodIMS="bins"}
#'
#'     }
#'
#'     Currently only Bruker DDA-PASEF experiments provide the data needed for the \code{"ms2"} approach.
#'
#'     \item \code{retRange} A \code{numeric} vector of length two that specifies the retention time range for the EICs.
#'     Data outside this range is excluded. Set to \code{NULL} to use the full range.
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
#'     \code{methodIMS="bins"}).
#'
#'   }
#'
#'   The following parameters are specifically for EICs from suspect data: \itemize{
#'
#'     \item \code{rtWindow},\code{mzwindow},\code{IMSWindow}: If retention times are present in the suspect list:
#'     specify the retention time, m/z and mobility (if \code{methodIMS="suspects"}) tolerance window to match features
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
#'     \item \code{mzWindow},\code{IMSWindow} The m/z and mobility (if \code{methodIMS="ms2"}) tolerance windows, used to \enumerate{
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
#'   This is especially important when both \code{methodMZ} and \code{methodIMS} are set to \code{"bins"}, as this may
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
findFeaturesPiek <- function(analysisInfo, genEICParams, peakParams, suspects = NULL, adduct = NULL,
                             minIntensityIMS = 25, verbose = TRUE)
{
    # UNDONE: add refs to docs, and highlight changes
    # UNDONE: use BP intensity?
    # UNDONE: test empties, eg no EICs
    
    if (!is.null(adduct))
        adduct <- checkAndToAdduct(adduct)
    
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    assertPiekGenEICParams(genEICParams, add = ac)
    assertFindPeakParams(peakParams, add = ac)
    if (genEICParams$methodMZ == "suspects")
        assertSuspectList(suspects, needsAdduct = is.null(adduct), skipInvalid = genEICParams$skipInvalid, null.ok = TRUE,
                          add = ac)
    checkmate::assertNumber(minIntensityIMS, lower = 0, finite = TRUE, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    maybePrintf <- \(...) if (verbose) printf(...)
    
    if (genEICParams$methodMZ == "suspects")
    {
        suspects <- prepareSuspectList(suspects, adduct = adduct, skipInvalid = genEICParams$skipInvalid,
                                       checkDesc = TRUE, prefCalcChemProps = genEICParams$prefCalcChemProps,
                                       neutralChemProps = genEICParams$neutralChemProps)
        suspectsOrig <- suspects
        suspects <- removeDuplicateFeatsSusps(suspects, FALSE, FALSE, "name")
        suspects <- suspects[order(mz)]
        suspsRM <- setdiff(suspectsOrig$name, suspects$name)
        maybePrintf("The following %d non-unique suspects were removed %s.\n", length(suspsRM),
                    getStrListWithMax(suspsRM, 10, ", "))
    }
    
    if (is.null(genEICParams[["retRange"]]))
        genEICParams$retRange <- c(0, 0)
    
    withIMS <- !is.null(genEICParams[["methodIMS"]])
    
    cacheDB <- openCacheDBScope()
    baseHash <- makeHash(genEICParams, peakParams, minIntensityIMS)
    anaHashes <- getMSFileHashesFromAvailBackend(analysisInfo, needIMS = withIMS)
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
        args <- list(backend, EICInfo$mzmin, EICInfo$mzmax, genEICParams$retRange[1], genEICParams$retRange[2],
                     EICInfo$mobmin, EICInfo$mobmax, mzExpIMSWindow = 0, minIntensityIMS = minIntensityIMS,
                     mode = mode, minEICIntensity = genEICParams$minEICIntensity,
                     minEICAdjTime = genEICParams$minEICAdjTime,
                     minEICAdjPoints = genEICParams$minEICAdjPoints,
                     minEICAdjIntensity = genEICParams$minEICAdjIntensity, topMost = topMost)
        ret <- do.call(if (mode == "test") getEICList else doGetEICsForAna, args)
        names(ret) <- EICInfo$EIC_ID
        if (mode != "test")
            ret <- pruneList(ret, checkZeroRows = TRUE, keepAttr = TRUE)
        return(ret)
    }
    
    fList <- list()
    if (nrow(anaInfoTBD) > 0)
    {
        fList <- applyMSData(anaInfoTBD, needIMS = withIMS, showProgress = FALSE, func = function(ana, path, backend)
        {
            maybePrintf("\n-------\nFinding features for '%s'...\n", ana)
            
            openMSReadBackend(backend, path)
         
            MS2Info <- NULL
            if (genEICParams$methodMZ == "ms2")
            {
                MS2Info <- if (identical(genEICParams$methodIMS, "ms2"))
                    getIsolationMZsAndMobs(backend, genEICParams$clusterMethod, genEICParams$mzWindow,
                                           genEICParams$IMSWindow, genEICParams$minTIC)
                else
                    getIsolationMZs(backend, genEICParams$clusterMethod, genEICParams$mzWindow, genEICParams$minTIC)
                setDT(MS2Info)
            }
            
            EICs <- EICInfo <- NULL
            EICInfoMZ <- getPiekEICsInfo(genEICParams, suspects, withIMS = FALSE, MS2Info = MS2Info)
            if (withIMS)
            {
                maybePrintf("Pre-checking %d m/z EICs... ", nrow(EICInfoMZ))
                testEICs <- getEICsAna(backend, EICInfoMZ, "test", genEICParams$topMostEICPre)
                testEICs <- unlist(testEICs)
                EICInfoMZ <- EICInfoMZ[EIC_ID %chin% names(testEICs)[testEICs]]
                maybePrintf("Done! Eliminated %d (%.2f%%) EICs\n", sum(!testEICs),
                            (sum(!testEICs) / length(testEICs)) * 100)
                
                # remove complete m/z bins that were filtered out before
                temp <- EICInfoMZ[, c("mzmin", "mzmax"), with = FALSE]
                setkeyv(temp, c("mzmin", "mzmax"))
                EICInfoMob <- getPiekEICsInfo(genEICParams, suspects, withIMS = TRUE, MS2Info = MS2Info)
                ov <- foverlaps(EICInfoMob, temp, type = "within", nomatch = NULL, which = TRUE)
                EICInfo <- EICInfoMob[ov$xid]
                
                maybePrintf("Loading %d m/z+mobility EICs... ", nrow(EICInfo))
                EICs <- getEICsAna(backend, EICInfo, "full", genEICParams$topMostEIC)
                maybePrintf("Done!\n")
            }
            else
            {
                maybePrintf("Loading %d m/z EICs... ", nrow(EICInfoMZ))
                EICs <- getEICsAna(backend, EICInfoMZ, "full_mz", genEICParams$topMostEIC)
                EICInfo <- EICInfoMZ[EIC_ID %chin% names(EICs)] # omit missing
                maybePrintf("Done!\n")
            }
            
            maybePrintf("Finding peaks in %d remaining EICs... ", length(EICs))
            peaks <- findPeaksInEICs(EICs, peakParams, withMobility = withIMS, calcStats = TRUE,
                                     logPath = file.path("log", "featEICs", paste0(ana, ".txt")), cacheDB = cacheDB)
            maybePrintf("Done! Found %d peaks.\n", nrow(peaks))

            # only keep those peaks with m/z in the "center" of the analyzed m/z and mobility range
            if (genEICParams$methodMZ == "bins")
            {
                peaks[, binMZStart := EICInfo[match(peaks$EIC_ID, EIC_ID)]$mzmin]
                peaks <- peaks[between(mz, binMZStart + genEICParams$mzStep/4, binMZStart + genEICParams$mzStep/4*3) == TRUE]
            }
            if (identical(genEICParams[["methodIMS"]], "bins"))
            {
                peaks[, binMobStart := EICInfo[match(peaks$EIC_ID, EIC_ID)]$mobmin]
                peaks <- peaks[between(mobility, binMobStart + genEICParams$mobStep/4, binMobStart + genEICParams$mobStep/4*3) == TRUE]
            }
            if (genEICParams$methodMZ == "suspects" && is.finite(genEICParams$rtWindow) && !is.null(suspects[["rt"]]))
            {
                # only keep peaks with at least one closely eluting suspect
                peaks[, keep := {
                    susp <- suspectsOrig[numLTE(abs(mz - omz), genEICParams$mzWindow) &
                                             numLTE(abs(rt - ort), genEICParams$rtWindow),
                                         env = list(omz = mz, ort = ret)]
                    if (identical(genEICParams[["methodIMS"]], "ms2") && !is.null(suspects[["mobility"]]))
                        susp <- susp[numLTE(abs(mobility - omob), genEICParams$IMSWindow), env = list(omob = mobility)]
                    nrow(susp) > 0
                }, by = .I]
                peaks <- peaks[keep == TRUE]
            }
            else if (genEICParams$methodMZ %in% "ms2" && is.finite(genEICParams$rtWindow))
            {
                # only keep peaks that elute closely to at least one MS2 spectrum
                peaks[, keep := {
                    MS2Peaks <- MS2Info[numLTE(abs(mz - omz), genEICParams$mzWindow), env = list(omz = mz)]
                    if (identical(genEICParams[["methodIMS"]], "ms2"))
                        MS2Peaks <- MS2Peaks[numLTE(abs(mobility - omob), genEICParams$IMSWindow), env = list(omob = mobility)]
                    allRet <- unlist(MS2Peaks$times)
                    any(numLTE(abs(ret - allRet), genEICParams$rtWindow))
                }, by = .I]
                peaks <- peaks[keep == TRUE]
            }
            
            peaks <- removeDuplicateFeatsSusps(peaks, TRUE, TRUE, "ID")
            peaks <- removeDTColumnsIfPresent(peaks, c("binMZStart", "binMobStart", "keep"))

            maybePrintf("%d peaks remain after filtering.\n", nrow(peaks))            
            
            return(peaks)
        })
        
        for (a in anaInfoTBD$analysis)
            saveCacheData("featuresPiek", fList[[a]], anaHashes[[a]], dbArg = cacheDB)
    }
    
    if (length(cachedData) > 0)
    {
        fList <- c(fList, cachedData)
        fList <- fList[analysisInfo$analysis] # put original order
    }
    
    if (verbose)
    {
        printf("\n-------\nDone!\n")
        printFeatStats(fList)
    }
    
    return(featuresPiek(analysisInfo = analysisInfo, features = fList, hasMobilities = withIMS))
}

#' @param methodMZ,methodIMS Sets the \code{methodMZ} and \code{methodIMS} parameters for the EIC generation. See the
#'   \verb{EIC formation parameters} section for details.
#' @param \dots Any additional parameters to be set in the returned parameter list. These will override the defaults.
#'   See the \verb{EIC formation parameters} section for details.
#'
#' @return \code{getPiekGenEICParams} returns a \code{list} of parameters for the EIC generation, which is used to set
#'   the \code{genEICParams} argument to \code{findFeaturesPiek}.
#'
#' @rdname findFeaturesPiek
#' @export
getPiekGenEICParams <- function(methodMZ, methodIMS = NULL, ...)
{
    checkmate::assertChoice(methodMZ, c("bins", "suspects", "ms2"))
    checkmate::assertChoice(methodIMS, c("bins", "suspects", "ms2"), null.ok = TRUE)
    
    if (methodMZ != "suspects" && identical(methodIMS, "suspects"))
        stop("methodIMS can only be 'suspects' if methodMZ is also set to 'suspects'", call. = FALSE)
    if (methodMZ != "ms2" && identical(methodIMS, "ms2"))
        stop("methodIMS can only be 'ms2' if methodMZ is also set to 'ms2'", call. = FALSE)
    
    ret <- list(methodMZ = methodMZ, methodIMS = methodIMS, retRange = NULL, minEICIntensity = 5000,
                minEICAdjTime = 5, minEICAdjPoints = 5, minEICAdjIntensity = 250, topMostEIC = 10000,
                topMostEICPre = 10000)
    
    if (methodMZ == "bins")
    {
        ret <- modifyList(ret, list(
            mzRange = c(80, 600),
            mzStep = 0.02
        ))
    }
    else if (methodMZ == "suspects")
    {
        ret <- modifyList(ret, list(
            rtWindow = defaultLim("retention", "medium"),
            mzWindow = defaultLim("mz", "medium"),
            # UNDONE these to separate param and also use elsewhere?
            skipInvalid = TRUE,
            prefCalcChemProps = TRUE,
            neutralChemProps = FALSE
        ))
    }
    else if (methodMZ == "ms2")
    {
        ret <- modifyList(ret, list(
            rtWindow = defaultLim("retention", "very_narrow"),
            mzWindow = defaultLim("mz", "narrow"),
            minTIC = 10000,
            clusterMethod = "distance"
        ))
    }
    
    if (!is.null(methodIMS))
    {
        if (methodIMS == "bins")
        {
            ret <- modifyList(ret, list(
                mobRange = c(0.5, 1.3),
                mobStep = 0.04
            ))
        }
        else if (methodIMS == "suspects")
        {
            ret <- modifyList(ret, list(
                IMSWindow = defaultLim("mobility", "medium")
            ))
        }
        else if (methodIMS == "ms2")
        {
            ret <- modifyList(ret, list(
                IMSWindow = defaultLim("mobility", "medium") # UNDONE: or narrow? If yes, update docs
            ))
        }
    }
    
    return(modifyList(ret, list(...), keep.null = TRUE))
}
