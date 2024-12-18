# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

# NOTE: No coverage calculation for Bruker tools as they cannot be run on CI
# nocov start

# haven't found a way to access constants via RCOM, so just list them here (obtained with visual studio)
DAConstants <- list(
    daMzXML = 8, daMzData = 9, daMzML = 10,
    daLine = 1, daProfile = 2,
    daTIC = 1, daEIC = 2, daBPC = 3,
    daAll = -3, daAllMS = -1, daAllMSn = -2,
    daMSFilterAll = 0, daMSFilterAllMSMS = 2, daMSFilterBBCID = 9, daMSFilterMS = 1, daMSFilterMSMS = 3,
    daBoth = 3, daNegative = 2, daPositive = 1,
    daBgrdTypeNone = 0, daBgrdTypeSpectral = 3,
    daCSV = 2, daASCII = 4
)

# DA needs backslashes
getBrukerAnalysisPath <- function(analysis, path) normalizePath(file.path(path, paste0(analysis, ".d")))

hideDAInScope <- withr::local_(function(x) getDAApplication()$Hide(), function(x) getDAApplication()$Show())

#' @details \code{showDataAnalysis} makes a hidden DataAnalysis window visible
#'   again. Most functions using DataAnalysis will hide the window during
#'   processing for efficiency reasons. If the window remains hidden
#'   (\emph{e.g.} because there was an error) this function can be used to make
#'   it visible again. This function can also be used to start DataAnalysis if
#'   it is not running yet.
#' @rdname bruker-utils
#' @export
showDataAnalysis <- function() getDAApplication()$Show()

getDAMaxIntMZAndFWHM <- function(spec)
{
    plist <- spec[["MSPeakList"]]
    psize <- plist[["Count"]]

    maxInt <- 0
    ret <- list(mz = 0, fwhm = 0)

    if (psize > 0)
    {
        for (i in seq_len(psize))
        {
            int <- plist[[i]][["Intensity"]]
            if (int > maxInt)
            {
                ret$mz <- plist[[i]][["m_over_z"]]
                ret$fwhm <- plist[[i]][["width"]]
                maxInt <- int
            }
        }
    }

    return(ret)
}

getDAApplication <- function()
{
    checkPackage("RDCOMClient")
    return(RDCOMClient::COMCreate("BDal.DataAnalysis.Application"))
}

# if path=NULL, analysis should be complete file path
getDAFileIndex <- function(DA, analysis, path, openFileIfClosed = TRUE)
{
    if (!is.null(path))
        analysis <- getBrukerAnalysisPath(analysis, path)

    getFIndex <- function()
    {
        openCount <- DA[["Analyses"]][["Count"]]

        for (i in seq_len(openCount))
        {
            if (DA[["Analyses"]][[i]][["Path"]] == analysis)
                return(i)
        }

        return(-1)
    }

    ret <- getFIndex()
    if (ret == -1 && openFileIfClosed) # File is not yet open?
    {

        DA[["Analyses"]]$Open(analysis)
        ret <- getFIndex()
    }

    if (ret == -1 && openFileIfClosed)
        warning(sprintf("Failed to open '%s'", analysis))

    return(ret)
}

closeSaveDAFile <- function(DA, DAFind, close, save)
{
    if (DAFind != -1)
    {
        if (save)
            DA[["Analyses"]][[DAFind]]$Save()
        if (close)
            DA[["Analyses"]][[DAFind]]$Close()
    }
    invisible(NULL)
}

#' @details \code{setDAMethod} Sets a given DataAnalysis method (\file{.m} file)
#'   to a set of analyses. \strong{NOTE}: as a workaround for a bug in
#'   DataAnalysis, this function will save(!), close and re-open any analyses
#'   that are already open prior to setting the new method. The \code{close}
#'   argument only controls whether the file should be closed after setting the
#'   method (files are always saved).
#'
#' @param method The full path of the DataAnalysis method.
#'
#' @rdname bruker-utils
#' @export
setDAMethod <- function(anaInfo, method, close = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    anaInfo <- assertAndPrepareAnaInfo(anaInfo, "bruker", add = ac)
    checkmate::assertDirectoryExists(method, add = ac)
    checkmate::assertFlag(close, add = ac)
    checkmate::reportAssertions(ac)

    DA <- getDAApplication()
    method <- normalizePath(method)
    hideDAInScope()

    for (i in seq_len(nrow(anaInfo)))
    {
        printf("Setting DA method of analysis '%s' to %s (%d/%d)...\n", anaInfo$analysis[i], method, i, nrow(anaInfo))

        # HACK: need to close file (if open) otherwise method is not applied
        ind <- getDAFileIndex(DA, anaInfo$analysis[i], anaInfo$path[i], openFileIfClosed = FALSE)
        closeSaveDAFile(DA, ind, TRUE, TRUE)

        ind <- getDAFileIndex(DA, anaInfo$analysis[i], anaInfo$path[i])
        if (ind == -1)
            next
        DA[["Analyses"]][[ind]]$LoadMethod(method)
        closeSaveDAFile(DA, ind, close, TRUE)
    }

    invisible(NULL)
}

#' @details \code{revertDAAnalyses} Reverts a given set of analyses to their
#'   unprocessed raw state.
#' @rdname bruker-utils
#' @export
revertDAAnalyses <- function(anaInfo, close = TRUE, save = close)
{
    ac <- checkmate::makeAssertCollection()
    anaInfo <- assertAndPrepareAnaInfo(anaInfo, "bruker", add = ac)
    assertDACloseSaveArgs(close, save, add = ac)
    checkmate::reportAssertions(ac)

    DA <- getDAApplication()
    hideDAInScope()

    for (i in seq_len(nrow(anaInfo)))
    {
        printf("Reverting DA analysis '%s' (%d/%d)...\n", anaInfo$analysis[i], i, nrow(anaInfo))

        ind <- getDAFileIndex(DA, anaInfo$analysis[i], anaInfo$path[i])
        if (ind == -1)
            next
        DA[["Analyses"]][[ind]]$LoadRawData()
        closeSaveDAFile(DA, ind, close, save)
    }

    invisible(NULL)
}

#' @details \code{recalibrarateDAFiles} Performs automatic mass recalibration of
#'   a given set of analyses. The current method settings for each analyses will
#'   be used.
#'
#' @rdname bruker-utils
#' @export
recalibrarateDAFiles <- function(anaInfo, close = TRUE, save = close)
{
    ac <- checkmate::makeAssertCollection()
    anaInfo <- assertAndPrepareAnaInfo(anaInfo, "bruker")
    assertDACloseSaveArgs(close, save, add = ac)
    checkmate::reportAssertions(ac)

    DA <- getDAApplication()
    hideDAInScope()

    for (i in seq_len(nrow(anaInfo)))
    {
        printf("Recalibrating analysis '%s' (%d/%d)... ", anaInfo$analysis[i], i, nrow(anaInfo))

        ind <- getDAFileIndex(DA, anaInfo$analysis[i], anaInfo$path[i])
        if (ind == -1 || !DA[["Analyses"]][[ind]]$RecalibrateAutomatically())
        {
            warning(sprintf("Failed to calibrate analysis '%s'", anaInfo$analysis[i]))
            cat("Failed!!\n")
        }
        else
        {
            # UNDONE: assume there is only one calibration item
            printf("Done! (std: %f ppm)\n", DA[["Analyses"]][[ind]][["CalibrationStatus"]][[1]]$StandardDeviation())
            closeSaveDAFile(DA, ind, close, save)
        }
    }

    invisible(NULL)
}

#' @details \code{getDACalibrationError} is used to obtain the standard
#'   deviation of the current mass calibration (in ppm).
#'
#' @return \code{getDACalibrationError} returns a \code{data.frame} with a
#'   column of all analyses (named \code{analysis}) and their mass error (named
#'   \code{error}).
#'
#' @rdname bruker-utils
#' @export
getDACalibrationError <- function(anaInfo)
{
    anaInfo <- assertAndPrepareAnaInfo(anaInfo, "bruker")

    DA <- getDAApplication()
    hideDAInScope()

    ret <- data.frame(analysis = anaInfo$analysis, error = numeric(nrow(anaInfo)), stringsAsFactors = FALSE)

    for (i in seq_len(nrow(anaInfo)))
    {
        printf("Getting error of analysis '%s' (%d/%d): ", anaInfo$analysis[i], i, nrow(anaInfo))

        ind <- getDAFileIndex(DA, anaInfo$analysis[i], anaInfo$path[i])
        if (ind == -1)
        {
            warning(sprintf("Failed to open analysis '%s'", anaInfo$analysis[i]))
            cat("Failed!!\n")
        }
        else
        {
            # UNDONE: assume there is only one calibration item
            ret[i, "error"] <- DA[["Analyses"]][[ind]][["CalibrationStatus"]][[1]]$StandardDeviation()
            cat("Done!\n")
        }
    }

    return(ret)
}

makeDAEIC <- function(mz, mzWindow, ctype = "EIC", mtype = "MS", polarity = "both", bgsubtr = FALSE, fragpath = "")
{
    trace <- switch(ctype,
                    TIC = RDCOMClient::COMCreate("DataAnalysis.TICChromatogramDefinition"),
                    BPC = RDCOMClient::COMCreate("DataAnalysis.BPCChromatogramDefinition"),
                    EIC = RDCOMClient::COMCreate("DataAnalysis.EICChromatogramDefinition"),
                    stop("Wrong ctype: should be TIC, BPC OR EIC")
    )

    trace[["MSFilter"]][["Type"]] <- switch (mtype,
                                             all = DAConstants$daMSFilterAll,
                                             MS = DAConstants$daMSFilterMS,
                                             MSMS = DAConstants$daMSFilterMSMS,
                                             allMSMS = DAConstants$daMSFilterAllMSMS,
                                             BBCID = DAConstants$daMSFilterBBCID,
                                             stop("Wrong mtype: should be all, MS, MSMS, allMSMS or BBCID")
    )

    trace[["Polarity"]] <- switch (polarity,
                                   both = DAConstants$daBoth,
                                   positive = DAConstants$daPositive,
                                   negative = DAConstants$daNegative,
                                   stop("Wrong polarity: should be all, positive or negative")
    )

    trace[["BackgroundType"]] <- if (bgsubtr) DAConstants$daBgrdTypeSpectral else DAConstants$daBgrdTypeNone
    trace[["MSFilter"]][["FragmentationPath"]] <- fragpath

    if (ctype != "TIC")
    {
        trace[["Range"]] <- mz
        trace[["WidthLeft"]] <- mzWindow[1]
        trace[["WidthRight"]] <- mzWindow[1]
    }

    return(trace)
}

#' @details \code{addDAEIC} adds an Extracted Ion Chromatogram (EIC) or other
#'   chromatographic trace to a given analysis which can be used directly with
#'   DataAnalysis.
#'
#' @param analysis Analysis name (without file extension).
#' @param path path of the analysis.
#' @param mz \emph{m/z} (Da) value used for the chromatographic trace (if
#'   applicable).
#' @param mtype MS filter for chromatographic trace. Valid values are:
#'   \code{"all"}, \code{"MS"}, \code{"MSMS"}, \code{"allMSMS"} and
#'   \code{"BBCID"}.
#' @param polarity Polarity filter for chromatographic trace. Valid values:
#'   \code{"both"}, \code{"positive"} and \code{"negative"}.
#' @param fragpath Precursor \emph{m/z} used for MS/MS traces (\code{""} for
#'   none).

#' @rdname bruker-utils
#' @export
addDAEIC <- function(analysis, path, mz, mzWindow = 0.005, ctype = "EIC", mtype = "MS", polarity = "both", bgsubtr = FALSE, fragpath = "",
                     name = NULL, hideDA = TRUE, close = FALSE, save = close)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertString(analysis, min.chars = 1, add = ac)
    checkmate::assertString(path, add = ac)
    checkmate::assertDirectoryExists(file.path(path, paste0(analysis, ".d")), add = ac)
    checkmate::assertNumber(mz, lower = 0, finite = TRUE, add = ac)
    checkmate::assertNumber(mzWindow, lower = 0, finite = TRUE, add = ac)
    checkmate::assertChoice(ctype, c("EIC", "TIC", "BPC"), add = ac)
    checkmate::assertChoice(mtype, c("MS", "MSMS", "allMSMS", "BBCID"), add = ac)
    checkmate::assertChoice(polarity, c("positive", "negative", "both"), add = ac)
    checkmate::assertFlag(bgsubtr, add = ac)
    checkmate::assert(
        checkmate::checkNumber(fragpath, lower = 0, finite = TRUE),
        checkmate::checkChoice(fragpath, "")
    )
    checkmate::assertString(name, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertFlag(hideDA, add = ac)
    assertDACloseSaveArgs(close, save, add = ac)
    checkmate::reportAssertions(ac)

    DA <- getDAApplication()

    if (hideDA)
        hideDAInScope()

    ind <- getDAFileIndex(DA, analysis, path)
    ret <- NULL
    if (ind != -1)
    {
        chroms <- DA[["Analyses"]][[ind]][["Chromatograms"]]
        oldEICCount <- chroms$Count()
        chroms$AddChromatogram(makeDAEIC(mz, mzWindow, ctype, mtype, polarity, bgsubtr, fragpath))
        ret <- chroms$Count()
        if (!is.null(name) && oldEICCount < ret)
            chroms[[ret]][["Name_"]] <- name
        closeSaveDAFile(DA, ind, close, save)
    }

    return(ret)
}

# UNDONE: make method?

#' @details \code{addAllDAEICs} adds Extracted Ion Chromatograms (EICs) for all
#'   features within a \code{\link{featureGroups}} object.
#'
#' @param fGroups The \code{\link{featureGroups}} object for which EICs should be made.
#' @param onlyPresent If \code{TRUE} then EICs are only generated for analyses
#'   where the feature was detected.
#'
#' @rdname bruker-utils
#' @export
addAllDAEICs <- function(fGroups, mzWindow = 0.005, ctype = "EIC", bgsubtr = FALSE, name = TRUE,
                         onlyPresent = TRUE, hideDA = TRUE, close = FALSE, save = close)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertNumber(mzWindow, lower = 0, finite = TRUE, add = ac)
    checkmate::assertChoice(ctype, c("EIC", "BPC"), add = ac)
    aapply(checkmate::assertFlag, . ~ bgsubtr + name + hideDA, fixed = list(add = ac))
    assertDACloseSaveArgs(close, save, add = ac)
    checkmate::reportAssertions(ac)

    gInfo <- groupInfo(fGroups)
    gNames <- names(fGroups)
    gCount <- length(fGroups)
    ftInd <- groupFeatIndex(fGroups)
    anaInfo <- analysisInfo(fGroups)
    anaCount <- nrow(anaInfo)
    fTable <- featureTable(fGroups)

    if (gCount == 0)
        return(invisible(NULL))

    DA <- getDAApplication()

    if (hideDA)
        hideDAInScope()

    printf("Adding EICs for %d feature groups in DataAnalysis...\n", gCount)
    prog <- openProgBar(0, anaCount)

    for (anai in seq_len(anaCount))
    {
        grpsInAna <- seq_len(gCount)
        if (onlyPresent)
            grpsInAna <- grpsInAna[sapply(grpsInAna, function(i) ftInd[[i]][anai] != 0)]

        if (length(grpsInAna) > 0)
        {
            ind <- getDAFileIndex(DA, anaInfo$analysis[anai], anaInfo$path[anai])
            if (ind != -1)
            {
                eics <- sapply(grpsInAna, function(grpi)
                {
                    makeDAEIC(gInfo[grpi, "mzs"], mzWindow, ctype, bgsubtr = bgsubtr)
                }, USE.NAMES = FALSE)
                
                chroms <- DA[["Analyses"]][[ind]][["Chromatograms"]]
                oldEICCount <- chroms$Count()
                chroms$AddChromatograms(eics)
                newEICCount <- chroms$Count()
                
                if (!is.null(name) && name)
                {
                    if ((newEICCount - oldEICCount) != length(grpsInAna))
                        warning("Failed to add some EICs, cannot set names.")
                    else
                    {
                        gn <- gNames[grpsInAna]
                        for (eici in seq(oldEICCount+1, newEICCount))
                            chroms[[eici]][["Name_"]] <- gn[eici - oldEICCount]
                    }
                }
                
                closeSaveDAFile(DA, ind, close, save)
            }
        }

        setTxtProgressBar(prog, anai)
    }

    setTxtProgressBar(prog, anaCount)
    close(prog)

    invisible(NULL)
}

clearDAChromsAndSpecs <- function(DA, gNames, DAFind)
{
    # check if s matches our naming format
    isFGroupName <- function(s) any(!is.na(pmatch(gNames, s)))

    chroms <- DA[["Analyses"]][[DAFind]][["Chromatograms"]]
    ccount <- chroms$Count()
    if (ccount > 0)
    {
        for (i in ccount:1) # reverse loop: indices should stay valid after deletion
        {
            if (isFGroupName(chroms[[i]][["Name"]]))
                chroms$DeleteChromatogram(i)
        }
    }

    specs <- DA[["Analyses"]][[DAFind]][["Spectra"]]
    scount <- specs$Count()
    if (scount > 0)
    {
        for (i in scount:1)
        {
            if (isFGroupName(specs[[i]][["Name"]]))
                specs$Delete(i)
        }
    }
}

generateDAEICsForPeakLists <- function(DA, ana, path, bgsubtr, MSMSType, gNames, featInfo, DAFind)
{
    cat("Adding EICs for spectra generation... ")

    # add general TIC MS chromatogram used for generating MS spectra
    MSEIC <- addDAEIC(ana, path, 0, 0.005, "TIC", "MS", bgsubtr = bgsubtr, name = "MS TIC", hideDA = FALSE)
    stopifnot(!is.null(MSEIC))

    # add all MSMS traces
    eics <- sapply(gNames, function(g)
    {
        fmz <- featInfo[group == g, mz]
        makeDAEIC(fmz, 0.005, "TIC", MSMSType, bgsubtr = bgsubtr, fragpath = fmz)
    }, simplify = TRUE, USE.NAMES = FALSE) # NOTE: simplify/USE.NAMES have to be this way to not get strange DCOM errors.

    chroms <- DA[["Analyses"]][[DAFind]][["Chromatograms"]]
    oldEICCount <- chroms$Count()
    chroms$AddChromatograms(eics)
    newEICCount <- chroms$Count()

    MSMSEICs <- list()
    if (newEICCount > oldEICCount)
    {
        gCount <- length(gNames)

        # not all MSMS EICs may have been added when no MSMS data exists, find back which were added
        curgrpi <- 1
        for (eic in seq(oldEICCount + 1, newEICCount))
        {
            eicMz <- as.numeric(chroms[[eic]][["Definition"]][["MSFilter"]][["FragmentationPath"]])

            while (curgrpi <= gCount && !numEQ(featInfo[group == gNames[curgrpi], mz], eicMz, tol = 5e-3))
                curgrpi <- curgrpi + 1

            if (curgrpi > gCount)
                break

            MSMSEICs[[gNames[curgrpi]]] <- eic
            chroms[[eic]][["Name_"]] <- sprintf("%s - %s", gNames[curgrpi], MSMSType)
            curgrpi <- curgrpi + 1
        }
    }

    cat("Done!\n")

    return(list(MSEIC = MSEIC, MSMSEICs = MSMSEICs))
}

generateDASpecsForPeakLists <- function(DA, maxMSRtWindow, MSMSType, gNames, featInfo, DAEICs, DAFind)
{
    chroms <- DA[["Analyses"]][[DAFind]][["Chromatograms"]]
    specs <- DA[["Analyses"]][[DAFind]][["Spectra"]]

    addSpectrum <- function(eic, grp, rt, rtmin, rtmax, mz, mtype)
    {
        if (!is.null(maxMSRtWindow) && diff(c(rtmin, rtmax) * 2) > maxMSRtWindow*2)
        {
            rtmin <- max(rtmin, rt - maxMSRtWindow)
            rtmax <- min(rtmax, rt + maxMSRtWindow)
        }

        chroms[[eic]]$ClearRangeSelections()
        chroms[[eic]]$AddRangeSelection(rtmin/60, rtmax/60, 0, 0) # divide by 60: seconds to minutes

        oldSpecCount <- specs$Count()
        chroms[[eic]]$AverageMassSpectrum(1, 0)
        newSpecCount <- specs$Count()

        if (oldSpecCount != newSpecCount)
        {
            # HACK HACK HACK: changing the name of a spectrum throws an error but actually seems to work,
            # bug in RDCOM-client?
            tryCatch(specs[[newSpecCount]][["Name"]] <- sprintf("%s - %s", grp, mtype), error = function(e) e)
            return(newSpecCount)
        }

        return(NA)
    }

    gCount <- length(gNames)
    printf("Adding spectra for %d feature groups...\n", gCount)
    prog <- openProgBar(0, gCount)

    MSSpecs <- list(); MSMSSpecs <- list()
    for (grpi in seq_along(gNames))
    {
        grp <- gNames[grpi]
        fi <- featInfo[group == grp]

        spec <- addSpectrum(DAEICs$MSEIC, grp, fi$ret, fi$retmin, fi$retmax, fi$mz, "MS")
        if (is.na(spec))
            warning(sprintf("Failed to add MS spectrum for group %s, analysis %s, m/z %f", grp,
                            DA[["Analyses"]][[DAFind]][["Spectra"]]$Name, fi$mz))
        else
            MSSpecs[[grp]] <- spec

        if (!is.null(DAEICs$MSMSEICs[[grp]]))
        {
            spec <- addSpectrum(DAEICs$MSMSEICs[[grp]], grp, fi$ret, fi$retmin, fi$retmax, fi$mz, MSMSType)
            if (!is.na(spec))
                MSMSSpecs[[grp]] <- spec
        }

        setTxtProgressBar(prog, grpi)
    }

    setTxtProgressBar(prog, gCount)
    close(prog)

    cat("Deconvoluting spectra ...")
    specs$Deconvolute()
    cat("Done!\n")

    return(list(MSSpecs = MSSpecs, MSMSSpecs = MSMSSpecs))
}

getDAPeakList <- function(findDA, ind, useFMF, getMSMS, minInt)
{
    DA <- getDAApplication()

    if (useFMF)
        specDA <- DA[["Analyses"]][[findDA]][["Compounds"]][[ind]][[if (getMSMS) 2 else 1]]
    else
        specDA <- DA[["Analyses"]][[findDA]][["Spectra"]][[ind]]

    if (specDA[["MSPeakList"]]$Count() == 0)
        return(NULL)

    pfile <- tempfile("peaklist", fileext = ".csv")
    specDA$ExportMassList(pfile, DAConstants$daCSV)

    cnames <- names(fread(pfile, nrows = 0))
    cnames <- intersect(cnames, c("m/z", "I", "Cmpnt."))
    hasCmpnt <- "Cmpnt." %in% cnames
    newNames <- c("mz", "intensity", if (hasCmpnt) "cmp" else character())

    # force colclasses to prevent warnings, make sure iso compounds designated as 'NA' are still interpreted as a value
    plist <- fread(pfile, colClasses = if (hasCmpnt) c(`Cmpnt.` = "character") else NULL, na.strings = NULL,
                   select = cnames, col.names = newNames)

    unlink(pfile) # amount of peaklists may be large, remove temp files straight away

    return(plist[intensity >= minInt])
}

checkDAFMFCompounds <- function(DA, featTable, analysisInd, verify)
{
    cmpds <- DA[["Analyses"]][[analysisInd]][["Compounds"]]
    cCount <- cmpds[["Count"]]
    ret <- TRUE

    if (cCount < nrow(featTable))
        ret <- sprintf("Number of DataAnalysis compounds is less than number of features for analysis %s!",
                       DA[["Analyses"]][[analysisInd]][["Name"]])
    else
    {
        for (fti in seq_len(nrow(featTable)))
        {
            ci <- featTable$ID[fti]
            if (!grepl(paste0("^Cmpd ", ci, ", MolFeature"), cmpds[[ci]][["Name"]]))
            {
                ret <- sprintf("Could not find back the DataAnalysis compound for feature %d in analysis %s!",
                               fti, DA[["Analyses"]][[analysisInd]][["Name"]])
                break
            }
        }
    }

    if (verify && !isTRUE(ret))
        stop(paste(ret, "Please re-run FMF by calling findFeatures with doFMF=\"force\"."))

    return(isTRUE(ret))
}

# nocov end
