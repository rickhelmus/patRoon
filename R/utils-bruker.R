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
getBrukerAnalysisPath <- function(analysis, path) gsub("/", "\\\\", file.path(path, paste0(analysis, ".d")))

hideDAInScope <- local_(function(x) getDAApplication()$Hide(), function(x) getDAApplication()$Show())

#' @details \code{showDataAnalysis} makes a hidden DataAnalysis window visible
#'   again. Most functions using DataAnalysis will hide the window during
#'   processing for efficiency reasons. If the window remains hidden
#'   (\emph{e.g.} because there was an error) this function can be used to make
#'   it visible again.
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
            if (plist[[i]][["Intensity"]] > maxInt)
            {
                ret$mz <- plist[[i]][["m_over_z"]]
                ret$fwhm <- plist[[i]][["width"]]
                maxInt <- spec$Intensity(i)
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

getDAFileIndex <- function(DA, analysis, path, openFileIfClosed = TRUE)
{
    getFIndex <- function()
    {
        openCount <- DA[["Analyses"]][["Count"]]

        if (openCount > 0)
        {
            afile <- getBrukerAnalysisPath(analysis, path)
            for (i in 1:openCount)
            {
                if (DA[["Analyses"]][[i]][["Path"]] == afile)
                    return(i)
            }
        }

        return(-1)
    }

    ret <- getFIndex()
    if (ret == -1 && openFileIfClosed) # File is not yet open?
    {
        DA[["Analyses"]]$Open(getBrukerAnalysisPath(analysis, path))
        ret <- getFIndex()
    }

    if (ret == -1 && openFileIfClosed)
        warning(sprintf("Failed to open '%s'", getBrukerAnalysisPath(analysis, path)))

    return(ret)
}

closeDAFile <- function(DA, analysis, path, save)
{
    ind <- getDAFileIndex(DA, analysis, path, FALSE)
    if (ind != -1)
    {
        if (save)
            DA[["Analyses"]][[ind]]$Save()
        DA[["Analyses"]][[ind]]$Close()
    }
}

#' @details \code{setDAMethod} Sets a given DataAnalysis method (\file{.m} file)
#'   to a set of analyses.
#'
#' @param method The full path of the DataAnalysis method.
#'
#' @rdname bruker-utils
#' @export
setDAMethod <- function(anaInfo, method)
{
    DA <- getDAApplication()
    hideDAInScope()

    for (i in seq_len(nrow(anaInfo)))
    {
        printf("Setting DA method of analysis '%s' to %s (%d/%d)...\n", anaInfo$analysis[i], method, i, nrow(anaInfo))

        # HACK: need to close file (if open) otherwise method is not applied
        closeDAFile(DA, anaInfo$analysis[i], anaInfo$path[i], TRUE)

        ind <- getDAFileIndex(DA, anaInfo$analysis[i], anaInfo$path[i])
        if (ind == -1)
            next
        DA[["Analyses"]][[ind]]$LoadMethod(method)
        DA[["Analyses"]][[ind]]$Save()
    }

    invisible(NULL)
}

#' @details \code{recalibrarateDAFiles} Performs automatic mass recalibration of
#'   a given set of analyses. The current method settings for each analyses will
#'   be used.
#'
#' @rdname bruker-utils
#' @export
recalibrarateDAFiles <- function(anaInfo)
{
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
            DA[["Analyses"]][[ind]]$Save()
            # UNDONE: assume there is only one calibration item
            printf("Done! (std: %f ppm)\n", DA[["Analyses"]][[ind]][["CalibrationStatus"]][[1]]$StandardDeviation())
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

#' @details \code{exportDAFiles} will export a set of analyses either in
#'   \file{.mzXML} or \file{.mzML} formats.
#'
#' @param format The output format of exported files. Should be either
#'   \code{"mzXML"}, \code{"mzML"} or \code{"mzData"}.
#' @param exportLine Export line spectra (\code{TRUE}) or profile spectra
#'   (\code{FALSE}). Usually line spectra are preferred, since profile spectra
#'   use signficantly more disk space and increase required memory during
#'   processing.
#' @param outPath Character vector of output paths for exported analyses. Will
#'   be recycled if necessary.
#' @param overWrite If \code{TRUE} existing files will be overwritten.
#'
#' @rdname bruker-utils
#' @export
exportDAFiles <- function(anaInfo, format = "mzML", exportLine = TRUE, outPath = anaInfo$path, overWrite = FALSE)
{
    outPath <- rep(outPath, length.out = length(anaInfo$path))

    if (!format %in% c("mzXML", "mzData", "mzML"))
        stop("Wrong export format!")
    else if (length(anaInfo$path) != length(outPath))
        stop("Size of output paths differs from input")

    expConstant <- if (format == "mzXML") DAConstants$daMzXML else if (format == "mzData") DAConstants$daMzData else DAConstants$daMzML
    expSpecConstant <- if (exportLine) DAConstants$daLine else DAConstants$daProfile

    DA <- getDAApplication()
    hideDAInScope()

    for (i in seq_len(nrow(anaInfo)))
    {
        printf("Exporting analysis '%s' (%d/%d)... ", anaInfo$analysis[i], i, nrow(anaInfo))

        outf <- file.path(outPath[i], paste0(anaInfo$analysis[i], ".", format))
        if (!overWrite && file.exists(outf))
        {
            cat("Skipped: already exists.\n")
            next
        }

        ind <- getDAFileIndex(DA, anaInfo$analysis[i], anaInfo$path[i])
        if (ind == -1)
        {
            cat("Failed!!")
            next
        }

        DA[["Analyses"]][[ind]]$Export(outf, expConstant, expSpecConstant)

        cat("Done!\n")
    }

    invisible(NULL)
}

makeDAEIC <- function(mz, mzWidth, ctype = "EIC", mtype = "MS", polarity = "both", bgsubtr = FALSE, fragpath = "")
{
    trace <- switch (ctype,
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
        trace[["WidthLeft"]] <- mzWidth[1]
        trace[["WidthRight"]] <- if (length(mzWidth) > 1) mzWidth[2] else mzWidth[1]
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
#' @param mzWidth \emph{m/z} width (in Da) used for the chromatographic trace
#'   (if applicable).
#' @param ctype Type of the chromatographic trace. Valid options are:
#'   \code{"EIC"} (extracted ion chromatogram), \code{"TIC"} (total ion
#'   chromatogram) and \code{"BPC"} (Base Peak Chromatogram).
#' @param mtype MS filter for chromatographic trace. Valid values are:
#'   \code{"all"}, \code{"MS"}, \code{"MSMS"}, \code{"allMSMS"} and
#'   \code{"BBCID"}.
#' @param polarity Polarity filter for chromatographic trace. Valid values:
#'   \code{"both"}, \code{"positive"} and \code{"negative"}.
#' @param bgsubtr If \code{TRUE} then background subtraction ('Spectral'
#'   algorithm) will be performed.
#' @param fragpath Precursor \emph{m/z} used for MS/MS traces (\code{""} for
#'   none).
#' @param name The name for the chromatographic trace. \code{NULL} for none.
#' @param hideDA Hides DataAnalysis while adding the chromatographic trace
#'   (faster).
#'
#' @rdname bruker-utils
#' @export
addDAEIC <- function(analysis, path, mz, mzWidth, ctype = "EIC", mtype = "MS", polarity = "both", bgsubtr = FALSE, fragpath = "",
                     name = NULL, hideDA = TRUE)
{
    DA <- getDAApplication()

    if (hideDA)
        hideDAInScope()

    ind <- getDAFileIndex(DA, analysis, path)
    ret <- NULL
    if (ind != -1)
    {
        chroms <- DA[["Analyses"]][[ind]][["Chromatograms"]]
        oldEICCount <- chroms$Count()
        chroms$AddChromatogram(makeDAEIC(mz, mzWidth, ctype, mtype, polarity, bgsubtr, fragpath))
        ret <- chroms$Count()
        if (!is.null(name))
            chroms[[ret]][["Name_"]] <- name
    }

    return(ret)
}

generateDACompounds <- function(fGroups, bgsubtr, maxRtMSWidth, clear, save, MSMSType)
{
    ftindex <- groupFeatIndex(fGroups)
    fTable <- featureTable(fGroups)
    anaInfo <- analysisInfo(fGroups)
    DA <- getDAApplication()

    getPrintFGroupName <- function(g) sprintf("grp #%d %s", g, colnames(ftindex)[g])
    isFGroupName <- function(s) grepl("grp #([0-9])+", s)

    addSpectrum <- function(find, eic, grp, rt, rtmin, rtmax, mz, mtype)
    {
        if (!is.null(maxRtMSWidth) && diff(c(rtmin, rtmax)) > maxRtMSWidth)
        {
            rtmin <- max(rtmin, rt - maxRtMSWidth/2)
            rtmax <- min(rtmax, rt + maxRtMSWidth/2)
        }

        chroms <- DA[["Analyses"]][[find]][["Chromatograms"]]
        chroms[[eic]]$ClearRangeSelections()
        chroms[[eic]]$AddRangeSelection(rtmin/60, rtmax/60, 0, 0) # divide by 60: seconds to minutes

        specs <- DA[["Analyses"]][[find]][["Spectra"]]
        oldSpecCount <- specs$Count()
        chroms[[eic]]$AverageMassSpectrum(1, 0)
        newSpecCount <- specs$Count()

        if (oldSpecCount != newSpecCount)
        {
            # HACK HACK HACK: changing the name of a spectrum throws an error but actually seems to work,
            # bug in RDCOM-client?
            tryCatch(specs[[newSpecCount]][["Name"]] <- sprintf("%s (rt: %.0f, m/z %.4f) - %s", getPrintFGroupName(grp), rt, mz, mtype),
                     error = function(e) e)
            return(newSpecCount)
        }

        return(NA)
    }

    hideDAInScope()

    compounds <- lapply(seq_along(anaInfo$analysis), function(find)
    {
        ana <- anaInfo$analysis[find]
        path <- anaInfo$path[find]

        hash <- makeHash(fGroups, ana, bgsubtr, maxRtMSWidth, MSMSType)
        cmp <- loadCacheData("compoundsDA", hash)
        if (!is.null(cmp))
            return(cmp)

        findDA <- getDAFileIndex(DA, ana, path)
        if (findDA != -1)
        {
            cmp <- data.table(group = colnames(ftindex), MSEIC = numeric(ncol(ftindex)), MSMSEIC = numeric(ncol(ftindex)),
                              MSSpec = numeric(ncol(ftindex)), MSMSSpec = numeric(ncol(ftindex)))

            # generate table with feature info for each group within current analysis
            ftinfo <- rbindlist(lapply(seq_len(ncol(ftindex)), function(g)
            {
                fi <- ftindex[[g]][find]
                if (fi == 0)
                    data.table(ret = NA) # dummy NA table, needs at least one column to be merged
                else
                    fTable[[ana]][fi, .(ret, mz, retmin, retmax)]
            }), fill = TRUE) # fill to TRUE: merge NA tables

            gcount <- ftinfo[, sum(!is.na(ret))] # nr of feature groups to process
            if (gcount == 0)
                return(cmp)

            chroms <- DA[["Analyses"]][[findDA]][["Chromatograms"]]

            if (clear)
            {
                cat("Clearing old chromatograms/spectra... ")

                ccount <- chroms$Count()
                if (ccount > 0)
                {
                    for (i in ccount:1) # reverse loop: indices should stay valid after deletion
                    {
                        if (isFGroupName(chroms[[i]][["Name"]]))
                            chroms$DeleteChromatogram(i)
                    }
                }

                specs <- DA[["Analyses"]][[findDA]][["Spectra"]]
                scount <- specs$Count()
                if (scount > 0)
                {
                    for (i in scount:1)
                    {
                        if (isFGroupName(specs[[i]][["Name"]]))
                            specs$Delete(i)
                    }
                }

                cat("Done!\n")
            }

            cat("Adding EICs for spectra generation... ")
            # add general TIC MS chromatogram used for generating MS spectra
            MSEIC <- addDAEIC(ana, path, -1, 0.005, "TIC", "MS", bgsubtr = bgsubtr, name = "MS TIC", hideDA = FALSE)
            stopifnot(!is.null(MSEIC))
            cmp$MSEIC <- MSEIC

            # add all MSMS traces
            egrps <- seq_along(ftindex)
            egrps <- egrps[ftindex[find] > 0] # filter out groups which are not present in current analysis
            if (length(egrps) > 0)
            {
                eics <- sapply(egrps, function(g)
                {
                    makeDAEIC(ftinfo[g, mz], 0.005, "TIC", MSMSType, bgsubtr = bgsubtr, fragpath = ftinfo[g, mz])
                }, USE.NAMES = FALSE)

                oldEICCount <- chroms$Count()
                chroms$AddChromatograms(eics)
                newEICCount <- chroms$Count()

                if (newEICCount > oldEICCount)
                {
                    # not all MSMS EICs may have been added when no MSMS data exists, find back which were added
                    curegrpi <- 1
                    egrpslen <- length(egrps)
                    for (eic in (oldEICCount+1):newEICCount)
                    {
                        eicMz <- as.numeric(chroms[[eic]][["Definition"]][["MSFilter"]][["FragmentationPath"]])
                        while (curegrpi <= egrpslen && !isTRUE(all.equal(ftinfo[egrps[curegrpi], mz], eicMz, tolerance = 1e-5)))
                            curegrpi <- curegrpi + 1

                        if (curegrpi > egrpslen)
                            break

                        set(cmp, egrps[curegrpi], "MSMSEIC", eic)
                        chroms[[eic]][["Name_"]] <- sprintf("%s (rt: %.0f, m/z %.4f) - %s", getPrintFGroupName(egrps[curegrpi]),
                                                            ftinfo[egrps[curegrpi], ret], ftinfo[egrps[curegrpi], mz], MSMSType)
                        curegrpi <- curegrpi + 1
                    }
                }
            }
            cat("Done!\n")

            pgroups <- 0 # nr of groups processed so far
            printf("Adding spectra for %d feature groups to analysis '%s'...\n", gcount, ana)
            prog <- txtProgressBar(0, gcount, style=3)

            for (g in seq_along(ftindex))
            {
                if (ftindex[[g]][find] == 0) # feature not present
                    next

                spec <- addSpectrum(findDA, MSEIC, g, ftinfo[g, ret], ftinfo[g, retmin], ftinfo[g, retmax], ftinfo[g, mz], "MS")
                if (is.na(spec))
                    warning(sprintf("Failed to add MS spectrum for group %s, analysis %s, m/z %f", colnames(ftindex)[g], ana, mz))
                else
                    set(cmp, g, "MSSpec", spec)

                freic <- cmp$MSMSEIC[g]
                if (freic > 0)
                {
                    spec <- addSpectrum(findDA, freic, g, ftinfo[g, ret], ftinfo[g, retmin], ftinfo[g, retmax],
                                        ftinfo[g, mz], MSMSType)

                    if (!is.na(spec))
                        set(cmp, g, "MSMSSpec", spec)
                }

                pgroups <- pgroups + 1
                setTxtProgressBar(prog, pgroups)
            }

            setTxtProgressBar(prog, gcount)
            close(prog)

            saveCacheData("compoundsDA", cmp, hash)

            if (save)
                DA[["Analyses"]][[findDA]]$Save()

            return(cmp)
        }
    })

    names(compounds) <- anaInfo$analysis

    return(compounds)
}

getDAPeakList <- function(findDA, ind, useFMF, getMSMS)
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

    # force colclasses to prevent warnings, make sure iso compounds designated as 'NA' are still interpreted as a value
    plist <- fread(pfile, colClasses = c(`Cmpnt.` = "character"), na.strings = NULL)[, c("m/z", "I", "Cmpnt.")]
    setnames(plist, c("mz", "intensity", "cmp"))
    unlink(pfile) # amount of peaklists may be large, remove temp files straight away

    return(plist)
}
