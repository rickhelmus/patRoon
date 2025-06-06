# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
NULL

#' Report feature group data (legacy interface)
#'
#' Functionality to report data produced by most workflow steps such as features, feature groups, calculated chemical
#' formulae and tentatively identified compounds. This is the legacy interface, for the updated interface see
#' \link{reporting}.
#'
#' These functions are usually called at the very end of the workflow. It is used to report various data on features and
#' feature groups. In addition, these functions may be used for reporting formulae and/or compounds that were generated
#' for the specified feature groups. Data can be reported in tabular form (\emph{i.e.} \file{.csv} files) by
#' \code{reportCSV} or graphically by \code{reportPDF} and \code{reportHTML}. The latter functions will plot for
#' instance chromatograms and annotated mass spectra, which are useful to get a graphical overview of results.
#'
#' All functions have a wide variety of arguments that influence the reporting process. Nevertheless, most parameters
#' are optional and only required to be given for fine tuning. In addition, only those objects (\emph{e.g.} formulae,
#' compounds, clustering) that are desired to be reported need to be specified.
#'
#' @param fGroups The \code{\link{featureGroups}} object that should be used for reporting data.
#' @param path The destination file path for files generated during reporting. Will be generated if needed.
#' @param reportFGroups If \code{TRUE} then feature group data will be reported.
#' @param formulas,compounds,compsCluster,components Further objects (\code{\link{formulas}}, \code{\link{compounds}},
#'   \code{\link{compoundsCluster}}, \code{\link{components}}) that should be reported. Specify \code{NULL} to skip
#'   reporting a particular object.
#' @param reportFormulaSpectra If \code{TRUE} then explained MS/MS spectra (if available) for candidate formulae will be
#'   reported. Specifying \code{formulas} and setting this argument to \code{FALSE} still allows further annotation of
#'   compound MS/MS spectra.
#' @param MSPeakLists A \code{\link{MSPeakLists}} object that is \emph{mandatory} when spectra for formulae and/or
#'   compounds will be reported.
#' @param retMin If \code{TRUE} then report retention times in minutes (otherwise seconds).
#' @param compoundsOnlyUsedScorings If \code{TRUE} then only scorings are plotted that actually have been used to rank
#'   data (see the \code{scoreTypes} argument to \code{\link{generateCompoundsMetFrag}} for more details).
#' @param formulasTopMost,compoundsTopMost Only this amount of top ranked candidate formulae/compounds are reported.
#'   Lower values may significantly speed up reporting. Set to \code{NULL} to ignore.
#' @param clearPath If \code{TRUE} then the destination path will be (recursively) removed prior to reporting.
#'
#' @template EICParams-arg
#'
#' @templateVar normParam compoundsNormalizeScores,formulasNormalizeScores
#' @templateVar excludeParam compoundsExclNormScores,formulasExclNormScores
#' @template norm-args
#'
#' @note Any formulae and compounds for feature groups which are not present within \code{fGroups} (\emph{i.e.} because
#'   it has been subset afterwards) will not be reported.
#'
#'   The \code{topMost}, \code{topMostByRGroup} and \code{onlyPresent} \link[=EICParams]{EIC parameters} may be ignored,
#'   \emph{e.g.}, when generating overview plots.
#'
#' @seealso reporting
#'
#' @name reporting-legacy
NULL


prepareReportPath <- function(path, clear)
{
    if (clear)
        unlink(path, TRUE)
    mkdirp(path)
}

optimizePngPlots <- function(plotFiles)
{
    plotsPerBlock <- 50
    blocks <- ceiling(length(plotFiles) / plotsPerBlock)
    pqcmd <- getExtDepPath("pngquant")
    mainArgs <- c("--skip-if-larger", "--force")
    
    # older pngquant versions may not support the --strip argument yet
    if (any(grepl("--strip", suppressWarnings(executeCommand(pqcmd, stderr = TRUE)), fixed = TRUE)))
        mainArgs <- c(mainArgs, "--strip")
    
    cmdQueue <- lapply(seq_len(blocks), function(bl)
    {
        startpl <- (bl - 1) * plotsPerBlock + 1
        endpl <- min(startpl + plotsPerBlock - 1, length(plotFiles))
        return(list(block = bl, startpl = startpl, endpl = endpl, command = pqcmd,
                    args = c(mainArgs, plotFiles[seq(startpl, endpl)])))
    })

    message("Optimizing plot sizes...")
    executeMultiProcess(cmdQueue, function(cmd)
    {
        # pngquant create new files as <filename_without_extension>-fs8.png --> rename them to original names
        plots <- plotFiles[seq(cmd$startpl, cmd$endpl)]
        qplots <- paste0(tools::file_path_sans_ext(plots), "-fs8.png")
        qexist <- file.exists(qplots)
        file.rename(qplots[qexist], plots[qexist])
    }, method = "classic")

    invisible(NULL)
}

makeCachedPlot <- function(out, plotFunc, plotArgs, w, h, bg = "white", parSettings = NULL, cacheDB)
{
    hash <- makeHash(do.call(paste0(plotFunc, "Hash"), plotArgs), w, h, bg, parSettings)
    cache <- loadCacheData("reportPlots", hash, cacheDB)

    if (!is.null(cache))
        writeBin(cache, out)
    else
    {
        withr::with_png(out, width = w, height = h, units = "in", res = 72, bg = bg, code = {
            if (!is.null(parSettings))
                withr::with_par(parSettings, do.call(plotFunc, plotArgs))
            else
                do.call(plotFunc, plotArgs)
        })
        saveCacheData("reportPlots", readBin(out, "raw", file.info(out)$size), hash, cacheDB)
    }
}

reportFGroupTable <- function(fGroups, path, retMin)
{
    printf("Exporting feature group tables...")

    if (length(fGroups) == 0)
    {
        printf("No feature groups!\n")
        return(invisible(NULL))
    }

    # UNDONE: inheritance...
    tbl <- if (isScreening(fGroups)) as.data.table(fGroups, collapseSuspects = NULL) else as.data.table(fGroups)
    if (retMin)
        tbl[, ret := ret / 60]
    
    fwrite(tbl, file.path(path, sprintf("%s.csv", class(fGroups))))

    printf("Done!\n")
}

reportFGroupPlots <- function(fGroups, path, plotGrid, EICParams, retMin, EICs)
{
    printf("Exporting feature group plots...\n")

    if (length(fGroups) == 0)
    {
        printf("No feature groups!\n")
        return(invisible(NULL))
    }

    gCount <- length(fGroups)
    gNames <- names(fGroups)

    pdfFile <- file.path(path, sprintf("%s.pdf", class(fGroups)))
    withr::local_pdf(pdfFile, paper = "a4", pointsize = 10, width = 8, height = 11)

    # all feature groups
    plotChroms(fGroups, EICParams = getDefEICParams(rtWindow = EICParams$rtWindow, mzExpWindow = EICParams$mzExpWindow,
                                                    topMost = 1), retMin = retMin, EICs = EICs, showPeakArea = TRUE,
               showFGroupRect = FALSE)

    plotsPerPage <- plotGrid[1] * plotGrid[2]
    prog <- openProgBar(0, gCount)
    scr <- NULL
    for (grpi in seq_len(gCount))
    {
        scrInd <- grpi %% plotsPerPage
        if (scrInd == 0)
            scrInd <- plotsPerPage
        else if (scrInd == 1)
        {
            if (!is.null(scr))
                close.screen(scr)
            scr <- split.screen(plotGrid)
        }
        
        screen(scr[scrInd])
        plotChroms(fGroups, groupName = gNames[grpi], EICParams = EICParams, retMin = retMin, EICs = EICs,
                   colourBy = "rGroups")
        setTxtProgressBar(prog, grpi)
    }

    setTxtProgressBar(prog, gCount)
    close(prog)
    if (!is.null(scr))
        close.screen(scr)

    # UNDONE: may fail sometimes (on AppVeyor)
    # tools::compactPDF(pdfFile)
    # tools::compactPDF(pdfFile, gs_quality = "ebook")

    invisible(NULL)
}

reportFeatureTable <- function(fGroups, path, retMin)
{
    printf("Exporting feature tables...\n")

    fTable <- featureTable(fGroups)

    fcount <- length(fTable)
    prog <- openProgBar(0, fcount)

    for (anai in seq_along(fTable))
    {
        ana <- names(fTable)[anai]
        out <- file.path(path, sprintf("%s-%s.csv", class(getFeatures(fGroups)), ana))
        tab <- copy(fTable[[ana]])
        if (retMin)
            tab[, ret := ret / 60]
        write.csv(tab, out)

        setTxtProgressBar(prog, anai)
    }

    setTxtProgressBar(prog, fcount)
    close(prog)
}

reportFormulaTable <- function(fGroups, path, formulas, normalizeScores, excludeNormScores)
{
    printf("Exporting formula table...")

    gNames <- names(fGroups)

    for (grp in groupNames(formulas))
    {
        if (grp %in% gNames && nrow(formulas[[grp]]) > 0)
        {
            out <- file.path(path, sprintf("%s-%s.csv", class(fGroups), grp))

            ft <- formulas[[grp]][, -getAllMergedConsCols("fragInfo", names(formulas[[grp]]),
                                                          mergedConsensusNames(formulas)), with = FALSE]
            if (normalizeScores != "none")
                ft <- normalizeAnnScores(ft, annScoreNames(formulas, TRUE), formulas@scoreRanges[[grp]],
                                         mergedConsensusNames(formulas), normalizeScores == "minmax", excludeNormScores)
            write.csv(ft, out)
        }
    }

    printf("Done!\n")
}

reportFormulaSpectra <- function(fGroups, path, formulas, topMost, normalizeScores, excludeNormScores,
                                 MSPeakLists, EICParams, retMin, EICs)
{
    printf("Exporting formula MS/MS spectra...\n")

    if (length(formulas) == 0)
    {
        printf("No formulas!\n")
        return(invisible(NULL))
    }

    if (length(formulas) == 0)
        return(invisible(NULL))

    if (!is.null(topMost))
        formulas <- filter(formulas, topMost = topMost)

    formGroups <- intersect(groupNames(formulas), names(fGroups))
    fcount <- length(formGroups)
    if (fcount == 0)
        return(invisible(NULL))
    
    prog <- openProgBar(0, fcount)

    for (grp in formGroups)
    {
        ft <- formulas[[grp]]

        if (nrow(ft) > 0)
        {
            out <- file.path(path, sprintf("%s-%s.pdf", class(fGroups), grp))
            # a4r: width=11.69, height=8.27
            withr::with_pdf(out, paper = "a4r", pointsize = 10, width = 11, height = 8, code = {
                plotChroms(fGroups, groupName = grp, EICParams = EICParams, retMin = retMin, EICs = EICs)
                
                for (ind in seq_len(nrow(ft)))
                {
                    # NOTE: layout/mfrow/mfcol doesn't work because of the legend positioning (thinks 2 plots are made...)
                    
                    scr <- split.screen(c(2, 1))
                    scr <- c(scr, split.screen(c(1, 2), screen = scr[2]))
                    
                    screen(scr[1])
                    
                    if (is.null(MSPeakLists[[grp]][["MSMS"]]))
                    {
                        # no MSMS spectrum, i.e. MS only formula
                        textPlot("No MS/MS data available.")
                    }
                    else
                        plotSpectrum(formulas, ind, grp, MSPeakLists = MSPeakLists)
                    
                    screen(scr[3])
                    plotScores(formulas, ind, grp, normalizeScores = normalizeScores,
                               excludeNormScores = excludeNormScores)
                    
                    screen(scr[4])
                    
                    textPlot(paste0(getFormInfoList(formulas[[grp]], ind, mergedConsensusNames(formulas), FALSE),
                                    collapse = "\n"))
                    
                    close.screen(scr)
                }
            })
        }

        setTxtProgressBar(prog, match(grp, formGroups))
    }

    setTxtProgressBar(prog, fcount)
    close(prog)
}

reportCompoundTable <- function(fGroups, path, compounds, normalizeScores, excludeNormScores, compsCluster)
{
    printf("Exporting identification results tables...")

    gInfo <- groupInfo(fGroups)
    gNames <- names(fGroups)
    anaInfo <- analysisInfo(fGroups)
    compTable <- annotations(compounds)
    mcn <- mergedConsensusNames(compounds)

    if (!is.null(compsCluster))
        cutcl <- cutClusters(compsCluster)

    if (normalizeScores != "none")
        compTable <- Map(compTable, compounds@scoreRanges, f = normalizeAnnScores,
                         MoreArgs = list(scoreCols = annScoreNames(compounds, TRUE), mConsNames = mcn,
                                         minMaxNormalization = normalizeScores == "minmax",
                                         exclude = excludeNormScores))

    for (grp in names(compTable))
    {
        if (grp %in% gNames && nrow(compTable[[grp]]) > 0)
        {
            out <- file.path(path, sprintf("%s-%s.csv", class(fGroups), grp))
            tab <- compTable[[grp]][, -getAllMergedConsCols("fragInfo", names(compTable[[grp]]),
                                                      mergedConsensusNames(compounds)), with = FALSE]
            if (!is.null(compsCluster) && !is.null(cutcl[[grp]]))
                tab[, cluster := cutcl[[grp]]]
            write.csv(tab, out)
        }
    }

    printf("Done!\n")
}

reportCompoundSpectra <- function(fGroups, path, MSPeakLists, compounds, compsCluster, formulas, EICParams, retMin,
                                  EICs, normalizeScores, exclNormScores, onlyUsedScorings, topMost)
{
    printf("Exporting compound identification MS/MS spectra...\n")

    gNames <- names(fGroups)
    gInfo <- groupInfo(fGroups)
    anaInfo <- analysisInfo(fGroups)
    pLists <- peakLists(MSPeakLists)

    if (length(compounds) == 0)
    {
        printf("No compounds!\n")
        return(invisible(NULL))
    }

    if (!is.null(topMost))
        compounds <- filter(compounds, topMost = topMost)

    if (!is.null(compsCluster))
        cutcl <- cutClusters(compsCluster)

    compTable <- annotations(compounds)
    idcount <- length(compTable)
    prog <- openProgBar(0, idcount)

    for (gi in seq_along(compTable))
    {
        grp <- names(compTable)[gi]
        if (nrow(compTable[[grp]]) > 0)
        {
            fgrpi <- match(grp, gNames)
            if (is.na(fgrpi))
                next

            out <- file.path(path, sprintf("%s-%s.pdf", class(fGroups), grp))
            # a4r: width=11.69, height=8.27
            withr::with_pdf(out, paper = "a4r", pointsize = 10, width = 11, height = 8, code = {
                plotChroms(fGroups, groupName = grp, EICParams = EICParams, retMin = retMin, EICs = EICs)
                
                for (idi in seq_len(nrow(compTable[[grp]])))
                {
                    # NOTE: layout/mfrow/mfcol doesn't work because of the legend positioning (thinks 2 plots are made...)
                    
                    scr <- split.screen(c(2, 1))
                    scr <- c(scr, split.screen(c(1, 2), screen = scr[2]))
                    
                    screen(scr[1])
                    plotSpectrum(compounds, idi, grp, MSPeakLists = MSPeakLists, formulas = formulas,
                                 plotStruct = TRUE)
                    
                    screen(scr[3])
                    plotScores(compounds, idi, grp, normalizeScores, exclNormScores, onlyUsedScorings)
                    
                    screen(scr[4])
                    
                    # draw text info
                    txt <- paste0(getCompInfoList(compTable[[grp]], idi, mergedConsensusNames(compounds), FALSE),
                                  collapse = "\n")
                    if (!is.null(compsCluster) && !is.null(cutcl[[grp]]))
                        txt <- paste(txt, sprintf("cluster: %d", cutcl[[grp]][idi]), sep = "\n")
                    
                    textPlot(txt)
                    
                    close.screen(scr)
                }
            })
        }
        setTxtProgressBar(prog, gi)
    }

    setTxtProgressBar(prog, idcount)
    close(prog)
}

reportCompoundClusters <- function(fGroups, compsCluster, path)
{
    printf("Exporting compound clusters...\n")

    cutcl <- cutClusters(compsCluster)
    cutcl <- cutcl[names(cutcl) %in% names(fGroups)]
    ls <- lengths(compsCluster)
    ccount <- length(cutcl)

    if (ccount == 0)
    {
        printf("No compound clusters!\n")
        return(invisible(NULL))
    }

    prog <- openProgBar(0, ccount)

    for (gi in seq_along(cutcl))
    {
        grp <- names(cutcl)[gi]
        if (length(cutcl[[grp]]) > 0)
        {
            out <- file.path(path, sprintf("%s-%s-clusters.pdf", class(fGroups), grp))
            withr::with_pdf(out, paper = "a4", code =
            {
                plot(compsCluster, groupName = grp)
                par(mfrow = c(3, 3))
                for (i in seq_len(ls[[grp]]))
                    plotStructure(compsCluster, groupName = grp, cluster = i)
            })
        }
        setTxtProgressBar(prog, gi)
    }

    setTxtProgressBar(prog, ccount)
    close(prog)
}

reportComponentTable <- function(components, path, retMin)
{
    printf("Exporting component table...")

    if (length(components) == 0)
        printf("No components!\n")
    else
    {
        cTable <- rbindlist(componentTable(components), idcol = "component")

        if (retMin)
            cTable[, ret := ret / 60]
        write.csv(cTable, file.path(path, "components.csv"))

        printf("Done!\n")
    }
}

reportComponentPlots <- function(fGroups, path, components, EICParams, retMin, EICs)
{
    printf("Exporting component plots...\n")

    if (length(components) == 0)
    {
        printf("No components!\n")
        return(invisible(NULL))
    }

    gInfo <- groupInfo(fGroups)
    anaInfo <- analysisInfo(fGroups)
    cInfo <- componentInfo(components)
    cTable <- componentTable(components)

    prog <- openProgBar(0, length(components))

    withr::local_pdf(file.path(path, "components.pdf"), paper = "a4", pointsize = 10, width = 8, height = 11)

    # HACK: this should be replaced some proper inheritance/methods at some point
    isHClust <- inherits(components, "componentsIntClust")

    if (isHClust)
    {
        plotHeatMap(components, interactive = FALSE)
        plot(components)
        clProps <- clusterProperties(components)
    }

    for (cmpi in seq_len(length(components)))
    {
        if (isHClust)
        {
            scr <- split.screen(c(3, 1))
            scr <- c(scr, split.screen(c(1, 2), scr[3]))
        }
        else
            scr <- split.screen(c(2, 1))

        screen(scr[1])
        plotChroms(components, cmpi, fGroups, title = sprintf("Component %d", cmpi),
                   EICParams = getDefEICParams(rtWindow = EICParams$rtWindow, mzExpWindow = EICParams$mzExpWindow),
                   retMin = retMin, EICs = EICs)

        screen(scr[2])
        plotSpectrum(components, cmpi, main = sprintf("ret: %.1f; m/z: %.4f - %.4f",
                                                      cInfo$ret[cmpi], min(cTable[[cmpi]]$mz),
                                                      max(cTable[[cmpi]]$mz)))

        if (isHClust)
        {
            screen(scr[4])
            plotInt(components, index = cmpi, plotArgs = list(main = "normalized"))

            screen(scr[5])
            fg <- fGroups[, unique(cTable[[cmpi]]$group)]
            plotInt(fg, average = clProps$average, plotArgs = list(main = "absolute"))
        }

        close.screen(scr)

        setTxtProgressBar(prog, cmpi)
    }

    close(prog)
}

#' @details \code{reportCSV} generates tabular data (\emph{i.e.} \file{.csv}
#'   files) for given data to be reported. This may also be useful to allow
#'   import by other tools for post processing.
#'
#' @param reportFeatures If set to \code{TRUE} then for each analysis a
#'   \file{.csv} file will be generated with information about its detected
#'   features.
#'
#' @rdname reporting-legacy
#' @aliases reportCSV
#' @export
setMethod("reportCSV", "featureGroups", function(fGroups, path, reportFeatures, formulas,
                                                 formulasNormalizeScores, formulasExclNormScores,
                                                 compounds, compoundsNormalizeScores, compoundsExclNormScores,
                                                 compsCluster, components, retMin, clearPath)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertPathForOutput(path, overwrite = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ reportFeatures + retMin + clearPath, fixed = list(add = ac))
    aapply(checkmate::assertClass, . ~ formulas + compounds + compsCluster + components,
           c("formulas", "compounds", "compoundsCluster", "components"),
           null.ok = TRUE, fixed = list(add = ac))
    aapply(assertNormalizationMethod, . ~ formulasNormalizeScores + compoundsNormalizeScores, fixed = list(add = ac))
    aapply(checkmate::assertCharacter, . ~ formulasExclNormScores + compoundsExclNormScores,
           min.chars = 1, null.ok = TRUE, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    if (length(fGroups) == 0)
    {
        cat("No feature groups, nothing to report...\n")
        return(invisible(NULL))
    }
    
    prepareReportPath(path, clearPath)

    reportFGroupTable(fGroups, path, retMin)

    if (reportFeatures)
    {
        p <- file.path(path, "features")
        mkdirp(p)
        reportFeatureTable(fGroups, p, retMin)
    }

    if (!is.null(formulas) && length(formulas) > 0)
    {
        p <- file.path(path, "formulas")
        mkdirp(p)
        reportFormulaTable(fGroups, p, formulas, formulasNormalizeScores, formulasExclNormScores)
    }

    if (!is.null(compounds))
    {
        p <- file.path(path, "compounds")
        mkdirp(p)
        reportCompoundTable(fGroups, p, compounds, compoundsNormalizeScores, compoundsExclNormScores, compsCluster)
    }

    if (!is.null(components))
        reportComponentTable(components, path, retMin)
})

#' @details \code{reportPDF} will report graphical data (\emph{e.g.} chromatograms and mass spectra) within PDF files.
#'   Compared to \code{reportHTML} this function may be faster and yield smaller report files, however, its
#'   functionality is a bit more basic and generated data is more 'scattered' around.
#'
#' @param EICGrid An integer vector in the form \code{c(columns, rows)} that is used to determine the plotting grid when
#'   reporting EICs in PDF files.
#'
#' @rdname reporting-legacy
#' @aliases reportPDF
#' @export
setMethod("reportPDF", "featureGroups", function(fGroups, path, reportFGroups,
                                                 formulas, formulasTopMost, formulasNormalizeScores,
                                                 formulasExclNormScores, reportFormulaSpectra,
                                                 compounds, compoundsNormalizeScores, compoundsExclNormScores,
                                                 compoundsOnlyUsedScorings, compoundsTopMost, compsCluster,
                                                 components, MSPeakLists, retMin, EICGrid, EICParams = getDefEICParams(),
                                                 clearPath)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertPathForOutput(path, overwrite = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ reportFGroups + reportFormulaSpectra + compoundsOnlyUsedScorings + retMin +
               clearPath, fixed = list(add = ac))
    aapply(checkmate::assertClass, . ~ formulas + compounds + compsCluster + components + MSPeakLists,
           c("formulas", "compounds", "compoundsCluster", "components", "MSPeakLists"),
           null.ok = TRUE, fixed = list(add = ac))
    aapply(assertNormalizationMethod, . ~ formulasNormalizeScores + compoundsNormalizeScores, fixed = list(add = ac))
    aapply(checkmate::assertCharacter, . ~ formulasExclNormScores + compoundsExclNormScores,
           min.chars = 1, null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertCount, . ~ formulasTopMost + compoundsTopMost, positive = TRUE, null.ok = TRUE,
           fixed = list(add = ac))
    checkmate::assertIntegerish(EICGrid, lower = 1, any.missing = FALSE, len = 2, add = ac)
    assertEICParams(EICParams, add = ac)
    checkmate::reportAssertions(ac)

    if ((!reportFGroups && is.null(formulas) && is.null(compounds) && is.null(components)) || length(fGroups) == 0)
    {
        cat("No feature groups, nothing to report...\n")
        return(invisible(NULL))
    }

    if (is.null(MSPeakLists) &&
        ((!is.null(formulas) && reportFormulaSpectra) || !is.null(compounds)))
        stop("MSPeakLists is NULL, please specify when reporting formula and/or compounds")

    prepareReportPath(path, clearPath)

    if (reportFGroups || !is.null(formulas) || !is.null(compounds) || !is.null(components))
    {
        cat("Loading all EICs... ")
        EICs <- getEICsForFGroups(fGroups, EICParams = EICParams)
        cat("Done!\n")
    }

    if (reportFGroups)
        reportFGroupPlots(fGroups, path, EICGrid, EICParams, retMin, EICs)

    if (reportFormulaSpectra && !is.null(formulas))
    {
        p <- file.path(path, "formulas")
        mkdirp(p)
        reportFormulaSpectra(fGroups, p, formulas, formulasTopMost, formulasNormalizeScores,
                             formulasExclNormScores, MSPeakLists, EICParams, retMin, EICs)
    }

    if (!is.null(compounds))
    {
        p <- file.path(path, "compounds")
        mkdirp(p)
        reportCompoundSpectra(fGroups, p, MSPeakLists, compounds, compsCluster, formulas, EICParams, retMin, EICs,
                              compoundsNormalizeScores, compoundsExclNormScores, compoundsOnlyUsedScorings,
                              compoundsTopMost)
    }

    if (!is.null(compsCluster))
    {
        p <- file.path(path, "compounds")
        mkdirp(p)
        reportCompoundClusters(fGroups, compsCluster, p)
    }

    if (!is.null(components))
        reportComponentPlots(fGroups, path, components, EICParams, retMin, EICs)
})


#' @details \code{reportHTML} will report graphical data (\emph{e.g.} chromatograms and mass spectra) and summary
#'   information in an easy browsable \code{HTML} file using \link{rmarkdown}, \link{flexdashboard} and \link{knitr}.
#'
#' @param reportPlots A character vector specifying what should be plotted. Valid options are: \code{"chord"},
#'   \code{"venn"}, \code{"upset"} (plot a chord, Venn and UpSet diagram, respectively), \code{"eics"} (plot EICs for
#'   individual feature groups) and \code{"formulas"} (plot annotated formula spectra). Set to \code{"none"} to plot
#'   none of these.
#' @param includeMFWebLinks A \code{character} specifying to which feature groups a web link should be added in the
#'   annotation page to \href{https://msbi.ipb-halle.de/MetFragBeta/index.xhtml}{MetFragWeb}. Options are:
#'   \code{"compounds"} (only to those with compounds results), \code{"MSMS"} (only to those with MSMS peak lists) or
#'   \code{"none"}.
#' @param interactiveHeat If \code{TRUE} an interactive heatmap HTML widget will be generated to display hierarchical
#'   clustering results. Set to \code{FALSE} for a 'regular' static plot.
#' @param TPs A \code{\link{transformationProducts}} object that should used for plotting hierarchies. Ignored if
#'   \code{TPs=NULL} or \code{components} is not a \code{componentsTPs} object.
#' @param TPGraphStructuresMax Maximum number of TP structures to plot in TP hierarchies, see the \code{TPs} argument.
#'   Sets the \code{structuresMax} argument to \code{\link[=plotGraph,featureGroups-method]{plotGraph}}.
#' @param selfContained If \code{TRUE} the output will be a standalone HTML file which contains all graphics and script
#'   dependencies. When \code{FALSE}, the latter will be placed in an additional directory (\file{report_files}) which
#'   should remain present when viewing the output file. Especially on Windows, a non-self contained output might be
#'   desirable when reporting large amounts of data to prevent \command{pandoc} from running out of memory.
#' @param optimizePng If \code{TRUE} then \command{pngquant} is used to reduce the size of generated graphics. A
#'   significant reduction in disk space usage may be seen, however, at the cost additional processing time. Multiple
#'   \command{pngquant} processes will be executed in parallel, which can be configured with
#'   \option{patRoon.MP.maxProcs} (parallelization will always happen with the \code{"classic"} method, see
#'   \link[=patRoon-package]{patRoon options}).
#' @param openReport If set to \code{TRUE} then the output report file will be opened with the system browser.
#' @param noDate If \code{TRUE} then the current date is not added to the report. This is mainly used for testing and
#'   its main purpose is to guarantees equal report files when \code{reportHTML()} is called multiple times with equal
#'   arguments.
#'
#' @template specSimParams-arg
#'
#' @templateVar what \code{reportHTML}
#' @template uses-multiProc
#'
#' @section Parallelization: Currently, \code{reportHTML} only uses \code{"classic"} multiprocessing, regardless of the
#'   \option{patRoon.MP.method} option.
#'
#' @references Creating MetFrag landing page URLs based on code from
#'   \href{https://github.com/Treutler/MetFamily}{MetFamily} R package. \cr\cr \addCitations{knitr}
#'
#' @rdname reporting-legacy
#' @aliases reportHTML
#' @export
setMethod("reportHTML", "featureGroups", function(fGroups, path, reportPlots, formulas, formulasTopMost,
                                                  formulasNormalizeScores, formulasExclNormScores,
                                                  compounds, compoundsNormalizeScores, compoundsExclNormScores,
                                                  compoundsOnlyUsedScorings, compoundsTopMost,
                                                  compsCluster, includeMFWebLinks, components, interactiveHeat,
                                                  MSPeakLists, specSimParams, TPs, retMin, EICParams,
                                                  TPGraphStructuresMax, selfContained, optimizePng, clearPath,
                                                  openReport, noDate)
{
    # UNDONE: mention pandoc win limits

    ac <- checkmate::makeAssertCollection()
    checkmate::assertPathForOutput(path, overwrite = TRUE, add = ac)
    reportPlots <- checkmate::matchArg(reportPlots, c("none", "chord", "venn", "upset", "eics", "formulas"),
                                       several.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ compoundsOnlyUsedScorings + interactiveHeat + retMin + selfContained +
               optimizePng + clearPath + openReport + noDate, fixed = list(add = ac))
    aapply(checkmate::assertClass, . ~ formulas + compounds + components + MSPeakLists + TPs,
           c("formulas", "compounds", "components", "MSPeakLists", "transformationProducts"),
           null.ok = TRUE, fixed = list(add = ac))
    aapply(assertNormalizationMethod, . ~ formulasNormalizeScores + compoundsNormalizeScores, fixed = list(add = ac))
    aapply(checkmate::assertCharacter, . ~ formulasExclNormScores + compoundsExclNormScores,
           min.chars = 1, null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertCount, . ~ formulasTopMost + compoundsTopMost, positive = TRUE, null.ok = TRUE,
           fixed = list(add = ac))
    assertEICParams(EICParams, add = ac)
    checkmate::assertCount(TPGraphStructuresMax, add = ac)
    checkmate::assertChoice(includeMFWebLinks, c("compounds", "MSMS", "none"), add = ac)
    assertSpecSimParams(specSimParams, add = ac)
    checkmate::reportAssertions(ac)

    if (length(fGroups) == 0)
    {
        cat("No feature groups, nothing to report...\n")
        return(invisible(NULL))
    }
    
    if ("none" %in% reportPlots)
        reportPlots <- ""

    if (is.null(MSPeakLists) &&
        ((!is.null(formulas) && "formulas" %in% reportPlots) || !is.null(compounds)))
        stop("MSPeakLists is NULL, please specify when reporting formula and/or compounds")

    prepareReportPath(path, clearPath)

    workPath <- file.path(tempdir(TRUE), "report")
    unlink(workPath, TRUE)
    mkdirp(workPath)

    file.copy(system.file("report-legacy", "main.Rmd", package = "patRoon"), workPath)
    file.copy(system.file("report-legacy", "featinfo.Rmd", package = "patRoon"), workPath)
    file.copy(system.file("report-legacy", "annotation.Rmd", package = "patRoon"), workPath)
    file.copy(system.file("report-legacy", "components.Rmd", package = "patRoon"), workPath)
    file.copy(system.file("report-legacy", "TPs.Rmd", package = "patRoon"), workPath)

    # rmarkdown needs absolute path as relative paths will be from the path of the Rmd
    if (!R.utils::isAbsolutePath(path))
        path <- R.utils::getAbsolutePath(path)

    EICs <- NULL

    # NOTE: loading EICs is done outside knitr environment: the current working
    # directory is changed so we loose access to the cache.

    # UNDONE: always load EICs? Need them for summary EIC plot
    # if (reportFGroups || (!is.null(formulas) && reportFormulaSpectra) ||
    #     !is.null(compounds) || !is.null(components))
    {
        cat("Loading all EICs... ")
        EICs <- getEICsForFGroups(fGroups, EICParams = EICParams)
        cat("Done!\n")
    }

    rmdVars <- list(outPath = path, fGroups = fGroups, groupNames = names(fGroups), gInfo = groupInfo(fGroups),
                    reportPlots = reportPlots, EICParams = EICParams, retMin = retMin, EICs = EICs,
                    compounds = compounds, compoundsTopMost = compoundsTopMost, compsCluster = compsCluster,
                    includeMFWebLinks = includeMFWebLinks, MSPeakLists = MSPeakLists, formulas = formulas,
                    formulasTopMost = formulasTopMost, formulasNormalizeScores = formulasNormalizeScores,
                    formulasExclNormScores = formulasExclNormScores,
                    compoundsNormalizeScores = compoundsNormalizeScores,
                    compoundsExclNormScores = compoundsExclNormScores,
                    compoundsOnlyUsedScorings = compoundsOnlyUsedScorings,
                    components = components, interactiveHeat = interactiveHeat, specSimParams = specSimParams,
                    TPs = TPs, TPGraphStructuresMax = TPGraphStructuresMax, selfContained = selfContained,
                    optimizePng = optimizePng, noDate = noDate)

    # HACK: not sure what exactly happens here, but... kableExtra adds latex
    # dependencies by default, which then may cause serious memory leakage when
    # rmarkdown::render() is called repeatedly. For now just remove them temporarily.
    knitMeta <- knitr::knit_meta("latex_dependency", clean = TRUE)
    on.exit(knitr::knit_meta_add(knitMeta), add = TRUE)

    outputFile <- file.path(path, "report.html")

    # normalize cache path so it can be used in report working directory
    withr::with_options(list(DT.warn.size = FALSE,
                             patRoon.cache.fileName = normalizePath(getOption("patRoon.cache.fileName")),
                             patRoon.progress.opts = list(file = stderr())),
                        rmarkdown::render(file.path(workPath, "main.Rmd"), output_file = outputFile,
                                          output_options = list(self_contained = selfContained),
                                          quiet = TRUE))

    if (openReport)
        utils::browseURL(paste0("file://", normalizePath(outputFile)))

    invisible(NULL)
})
