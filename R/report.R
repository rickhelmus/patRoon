#' @include main.R
NULL

#' Report feature group data
#'
#' Functionality to report data produced by most workflow steps such as
#' features, feature groups, calculated chemical formulae and tentatively
#' identified compounds.
#'
#' These functions are usually called at the very end of the workflow. It is
#' used to report various data on features and feature groups. In addition,
#' these functions may be used for reporting formulae and/or compounds that were
#' generated for the specified feature groups. Data can be reported in tabular
#' form (\emph{i.e.} \file{.csv} files) by \code{reportCSV} or graphically by
#' \code{reportPDF} and \code{reportMD}. The latter functions will plot for
#' instance chromatograms and annotated mass spectra, which are useful to get a
#' graphical overview of results.
#'
#' All functions have a wide variety of arguments that influence the reporting
#' process. Nevertheless, most parameters are optional and only required to be
#' given for fine tuning. In addition, only those objects (\emph{e.g.} formulae,
#' compounds, clustering) that are desired to be reported need to be specified.
#'
#' @param fGroups The \code{\link{featureGroups}} object that should be used for
#'   reporting data.
#' @param path The destination file path files generated during reporting.
#' @param reportFGroups If \code{TRUE} then feature group data will be reported.
#' @param formConsensus,compounds,components Further objects
#'   (\code{\link{formulaConsensus}}, \code{\link{compounds}},
#'   \code{\link{components}}) that should be reported. Specify \code{NULL} to
#'   skip reporting a particular object.
#' @param reportFormulaSpectra If \code{TRUE} then explained MS/MS spectra (if
#'   available) for candidate formulae will be reported. Specifying
#'   \code{formConsensus} and setting this argument to \code{FALSE} still allows
#'   further annotation of compound MS/MS spectra.
#' @param compoundNormalizeScores If \code{TRUE} compound identification scores
#'   will be normalized to maximum values.
#' @param cInfo If not \code{NULL} a \code{\link{clusterInfo}} object to be used
#'   for reporting hierarchical clustering results.
#' @param clusterK Number of branches the hierarchical cluster should be cut
#'   into. If specified (\emph{i.e.} not \code{NULL}) then individual clusters
#'   will be reported.
#' @param silInfo A \code{\link{silhouetteInfo}} object to be used for reporting
#'   a silhouette plot.
#' @param clusterMaxLabels Maximum amount of clusters when labels are still
#'   drawn within a dendrogram.
#' @param MSPeakLists A \code{\link{MSPeakLists}} object that is
#'   \emph{mandatory} when spectra for formulae and/or compounds will be
#'   reported.
#' @param retMin If \code{TRUE} then report retention times in minutes
#'   (otehrwise seconds).
#' @param EICRtWindow,EICMzWindow,EICTopMost,EICOnlyPresent Plotting parameters
#'   passed to \code{\link{plotEIC}} (\emph{i.e.} \code{rtWindow},
#'   \code{mzWindow}, \code{topMost} and \code{onlyPresent} arguments).
#' @param clearPath If \code{TRUE} then the destination path will be
#'   (recursively) removed prior to reporting.
#'
#' @note Any formulae and compounds for feature groups which are not present
#'   within \code{fGroups} (\emph{i.e.} because it has been subset afterwards)
#'   will not be reported.
#'
#' @name reporting
NULL


prepareReportPath <- function(path, clear)
{
    if (clear)
        unlink(path, TRUE)
    mkdirp(path)
}

optimizePngPlots <- function(plotFiles, progressOut, maxProcAmount)
{
    plotsPerBlock <- 50
    blocks <- ceiling(length(plotFiles) / plotsPerBlock)
    pqcmd <- getCommandWithOptPath("pngquant", "pngquant")
    mainArgs <- c("--skip-if-larger", "--strip", "--force")
    cmdQueue <- lapply(seq_len(blocks), function(bl)
    {
        startpl <- (bl - 1) * plotsPerBlock + 1
        endpl <- startpl + plotsPerBlock - 1
        return(list(block = bl, startpl = startpl, endpl = endpl, command = pqcmd,
                    args = c(mainArgs, plotFiles[seq(startpl, endpl)])))
    })

    message("Optimizing plot sizes...")
    executeMultiProcess(cmdQueue, function(cmd, ...)
    {
        # pngquant create new files as <filename_without_extension>-fs8.png --> rename them to original names
        plots <- plotFiles[seq(cmd$startpl, cmd$endpl)]
        qplots <- paste0(tools::file_path_sans_ext(plots), "-fs8.png")
        qexist <- file.exists(qplots)
        file.rename(qplots[qexist], plots[qexist])
    }, progressOut = progressOut, maxProcAmount = maxProcAmount)

    invisible(NULL)
}

reportFGroupTable <- function(fGroups, path, fGroupsAsRows, reportAnalysisInfo, reportRetMz, retMin)
{
    printf("Exporting feature group table...")

    gTable <- copy(groups(fGroups))
    gInfo <- groupInfo(fGroups)
    anaInfo <- analysisInfo(fGroups)

    if (reportAnalysisInfo)
    {
        gTable <- insertDTColumn(gTable, "ref", anaInfo$ref, 1)
        gTable <- insertDTColumn(gTable, "replicate_group", anaInfo$group, 1)
    }

    if (fGroupsAsRows)
    {
        gTable <- transpose(gTable)
        setnames(gTable, anaInfo$analysis)
        initbl <- if (reportAnalysisInfo) rep("", 2) else c()

        rnames <- rownames(gInfo)
        if (reportAnalysisInfo)
            rnames <- c("replicate_group", "ref", rnames)

        if (reportRetMz)
        {
            gTable <- insertDTColumn(gTable, "m/z", c(initbl, gInfo$mzs), 1)
            rts <- if (retMin) gInfo$rts / 60 else gInfo$rts
            gTable <- insertDTColumn(gTable, "retention", c(initbl, rts), 1)
        }
    }
    else
    {
        if (reportRetMz)
        {
            gi <- setnames(transpose(gInfo[, c("rts", "mzs")]), rownames(gInfo))
            if (reportAnalysisInfo) # add two dummy columns
                gi <- cbind(replicate_group = "", ref = "", gi)

            gTable <- rbind(gi, gTable)
            rnames <- c("retention", "mz", anaInfo$analysis)
        }
    }

    write.csv(gTable, file.path(path, sprintf("%s.csv", class(fGroups))), row.names = rnames)

    printf("Done!\n")
}

reportFGroupPlots <- function(fGroups, path, plotGrid, rtWindow, mzWindow, retMin, topMost, EICs)
{
    printf("Exporting feature group plots...\n")

    gTable <- groups(fGroups)
    gCount <- length(fGroups)
    gInfo <- groupInfo(fGroups)
    anaInfo <- analysisInfo(fGroups)

    pdfFile <- file.path(path, sprintf("%s.pdf", class(fGroups)))
    pdf(pdfFile, paper = "a4", pointsize = 10, width = 8, height = 11)

    # all feature groups
    plotEIC(fGroups, rtWindow, mzWindow, retMin, 1, EICs, TRUE, FALSE)

    par(mfrow = plotGrid)
    prog <- txtProgressBar(0, gCount, style=3)
    for (grpi in seq_len(gCount))
    {
        plotEIC(fGroups[, grpi], rtWindow, mzWindow, retMin, topMost, EICs)
        setTxtProgressBar(prog, grpi)
    }

    setTxtProgressBar(prog, gCount)
    close(prog)
    dev.off()

    tools::compactPDF(pdfFile)
    # tools::compactPDF(pdfFile, gs_quality = "ebook")

    invisible(NULL)
}

reportFeatureTable <- function(fGroups, path, retMin)
{
    printf("Exporting feature tables...\n")

    fTable <- featureTable(fGroups)

    fcount <- length(fTable)
    prog <- txtProgressBar(0, fcount, style = 3)

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

reportFormulaTable <- function(fGroups, formConsensus, path, retMin)
{
    printf("Exporting formula table...")

    forms <- copy(formulaTable(formConsensus))
    if (retMin)
        forms[, ret := ret / 60]

    forms <- forms[group %in% names(fGroups)]

    write.csv(forms, file.path(path, "formulas.csv"))

    printf("Done!\n")
}

reportFormulaSpectra <- function(fGroups, path, formConsensus, MSPeakLists, EICRtWindow, EICMzWindow, retMin,
                                 EICTopMost, EICs)
{
    printf("Exporting formula MS/MS spectra...\n")

    gInfo <- groupInfo(fGroups)
    anaInfo <- analysisInfo(fGroups)
    # specs <- loadAllSpectra(anaInfo$analysis, anaInfo$path)
    forms <- formulaTable(formConsensus)
    pLists <- peakLists(MSPeakLists)

    forms <- forms[byMSMS == TRUE & group %in% rownames(gInfo)]

    if (nrow(forms) == 0)
        return()

    formGroups <- unique(forms$group)
    fcount <- length(formGroups)
    prog <- txtProgressBar(0, fcount, style = 3)

    for (grp in formGroups)
    {
        grpi <- match(grp, rownames(gInfo))

        out <- file.path(path, sprintf("%s-%s.pdf", class(fGroups), grp))
        pdf(out, paper = "a4", pointsize = 10, width = 8, height = 11)

        plotEIC(fGroups[, grp], EICRtWindow, EICMzWindow, retMin, EICTopMost, EICs)

        layout(matrix(1:16, 4, 4, byrow = TRUE), widths = c(2, 1, 2, 1))

        for (precursor in unique(forms[group == grp, formula]))
        {
            plotSpec(formConsensus, precursor, grp, MSPeakLists)

            oldp <- par(mar = c(0, 2, 0, 0))
            plot(1:10, 1:10, ann = FALSE, xaxt = "n", yaxt = "n", xlab = "", ylab = "", type = "n", adj = 1, bty = "n") # empty dummy plot
            text(1, 5, paste0(getFormInfoList(formConsensus, precursor, grp), collapse = "\n"), adj = 0, cex = 0.8)
            par(oldp)
        }
        dev.off()

        setTxtProgressBar(prog, match(grp, formGroups))
    }

    setTxtProgressBar(prog, fcount)
    close(prog)
}

reportCompoundTable <- function(fGroups, path, compounds, normalizeScores)
{
    printf("Exporting identification results tables...")

    gInfo <- groupInfo(fGroups)
    gNames <- names(fGroups)
    anaInfo <- analysisInfo(fGroups)
    compTable <- compoundTable(compounds)

    mcn <- mergedCompoundNames(compounds)

    if (normalizeScores)
        compTable <- sapply(compTable, normalizeCompScores, mCompNames = mcn, simplify = FALSE)

    for (grp in names(compTable))
    {
        if (grp %in% gNames && nrow(compTable[[grp]]) > 0)
        {
            out <- file.path(path, sprintf("%s-%s.csv", class(fGroups), grp))
            write.csv(compTable[[grp]][, -getAllCompCols("fragInfo", names(compTable[[grp]]), mergedCompoundNames(compounds)),
                                       with = FALSE], out)
        }
    }

    printf("Done!\n")
}

reportCompoundSpectra <- function(fGroups, path, MSPeakLists, compounds, formConsensus, EICRtWindow, EICMzWindow, retMin,
                                  EICTopMost, EICs, normalizeScores)
{
    printf("Exporting compound identification MS/MS spectra...\n")

    gNames <- names(fGroups)
    gInfo <- groupInfo(fGroups)
    anaInfo <- analysisInfo(fGroups)
    # specs <- loadAllSpectra(anaInfo$analysis, anaInfo$path)
    compTable <- compoundTable(compounds)
    pLists <- peakLists(MSPeakLists)

    if (!is.null(formConsensus))
        fTable <- formulaTable(formConsensus)

    idcount <- length(compTable)
    prog <- txtProgressBar(0, idcount, style = 3)

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
            pdf(out, paper = "a4r", pointsize = 10, width = 11, height = 8)

            plotEIC(fGroups[, fgrpi], EICRtWindow, EICMzWindow, retMin, EICTopMost, EICs)

            ncols <- 2
            for (idi in seq_len(nrow(compTable[[grp]])))
            {
                # NOTE: layout/mfrow/mfcol doesn't work because of the legend positioning (thinks 2 plots are made...)

                scr <- split.screen(c(2, 1))
                scr <- c(scr, split.screen(c(1, 2), screen = scr[2]))

                screen(scr[1])
                plotSpec(compounds, idi, grp, MSPeakLists = MSPeakLists)

                screen(scr[3])
                plotScores(compounds, idi, grp, normalizeScores)

                screen(scr[4])
                # draw text info
                oldp <- par(mar = c(0, 2, 0, 0))
                plot(1:10, 1:10, ann = FALSE, xaxt = "n", yaxt = "n", xlab = "", ylab = "", type = "n", adj = 1, bty = "n") # empty dummy plot
                text(1, 5, paste0(getCompInfoList(compTable[[grp]], idi, FALSE, mergedCompoundNames(compounds)),
                                  collapse = "\n"), adj = 0, cex = 0.8)
                par(oldp)

                close.screen(scr)
            }

            dev.off()
        }
        setTxtProgressBar(prog, gi)
    }

    setTxtProgressBar(prog, idcount)
    close(prog)
}

reportComponentTable <- function(components, path, retMin)
{
    printf("Exporting component table...")

    cTable <- rbindlist(componentTable(components), idcol = "component")

    if (retMin)
        cTable[, rt := rt / 60]
    write.csv(cTable, file.path(path, "components.csv"))

    printf("Done!\n")
}

reportComponentSpectra <- function(fGroups, path, components, EICRtWindow, EICMzWindow, retMin, EICs)
{
    printf("Exporting component spectra...\n")

    gInfo <- groupInfo(fGroups)
    anaInfo <- analysisInfo(fGroups)
    cInfo <- componentInfo(components)
    cTable <- componentTable(components)

    prog <- txtProgressBar(0, length(components), style = 3)

    pdf(file.path(path, "components.pdf"), paper = "a4", pointsize = 10, width = 8, height = 11)

    for (cmpi in seq_len(length(components)))
    {
        scr <- split.screen(c(2, 1))

        screen(scr[1])
        plotEIC(components, cmpi, fGroups, title = sprintf("Component %d", cmpi), rtWindow = EICRtWindow,
                mzWindow = EICMzWindow, retMin = retMin, EICs = EICs)

        screen(scr[2])
        plotSpec(components, cmpi, main = sprintf("ret: %.1f; m/z: %.4f - %.4f",
                                                  cInfo$ret[cmpi], min(cTable[[cmpi]]$mz), max(cTable[[cmpi]]$mz)))

        close.screen(scr)

        setTxtProgressBar(prog, cmpi)
    }

    dev.off()
    close(prog)
}

reportHClusterTable <- function(fGroups, path, cInfo, k)
{
    printf("Exporting hierarchical clustering table...")

    ct <- cutree(cInfo@clust, k)
    gInfo <- groupInfo(fGroups)

    gTable <- setnames(as.data.table(gInfo), c("rts", "mzs"), c("ret", "mz"))
    gTable[, group := rownames(gInfo)]
    gTable[, cluster := ct[group]]

    write.csv(gTable, file.path(path, "hcluster.csv"))
}

# UNDONE: eawag refs, call also from report?
reportHCluster <- function(fGroups, path, cInfo, k, silInfo = NULL)
{
    printf("Exporting hierarchical clustering plots...")

    out <- file.path(path, sprintf("%s-hcluster.pdf", class(fGroups)))
    pdf(out, paper = "a4", pointsize = 10, width = 8, height = 11)

    drawHeatMap(cInfo)

    if (!is.null(silInfo))
        plot(silInfo)

    plot(cInfo, k, cex = 0.3,
         labels = if (length(cInfo) <= clusterMaxLabels) NULL else FALSE)

    if (!is.null(k))
    {
        par(mfrow = c(2, 2))
        for (cl in seq_len(k))
        {
            plotInt(cInfo, k, cl, main = sprintf("Cluster %d - normalized", cl))
            plotInt(hClusterFilter(fGroups, cInfo, k, cl), TRUE, main = sprintf("Cluster %d - absolute", cl))
        }
    }

    dev.off()

    printf("Done!\n")
}

#' @details \code{reportCSV} generates tabular data (\emph{i.e.} \file{.csv}
#'   files) for given data to be reported. This may also be useful to allow
#'   import by other tools for post processing.
#'
#' @param reportFGroupsAsRows Report feature groups as rows (instead of columns)
#'   within the resulting \file{.csv} file.
#' @param reportFGroupsAnalysisInfo Include analyses information (reference and
#'   replicate groups) in the reported feature groups table \file{.csv} file.
#' @param reportFGroupsRetMz Include feature group information (retention time
#'   and \emph{m/z}) within the reported feature groups table \file{.csv} file.
#' @param reportFeatures If set to \code{TRUE} then for each analysis a
#'   \file{.csv} file will be generated with information about its detected
#'   features.
#'
#' @rdname reporting
#' @aliases reportCSV
#' @export
setMethod("reportCSV", "featureGroups", function(fGroups, path, reportFGroupsAsRows, reportFGroupsAnalysisInfo,
                                                 reportFGroupsRetMz, reportFeatures,
                                                 formConsensus, compounds, compoundNormalizeScores,
                                                 components, cInfo, clusterK, retMin, clearPath)
{
    prepareReportPath(path, clearPath)

    reportFGroupTable(fGroups, path, reportFGroupsAsRows, reportFGroupsAnalysisInfo, reportFGroupsRetMz, retMin)

    if (reportFeatures)
    {
        p <- file.path(path, "features")
        mkdirp(p)
        reportFeatureTable(fGroups, p, retMin)
    }

    if (!is.null(formConsensus) && length(formConsensus) > 0)
        reportFormulaTable(fGroups, formConsensus, path, retMin)

    if (!is.null(compounds))
    {
        p <- file.path(path, "compounds")
        mkdirp(p)
        reportCompoundTable(fGroups, p, compounds, compoundNormalizeScores)
    }

    if (!is.null(components))
        reportComponentTable(components, path, retMin)

    if (!is.null(cInfo))
        reportHClusterTable(fGroups, path, cInfo, clusterK)

})

#' @details \code{reportPDF} will report graphical data (\emph{e.g.}
#'   chromatograms and mass spectra) within PDF files. Compared
#'   to \code{reportMD} this function may be faster and yield smaller report
#'   files, however, its functionality is a bit more basic and generated data is
#'   more 'scattered' around.
#'
#' @param EICGrid An integer vector in the form \code{c(columns, rows)} that is
#'   used to determine the plotting grid when reporting EICs in PDF files.
#'
#' @rdname reporting
#' @aliases reportPDF
#' @export
setMethod("reportPDF", "featureGroups", function(fGroups, path, reportFGroups,
                                                 formConsensus, reportFormulaSpectra,
                                                 compounds, compoundNormalizeScores,
                                                 components, cInfo, clusterK, silInfo, clusterMaxLabels,
                                                 MSPeakLists, retMin, EICGrid, EICRtWindow, EICMzWindow,
                                                 EICTopMost, EICOnlyPresent, clearPath)
{
    if (!reportFGroups && is.null(formConsensus) && is.null(compounds) && is.null(components) &&
        is.null(cInfo))
    {
        cat("Nothing to report...\n")
        return(NULL)
    }
    
    if (is.null(MSPeakLists) &&
        ((!is.null(formConsensus) && reportFormulaSpectra) || !is.null(compounds)))
        stop("MSPeakLists is NULL, please specify when reporting formula and/or compounds")

    prepareReportPath(path, clearPath)

    if (reportFGroups || !is.null(formConsensus) || !is.null(compounds) || !is.null(components))
    {
        cat("Loading all EICs... ")
        EICs <- loadXCMSEICForFGroups(fGroups, EICRtWindow, EICMzWindow, EICTopMost, EICOnlyPresent)
        cat("Done!\n")
    }

    if (reportFGroups)
        reportFGroupPlots(fGroups, path, EICGrid, EICRtWindow, EICMzWindow, retMin, EICTopMost, EICs)

    if (reportFormulaSpectra && !is.null(formConsensus) && length(formConsensus) > 0)
    {
        p <- file.path(path, "formulas")
        mkdirp(p)
        reportFormulaSpectra(fGroups, p, formConsensus, MSPeakLists, EICRtWindow, EICMzWindow, retMin, EICTopMost, EICs)
    }

    if (!is.null(compounds))
    {
        p <- file.path(path, "compounds")
        mkdirp(p)
        reportCompoundSpectra(fGroups, p, MSPeakLists, compounds, formConsensus, EICRtWindow, EICMzWindow, retMin,
                              EICTopMost, EICs, compoundNormalizeScores)
    }

    if (!is.null(components))
        reportComponentSpectra(fGroups, path, components, EICRtWindow, EICMzWindow, retMin, EICs)

    if (!is.null(cInfo))
    {
        if (is.null(clusterK))
            stop("Please specify clusterK argument!")
        reportHCluster(fGroups, path, cInfo, clusterK, silInfo, clusterMaxLabels)
    }
})

#' @details \code{reportMD} will report graphical data (\emph{e.g.}
#'   chromatograms and mass spectra) and summary information in an easy
#'   browsable \code{HTML} file using \link{rmarkdown}, \link{flexdashboard} and
#'   \link{knitr}.
#'
#' @param reportChord If \code{TRUE} then a chord diagram for all feature groups
#'   is plotted (data will be averaged among replicates).
#' @param interactiveHeat If \code{TRUE} an interactive heatmap HTML widget will
#'   be generated to display hierarchical clustering results. Set to
#'   \code{FALSE} for a 'regular' static plot.
#' @param selfContained If \code{TRUE} the output will be a standalone HTML file
#'   which contains all graphics and script dependencies. When \code{FALSE}, the
#'   latter will be placed in an additional directory (\file{report_files})
#'   which should remain present when viewing the output file. Especially on
#'   Windows, a non-self contained output might be desirable when reporting
#'   large amounts of data to prevent \command{pandoc} from running out of
#'   memory.
#' @param optimizePng If \code{TRUE} then \command{pngquant} is used to reduce
#'   the size of generated graphics. A signficant reduction in disk space usage
#'   may be seen, however, at the cost additional processing time.
#' @param maxProcAmount Maximum amount of \command{pngquant} commands to run in
#'   parallel. Higher numbers will decrease processing time, with an optimum
#'   usually close to the amount of CPU cores.
#'
#' @references \addCitations{knitr}{2} \cr\cr
#'   \addCitations{knitr}{3}
#'
#' @rdname reporting
#' @aliases reportMD
#' @export
setMethod("reportMD", "featureGroups", function(fGroups, path, reportChord, reportFGroups,
                                                formConsensus, reportFormulaSpectra,
                                                compounds, compoundNormalizeScores,
                                                components, cInfo, clusterK, silInfo,
                                                interactiveHeat, clusterMaxLabels, MSPeakLists, retMin,
                                                EICRtWindow, EICMzWindow, EICTopMost, EICOnlyPresent,
                                                selfContained, optimizePng, maxProcAmount, clearPath)
{
    # UNDONE: mention pandoc win limits
    
    if (is.null(MSPeakLists) &&
        ((!is.null(formConsensus) && reportFormulaSpectra) || !is.null(compounds)))
        stop("MSPeakLists is NULL, please specify when reporting formula and/or compounds")

    prepareReportPath(path, clearPath)

    workPath <- file.path(tempdir(), "report")
    unlink(workPath, TRUE)
    mkdirp(workPath)

    file.copy(system.file("report", "main.Rmd", package = "patRoon"), workPath)
    file.copy(system.file("report", "spectra.Rmd", package = "patRoon"), workPath)
    file.copy(system.file("report", "components.Rmd", package = "patRoon"), workPath)
    file.copy(system.file("report", "cluster.Rmd", package = "patRoon"), workPath)

    # rmarkdown needs absolute path as relative paths will be from the path of the Rmd
    if (!R.utils::isAbsolutePath(path))
        path <- R.utils::getAbsolutePath(path)

    EICs <- NULL

    # NOTE: loading ECIs is done outside knitr environment: the current working
    # directory is changed so we loose access to the cache.

    # UNDONE: always load EICs? Need them for summary EIC plot
    # if (reportFGroups || (!is.null(formConsensus) && reportFormulaSpectra) ||
    #     !is.null(compounds) || !is.null(components))
    {
        cat("Loading all EICs... ")
        EICs <- loadXCMSEICForFGroups(fGroups, EICRtWindow, EICMzWindow, EICTopMost, EICOnlyPresent)
        cat("Done!\n")
    }

    rmdVars <- list(outPath = path, fGroups = fGroups, groupNames = names(fGroups), gInfo = groupInfo(fGroups), reportChord = reportChord,
                    reportFGroups = reportFGroups, EICRtWindow = EICRtWindow, EICMzWindow = EICMzWindow,
                    retMin = retMin, EICTopMost = EICTopMost, EICOnlyPresent = EICOnlyPresent, EICs = EICs,
                    compounds = compounds, MSPeakLists = MSPeakLists, formConsensus = formConsensus,
                    compoundNormalizeScores = compoundNormalizeScores, components = components,
                    cInfo = cInfo, clusterK = clusterK, silInfo = silInfo, interactiveHeat = interactiveHeat,
                    clusterMaxLabels = clusterMaxLabels, selfContained = selfContained,
                    optimizePng = optimizePng, maxProcAmount = maxProcAmount)

    rmarkdown::render(file.path(workPath, "main.Rmd"), output_file = file.path(path, "report.html"),
                      output_options = list(self_contained = selfContained),
                      quiet = TRUE)
})
