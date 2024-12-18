# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
NULL


generateHTMLReportPlots <- function(fGroups, MSPeakLists, formulas, compounds, compsCluster, components, TPs, settings,
                                    outPath, EICs, EICParams, specSimParams, parallel)
{
    ret <- list()
    
    cat("Genarate summary plots...")
    
    ret$overview$chroms <- makeHTMLReportPlot("chroms", outPath, "plotChroms",
                                              list(fGroups, retMin = settings$features$retMin, EICs = EICs,
                                                   EICParams = modifyList(EICParams, list(topMost = 1,
                                                                                          topMostByRGroup = FALSE,
                                                                                          onlyPresent = TRUE)),
                                                   showPeakArea = TRUE, showFGroupRect = FALSE, colourBy = "fGroups",
                                                   showLegend = FALSE, intMax = settings$features$chromatograms$intMax),
                                              parParams = list(mai = c(0.9, 0.8, 0.6, 0.1)), width = 10, height = 4)
    
    ret$overview$retMZ <- makeHTMLReportPlot("retmz", outPath, "plot",
                                             list(fGroups, colourBy = "fGroups", showLegend = FALSE,
                                                  retMin = settings$features$retMin),
                                             parParams = list(mai = c(0.9, 0.8, 0.1, 0.1)), width = 10, height = 4)
    
    rGroupLenNonEmpty <- length(replicateGroups(removeEmptyAnalyses(fGroups)))
    rGroupLen <- length(replicateGroups(fGroups))
    anyOverlap <- rGroupLen > 1 &&
        length(unique(fGroups, which = replicateGroups(fGroups), outer = TRUE)) < length(fGroups)
    
    ret$overview$chord <- ret$overview$venn  <- ret$overview$UpSet <- NULL
    if (anyOverlap && rGroupLenNonEmpty > 1)
    {
        if ("chord" %in% settings$summary && rGroupLenNonEmpty > 2)
        {
            ret$overview$chord <- makeHTMLReportPlot("chord", outPath, "plotChord", list(fGroups, average = TRUE),
                                                     width = 7, height = 7)
        }
        if ("venn" %in% settings$summary && rGroupLen < 6)
        {
            ret$overview$venn <- makeHTMLReportPlot("venn", outPath, "plotVenn", list(fGroups), width = 7, height = 7)
        }
        
        if ("upset" %in% settings$summary)
            ret$overview$UpSet <- makeHTMLReportPlot("upset", outPath, "plotUpSet", list(fGroups), doPrint = TRUE,
                                                     width = 7, height = 7)
    }
    cat(" Done!\n")
    
    ret$chromsLarge <- genHTMLReportPlotsChromsLarge(fGroups, settings, outPath, EICs, EICParams, parallel)
    ret$chromsSmall <- genHTMLReportPlotsChromsSmall(fGroups, settings, outPath, EICs, EICParams, parallel)
    ret$chromsFeatures <- genHTMLReportPlotsChromsFeatures(fGroups, settings, outPath, EICs, EICParams, parallel)
    
    ret$intPlots <- genHTMLReportPlotsIntPlots(fGroups, settings, outPath, parallel)
    
    gNames <- names(fGroups)
    if (!is.null(compounds))
        compounds <- compounds[gNames]
    
    ret$structs <- genHTMLReportPlotsStructs(fGroups, compounds, settings, outPath, parallel)
    
    if (!is.null(MSPeakLists))
        ret$MSPeakLists <- genHTMLReportPlotsMSPeakLists(MSPeakLists[, gNames], settings, outPath, parallel)    
    if (!is.null(formulas))
        ret$formulas <- genHTMLReportPlotsFormulas(formulas[gNames], MSPeakLists, settings, outPath, parallel)
    if (!is.null(compounds))
        ret$compounds <- genHTMLReportPlotsCompounds(compounds, MSPeakLists, formulas, settings, outPath, parallel)
    if (!is.null(compsCluster))
        ret$compsCluster <- genHTMLReportPlotsCompsCluster(compsCluster[gNames], settings, outPath, parallel)
    if (!is.null(components))
    {
        if (!inherits(components, "componentsTPs"))
            ret$components <- genHTMLReportPlotsComponents(fGroups, components, settings, outPath, EICs, EICParams,
                                                           parallel)
        else
            ret$TPs <- genHTMLReportPlotsTPs(fGroups, components, MSPeakLists, formulas, compounds, settings,
                                             specSimParams, outPath, EICs, parallel)
    }
    return(ret)
}


reportHTMLUtils <- setRefClass("reportHTMLUtils",
                               fields = list(objects = "list", EICs = "list", plots = "list", settings = "list"))

reportHTMLUtils$methods(
    hasSuspects = function() isScreening(objects$fGroups) && nrow(screenInfo(objects$fGroups)) > 0,
    hasSuspAnn = function() isSuspAnnotated(objects$fGroups) && nrow(screenInfo(objects$fGroups)) > 0,
    hasInternalStandards = function() nrow(internalStandards(objects$fGroups)) > 0,
    hasSets = function() isFGSet(objects$fGroups),
    hasConcs = function() nrow(concentrations(objects$fGroups)) > 0,
    hasTox = function() nrow(toxicities(objects$fGroups)) > 0,
    hasFQualities = function() hasFGroupScores(objects$fGroups),
    hasObj = function(name) !is.null(objects[[name]]) && length(intersect(names(objects$fGroups), groupNames(objects[[name]]))) > 0,
    hasComponents = function() hasObj("components") && !inherits(objects$components, "componentsTPs"),
    hasComponentsIntClust = function() hasComponents() && inherits(objects$components, "componentsIntClust"),
    hasComponentsSpecClust = function() hasComponents() && inherits(objects$components, "componentsSpecClust"),
    hasComponentsNT = function() hasComponents() && inherits(objects$components, c("componentsNT", "componentsNTSet")),
    hasComponentInfo = function() hasComponents() && !hasComponentsIntClust() && !hasComponentsSpecClust(),
    hasComponentsTPs = function() hasObj("components") && inherits(objects$components, "componentsTPs"),
    hasComponentsFromTPs = function() hasComponentsTPs() && objects$components@fromTPs,
    hasTPSims = function() hasComponentsTPs() && any(c("specSimilarity", "fragmentMatches") %in% names(as.data.table(objects$components))),
    hasTPs = function() !is.null(objects[["TPs"]]) && hasComponentsTPs(),
    hasTPGraphs = function() hasTPs() && inherits(objects$TPs, c("transformationProductsStructure", "transformationProductsFormula")) && settings$TPs$graphs,
    hasMSPL = function() hasObj("MSPeakLists"),
    hasFormulas = function() hasObj("formulas") && settings$formulas$include,
    hasCompounds = function() hasObj("compounds"),
    hasCompsCluster = function() hasObj("compsCluster"),
    
    getFGSets = function() sets(objects$fGroups),
    
    plotImg = function(p) paste0("<img src='", p, "'></img>")
)


doReportHTML <- function(fGroups, MSPeakLists, formulas, compounds, compsCluster, components, TPs, settings, path,
                         EICParams, specSimParams, openReport, parallel)
{
    workPath <- file.path(tempdir(TRUE), "report")
    unlink(workPath, TRUE)
    mkdirp(workPath)
    
    file.copy(system.file("report", c("report.Rmd", "details.Rmd", "istds.Rmd", "components_int.Rmd",
                                      "components_spec.Rmd", "components_nt.Rmd", "modal.html"), package = "patRoon"),
              workPath)
    
    # rmarkdown needs absolute path as relative paths will be from the path of the Rmd
    path <- normalizePath(path)

    cat("Loading all EICs... ")
    EICParams$onlyPresent <- settings$features$chromatograms$features != "all"
    if (!settings$features$chromatograms$large && isFALSE(settings$features$chromatograms$features) &&
        settings$features$chromatograms$small)
    {
        # plot only small chromatograms (and summary overview), get minimum set of EICs
        EICParams$topMost <- 1
        EICParams$topMostByRGroup <- FALSE
    }
    else
        EICParams["topMost"] <- list(NULL)
    
    EICs <- getEICsForFGroups(fGroups, EICParams = EICParams)
    cat("Done!\n")
    
    allPlots <- generateHTMLReportPlots(fGroups, MSPeakLists, formulas, compounds, compsCluster, components, TPs,
                                        settings, path, EICs, EICParams, specSimParams, parallel)
    
    reportEnv <- new.env()
    
    reportEnv$settings <- settings
    
    reportEnv$plots <- rapply(allPlots, function(p)
    {
        wh <- which(nzchar(p))
        if (length(wh) == 0)
        {} # nothing
        else if (settings$general$selfContained)
            p[wh] <- sapply(p[wh], knitr::image_uri)
        else
            p[wh] <- file.path("report_files", "plots", basename(p[wh])) # make paths relative for correct HTML links
        return(p)
    }, how = "replace")
    
    reportEnv$utils <- reportHTMLUtils$new(objects = list(fGroups = fGroups, MSPeakLists = MSPeakLists,
                                                          formulas = formulas, compounds = compounds,
                                                          compsCluster = compsCluster, components = components,
                                                          TPs = TPs),
                                           EICs = EICs, plots = reportEnv$plots, settings = settings)
    reportEnv$EICs <- EICs
    
    reportEnv$objectsShow <- paste0(utils::capture.output({
        for (o in pruneList(list(fGroups, MSPeakLists, formulas, compounds, components, TPs)))
        {
            show(o)
            cat("\n")
        }
    }), collapse = "\n")
    
    # HACK: not sure what exactly happens here, but... kableExtra adds latex
    # dependencies by default, which then may cause serious memory leakage when
    # rmarkdown::render() is called repeatedly. For now just remove them temporarily.
    knitMeta <- knitr::knit_meta("latex_dependency", clean = TRUE)
    on.exit(knitr::knit_meta_add(knitMeta), add = TRUE)
    
    outputFile <- file.path(path, "report.html")
    
    # normalize cache path so it can be used in report working directory
    withr::with_options(list(patRoon.cache.fileName = normalizePath(getOption("patRoon.cache.fileName")),
                             patRoon.progress.opts = list(file = stderr())),
                        rmarkdown::render(file.path(workPath, "report.Rmd"), output_file = outputFile,
                                          output_options = list(self_contained = settings$general$selfContained),
                                          quiet = TRUE, envir = reportEnv))
    
    if (openReport)
        utils::browseURL(paste0("file://", normalizePath(outputFile)))
    
    usedPlotFiles <- unlist(allPlots)
    usedPlotFiles <- usedPlotFiles[nzchar(usedPlotFiles)]
    Sys.setFileTime(usedPlotFiles, Sys.time()) # update file times in case plot files already exist and were re-used
    
    allPlotFiles <- normalizePath(list.files(file.path(path, "report_files", "plots"), pattern = "\\.svg$",
                                             full.names = TRUE))
    oldPlotFiles <- setdiff(allPlotFiles, usedPlotFiles)
    opfDates <- lapply(file.info(oldPlotFiles)$mtime, as.Date)
    opfAge <- lapply(opfDates, difftime, time1 = Sys.Date(), units = "days")
    opfAge <- sapply(opfAge, as.numeric)
    oldPlotFiles <- oldPlotFiles[opfAge >= settings$general$keepUnusedPlots]
    if (length(oldPlotFiles) > 0)
    {
        file.remove(oldPlotFiles)
        printf("Removed %d old plot files.\n", length(oldPlotFiles))
    }
    
    invisible(NULL)
}
