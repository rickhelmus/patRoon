#' @include main.R
NULL

reportHTMLUtils <- setRefClass("reportHTMLUtils",
                               fields = list(objects = "list", EICs = "list", plots = "list", settings = "list"))

reportHTMLUtils$methods(
    hasSuspects = function() isScreening(objects$fGroups) && nrow(screenInfo(objects$fGroups)) > 0,
    hasSuspAnn = function() isSuspAnnotated(objects$fGroups) && nrow(screenInfo(objects$fGroups)) > 0,
    hasInternalStandards = function() nrow(internalStandards(objects$fGroups)) > 0,
    hasSets = function() isFGSet(objects$fGroups),
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
    hasTPGraphs = function() hasTPs() && inherits(objects$TPs, c("transformationProductsStructure", "transformationProductsFormula")),
    hasMSPL = function() hasObj("MSPeakLists"),
    hasFormulas = function() hasObj("formulas"),
    hasCompounds = function() hasObj("compounds"),
    hasCompsCluster = function() hasObj("compsCluster"),
    
    getFGSets = function() sets(objects$fGroups)
)

# UNDONE: method
#' @export
reportHTMLNew <- function(fGroups, MSPeakLists = NULL, formulas = NULL, compounds = NULL,
                          compsCluster = NULL, components = NULL, TPs = NULL,
                          settingsFile = system.file("report", "settings.yml", package = "patRoon"),
                          path = NULL, EICParams = getDefEICParams(topMost = 1, topMostByRGroup = TRUE),
                          clearPath = FALSE, openReport = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    if (!is.null(path))
        checkmate::assertPathForOutput(path, overwrite = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ clearPath + openReport, fixed = list(add = ac))
    aapply(checkmate::assertClass, . ~ MSPeakLists + formulas + compounds + compsCluster + components + TPs,
           c("MSPeakLists", "formulas", "compounds", "compoundsCluster", "components", "transformationProducts"),
           null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertFileExists(settingsFile, "r", add = ac)
    # assertSpecSimParams(specSimParams, add = ac)
    checkmate::reportAssertions(ac)
    
    if (is.null(MSPeakLists) && (!is.null(formulas) || !is.null(compounds)))
        stop("MSPeakLists is NULL, please specify when reporting formula and/or compounds")

    settings <- readYAML(settingsFile)
    assertReportSettings(settings)
    
    if (is.null(path))
        path <- settings$general$path
    
    if (length(fGroups) == 0)
    {
        cat("No feature groups, nothing to report...\n")
        return(invisible(NULL))
    }
    
    prepareReportPath(path, clearPath)
    
    workPath <- file.path(tempdir(TRUE), "report")
    unlink(workPath, TRUE)
    mkdirp(workPath)
    
    file.copy(system.file("report", c("report.Rmd", "details.Rmd", "istds.Rmd", "components_int.Rmd",
                                      "components_spec.Rmd", "components_nt.Rmd", "modal.html"), package = "patRoon"),
              workPath)
    
    # rmarkdown needs absolute path as relative paths will be from the path of the Rmd
    path <- normalizePath(path)

    cat("Loading all EICs... ")
    # UNDONE: make onlyPresent configurable, check if feature EICs are needed
    EICs <- getEICsForFGroups(fGroups, EICParams = modifyList(EICParams, list(topMost = NULL, onlyPresent = FALSE),
                                                              keep.null = TRUE))
    cat("Done!\n")
    
    reportEnv <- new.env()
    
    reportEnv$settings <- settings
    reportEnv$plots <- generateHTMLReportPlots(fGroups, MSPeakLists, formulas, compounds, compsCluster, components,
                                               TPs, settings, path, EICs, EICParams)
    reportEnv$utils <- reportHTMLUtils$new(objects = list(fGroups = fGroups, MSPeakLists = MSPeakLists,
                                                          formulas = formulas, compounds = compounds,
                                                          compsCluster = compsCluster, components = components,
                                                          TPs = TPs),
                                           EICs = EICs, plots = reportEnv$plots,
                                           settings = settings)
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
    
    myPlotFiles <- normalizePath(file.path(path, unlist(reportEnv$plots)))
    allPlotFiles <- normalizePath(list.files(file.path(path, "report_files", "plots"), pattern = "\\.svg$",
                                             full.names = TRUE))
    oldPlotFiles <- setdiff(allPlotFiles, myPlotFiles)
    opfDates <- lapply(file.info(oldPlotFiles)$mtime, as.Date)
    opfAge <- lapply(opfDates, difftime, time1 = Sys.Date(), units = "days")
    opfAge <- sapply(opfAge, as.numeric)
    oldPlotFiles <- oldPlotFiles[opfAge >= settings$general$plotFileRetention]
    if (length(oldPlotFiles) > 0)
    {
        file.remove(oldPlotFiles)
        printf("Removed %d old plot files.\n", length(oldPlotFiles))
    }
    
    invisible(NULL)
}
