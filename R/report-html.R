#' @include main.R
NULL

reportHTMLUtils <- setRefClass("reportHTMLUtils",
                               fields = list(objects = "list", EICs = "list", plots = "list", properties = "list"))

reportHTMLUtils$methods(
    hasSuspects = function() isScreening(objects$fGroups),
    hasSets = function() isFGSet(objects$fGroups),
    hasFQualities = function() hasFGroupScores(objects$fGroups),
    hasComponents = function() !is.null(objects[["components"]]) && !inherits(objects$components, "componentsTPs"),
    hasComponentsIntClust = function() hasComponents() && inherits(objects$components, "componentsIntClust"),
    hasComponentsSpecClust = function() hasComponents() && inherits(objects$components, "componentsSpecClust"),
    hasComponentsNT = function() hasComponents() && inherits(objects$components, "componentsNT"),
    hasTPs = function() !is.null(objects[["components"]]) && inherits(objects$components, "componentsTPs"),
    hasFormulas = function() !is.null(objects[["formulas"]]),
    hasCompounds = function() !is.null(objects[["compounds"]]),
    hasCompsCluster = function() !is.null(objects[["compsCluster"]])
)

# UNDONE: method
#' @export
reportHTMLNew <- function(fGroups, path = "report", MSPeakLists = NULL, formulas = NULL, compounds = NULL,
                          compsCluster = NULL, components = NULL, TPs = NULL, selfContained = FALSE, clearPath = FALSE,
                          openReport = TRUE, noDate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertPathForOutput(path, overwrite = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ selfContained + clearPath + openReport + noDate,
           fixed = list(add = ac))
    aapply(checkmate::assertClass, . ~ MSPeakLists + formulas + compounds + compsCluster + components + TPs,
           c("MSPeakLists", "formulas", "compounds", "compoundsCluster", "components", "transformationProducts"),
           null.ok = TRUE, fixed = list(add = ac))
    # assertSpecSimParams(specSimParams, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(fGroups) == 0)
    {
        cat("No feature groups, nothing to report...\n")
        return(invisible(NULL))
    }
    
    if (is.null(MSPeakLists) && (!is.null(formulas) || !is.null(compounds)))
        stop("MSPeakLists is NULL, please specify when reporting formula and/or compounds")
    
    prepareReportPath(path, clearPath)
    
    workPath <- file.path(tempdir(TRUE), "report")
    unlink(workPath, TRUE)
    mkdirp(workPath)
    
    file.copy(system.file("report", "report.Rmd", package = "patRoon"), workPath)
    file.copy(system.file("report", "details.Rmd", package = "patRoon"), workPath)
    file.copy(system.file("report", "components_int.Rmd", package = "patRoon"), workPath)
    file.copy(system.file("report", "components_spec.Rmd", package = "patRoon"), workPath)
    file.copy(system.file("report", "components_nt.Rmd", package = "patRoon"), workPath)
    file.copy(system.file("report", "modal.html", package = "patRoon"), workPath)
    
    # rmarkdown needs absolute path as relative paths will be from the path of the Rmd
    path <- normalizePath(path)

    cat("Loading all EICs... ")
    # UNDONE: params
    # EICs <- getEICsForFGroups(fGroups, 30, 0.002, 1, TRUE, TRUE)
    EICs <- getEICsForFGroups(fGroups, 30, 0.002, NULL, FALSE, FALSE)
    cat("Done!\n")
    
    reportEnv <- new.env()
    
    reportEnv$properties <- list(noDate = noDate)
    reportEnv$plots <- generateHTMLReportPlots(fGroups, MSPeakLists, formulas, compounds, compsCluster, components,
                                               TPs, path, EICs, selfContained)
    reportEnv$utils <- reportHTMLUtils$new(objects = list(fGroups = fGroups, MSPeakLists = MSPeakLists,
                                                          formulas = formulas, compounds = compounds,
                                                          compsCluster = compsCluster, components = components,
                                                          TPs = TPs),
                                           EICs = EICs, plots = reportEnv$plots,
                                           properties = list(selfContained = selfContained))
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
                                          output_options = list(self_contained = selfContained),
                                          quiet = TRUE, envir = reportEnv))
    
    if (openReport)
        utils::browseURL(paste0("file://", normalizePath(outputFile)))
    
    invisible(NULL)
}
