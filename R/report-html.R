#' @include main.R
NULL

makeHTMLReportPlot <- function(out, outPath, selfContained, ...)
{
    destPath <- if (selfContained) "." else file.path(outPath, "report_files", "plots")
    mkdirp(destPath)
    out <- file.path(destPath, out)
    withr::with_svg(out, ...)
    return(out)
}

generateReportPlots <- function(fGroups, MSPeakLists, formulas, compounds, components, TPs, outPath, EICs,
                                selfContained)
{
    ret <- list()
    
    ret$overview$chroms <- makeHTMLReportPlot("chroms.svg", outPath, selfContained, {
        par(mai = c(0.9, 0.8, 0.6, 0.1))
        # UNDONE: params
        plotChroms(fGroups, 30, 0.005, TRUE, 1, FALSE, EICs, TRUE, FALSE, colourBy = "fGroups", showLegend = FALSE,
                   onlyPresent = TRUE)
    })
    
    return(ret)
}

# UNDONE: method
#' @export
reportHTMLNew <- function(fGroups, path = "report", MSPeakLists = NULL, formulas = NULL, compounds = NULL,
                          components = NULL, TPs = NULL, selfContained = FALSE, clearPath = FALSE, openReport = TRUE,
                          noDate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertPathForOutput(path, overwrite = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ selfContained + clearPath + openReport + noDate,
           fixed = list(add = ac))
    aapply(checkmate::assertClass, . ~ MSPeakLists + formulas + compounds + components + TPs,
           c("MSPeakLists", "formulas", "compounds", "components", "transformationProducts"),
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
    
    # rmarkdown needs absolute path as relative paths will be from the path of the Rmd
    path <- normalizePath(path)

    cat("Loading all EICs... ")
    # UNDONE: params
    EICs <- getEICsForFGroups(fGroups, 30, 0.002, 1, TRUE, TRUE)
    cat("Done!\n")
    
    reportEnv <- new.env()
    
    reportEnv$properties <- list(noDate = noDate)
    reportEnv$plots <- generateReportPlots(fGroups, MSPeakLists, formulas, compounds, components, TPs, path, EICs,
                                           selfContained)
    
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
