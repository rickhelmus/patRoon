#' @include main.R
NULL

makeHTMLReportPlot <- function(out, outPath, selfContained, code, ...)
{
    if (FALSE)
    {
        if (selfContained)
        {
            svgstr <- svglite::svgstring(standalone = FALSE, fix_text_size = FALSE, ...)
            on.exit(dev.off(), add = TRUE)
            force(code)
            ret <- as.character(svgstr())
            
            # replace fixed width/height properties to allow proper scaling (see https://stackoverflow.com/a/45144857).
            # NOTE: use sub so that only header (first occurrence) is modified. Furthermore, note that the svglite css class
            # is changed in report.Rmd.
            # ret <- sub("width=\\'[[:graph:]]+\\'", "width='100%'", ret)
            # ret <- sub("height=\\'[[:graph:]]+\\'", "height='auto'", ret)
            return(ret)
        }
    }
    else if (selfContained)
    {
        # UNDONE: while embedding the SVG directly would be nice, this seems to give major headaches with scaling,
        # especially with Firefox... For now just base64 it :(
        withSVGLite(out, standalone = TRUE, code = code, ...)
        return(paste0("<img src=", knitr::image_uri(out), "></img>"))
    }
    
    destPath <- file.path(outPath, "report_files", "plots")
    mkdirp(destPath)
    out <- file.path(destPath, out)
    withSVGLite(out, standalone = TRUE, code = code, ...)
    return(paste0("<img src='", out, "'></img>"))
    
    # UNDONE: object tag makes text selectable but messes up layout...
    # return(paste0("<object data='", out, "' type='image/svg+xml' width=500 height=300></object>"))
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
    }, width = 10, height = 4)
    
    ret$overview$retMZ <- makeHTMLReportPlot("retmz.svg", outPath, selfContained, {
        par(mai = c(0.9, 0.8, 0.1, 0.1))
        # UNDONE: params
        plot(fGroups, colourBy = "fGroups", showLegend = FALSE, retMin = TRUE)
    }, width = 10, height = 4)
    
    rGroupLenNonEmpty <- length(replicateGroups(removeEmptyAnalyses(fGroups)))
    rGroupLen <- length(replicateGroups(fGroups))
    anyOverlap <- rGroupLen > 1 &&
        length(unique(fGroups, which = replicateGroups(fGroups), outer = TRUE)) < length(fGroups)
    
    ret$overview$chord <- ret$overview$venn  <- ret$overview$UpSet <- NULL
    if (anyOverlap && rGroupLenNonEmpty > 1)
    {
        if (rGroupLenNonEmpty > 2)
        {
            ret$overview$chord <- makeHTMLReportPlot("chord.svg", outPath, selfContained, {
                # UNDONE: params(?)
                plotChord(fGroups, average = TRUE)
            }, width = 7, height = 7)
        }
        if (rGroupLen < 6)
        {
            ret$overview$venn <- makeHTMLReportPlot("venn.svg", outPath, selfContained, {
                # UNDONE: params(?)
                plotVenn(fGroups)
            }, width = 7, height = 7)
        }
        
        # UpSet
        ret$overview$UpSet <- makeHTMLReportPlot("upset.svg", outPath, selfContained, {
            # UNDONE: params(?)
            print(plotUpSet(fGroups))
        }, width = 7, height = 7)
    }
    
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
