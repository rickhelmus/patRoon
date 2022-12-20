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
