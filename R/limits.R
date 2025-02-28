getLimitsFile <- function()
{
    for (p in c(getOption("patRoon.path.limits", ""), getwd(), system.file("misc", package = "patRoon")))
    {
        if (!nzchar(p))
            next
        fp <- file.path(p, "limits.yml")
        if (file.exists(fp))
            return(fp)
    }
}

#' @export
defaultLim <- function(category, level)
{
    limits <- readYAML(getLimitsFile())
    limits <- assertAndPrepareLimits(limits)
    
    checkmate::assertChoice(category, setdiff(names(limits), "general"))
    checkmate::assertChoice(level, names(limits[[category]]))
    
    return(limits[[category]][[level]])
}

#' @export
genLimitsFile <- function(outPath = normalizePath("."), IMS = "bruker")
{
    checkmate::assertString(outPath, min.chars = 1)
    checkmate::assertChoice(IMS, c("bruker", "agilent"))
    
    fp <- file.path(outPath, "limits.yml")
    checkmate::assertPathForOutput(fp, overwrite = TRUE)
    
    limits <- readYAML(system.file("misc", "limits.yml", package = "patRoon"))
    limits$general$IMS <- IMS
    writeYAML(limits, fp)

    return(invisible(NULL))    
}
