#' @include main.R
#' @include components-features.R
NULL

#' @export
componentsOpenMS <- setClass("componentsOpenMS", contains = "componentsFeatures")

setMethod("initialize", "componentsOpenMS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "openms", ...))

#' @rdname component-generation
#' @export
setMethod("generateComponentsOpenMS", "featureGroups", function(fGroups, ionization, chargeMax = 1,
                                                                mzWindow = 0.005, extraOpts = NULL)
{
    # UNDONE: all features are currently annotated (ie including not in a group), should be fine once featng is merged
    # UNDONE: convert adduct format
    # UNDONE: get neutral masses
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertChoice(ionization, c("positive", "negative"), add = ac)
    checkmate::reportAssertions(ac)
    
    anaInfo <- analysisInfo(fGroups)
    featTable <- featureTable(fGroups)
    
    printf("Annotating all features with OpenMS for %d analyses ...\n", nrow(anaInfo))
    
    params <- list(ionization = ionization, chargeMax = chargeMax, mzWindow = mzWindow, extraOpts = extraOpts)
    paramsHash <- makeHash(params)
    
    cmdQueue <- lapply(seq_len(nrow(anaInfo)), function(anai)
    {
        hash <- makeHash(featTable[[anai]], paramsHash)
        logf <- paste0("mad-", anaInfo$analysis[anai], ".txt")
        return(list(hash = hash, fTable = featTable[[anai]], logFile = logf))
    })
    names(cmdQueue) <- anaInfo$analysis
    
    if (length(cmdQueue) == 0)
        return(componentsOpenMS())
    
    featComponents <- executeMultiProcess(cmdQueue, function(cmd, ...)
    {
        fcmp <- patRoon:::parseAdductConsXMLFile(cmd$outFile)
        unlink(cmd[c("inFile", "outFile")]) # remove temporary files, as their size may be considerable
        return(fcmp)
    }, prepareHandler = function(cmd)
    {
        inFile <- tempfile(fileext = ".featureXML")
        writeFeatureXML(cmd$fTable, inFile)
        outFile <- tempfile(fileext = ".consensusXML")
        cmdMAD <- do.call(patRoon:::getOpenMSMADCommand, c(list(inFile = inFile, outFile = outFile), params))
        return(c(cmd, list(outFile = outFile), cmdMAD))
    }, logSubDir = "openms", cacheName = "componentsOpenMS")
        
    return(componentsOpenMS(fGroups = fGroups, featureComponents = featComponents))
})


getOpenMSMADCommand <- function(inFile, outFile, ionization, chargeMax, mzWindow, extraOpts)
{
    boolToChr <- function(b) if (b) "true" else "false"
    
    settings <- list("-algorithm:MetaboliteFeatureDeconvolution:negative_mode" = boolToChr(ionization == "negative"),
                     "-algorithm:MetaboliteFeatureDeconvolution:charge_max" = chargeMax,
                     "-algorithm:MetaboliteFeatureDeconvolution:mass_max_diff" = mzWindow)
    
    if (!is.null(extraOpts))
        settings <- modifyList(settings, extraOpts)
    
    return(list(command = getCommandWithOptPath("MetaboliteAdductDecharger", "OpenMS"),
                args = c(OpenMSArgListToOpts(settings), "-in", inFile, "-out_cm", outFile)))
}
