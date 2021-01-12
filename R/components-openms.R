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
                                                                mzWindow = 0.005, minSize = 2,
                                                                relMinAdductAbundance = 0.75,
                                                                extraOpts = NULL)
{
    # UNDONE: all features are currently annotated (ie including not in a group), should be fine once featng is merged
    # UNDONE: keep charge column?
    # UNDONE: more parameters and proper defaults
    
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertChoice(ionization, c("positive", "negative"), add = ac)
    aapply(checkmate::assertCount, . ~ chargeMax + minSize, positive = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ mzWindow + relMinAdductAbundance, finite = TRUE, lower = 0,
           fixed = list(add = ac))
    checkmate::assertList(extraOpts, any.missing = FALSE, names = "unique", null.ok = TRUE, add = ac)
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
        # prune unassigned features
        fcmp <- Filter(function(cmp) nrow(cmp) > 1 || nzchar(cmp$adduct), fcmp)
        # convert to data.tables and fix adducts
        fcmp <- lapply(fcmp, function(cmp)
        {
            setDT(cmp)
            cmp[, adduct := mapply(adduct, charge,
                                   FUN = function(a, c) as.character(as.adduct(a, format = "openms", charge = c)))]
        })
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
    
    return(componentsOpenMS(fGroups = fGroups, mzWindow = mzWindow, minSize = minSize,
                            relMinAdductAbundance = relMinAdductAbundance, featureComponents = featComponents))
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
