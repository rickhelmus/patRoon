#' @include main.R
#' @include components-features.R
NULL

#' @export
componentsOpenMS <- setClass("componentsOpenMS", contains = "componentsFeatures")

setMethod("initialize", "componentsOpenMS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "openms", ...))

#' @rdname component-generation
#' @export
setMethod("generateComponentsOpenMS", "featureGroups", function(fGroups, ionization, chargeMin = 1,
                                                                chargeMax = 1, chargeSpan = 3,
                                                                qTry = "heuristic",
                                                                potentialAdducts = defaultOpenMSAdducts(ionization),
                                                                minRTOverlap = 0.66, retWindow = 1,
                                                                mzWindow = 0.005, minSize = 2,
                                                                relMinAdductAbundance = 0.75,
                                                                extraOpts = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertChoice(ionization, c("positive", "negative"), add = ac)
    aapply(checkmate::assertCount, . ~ chargeMin + chargeMax + chargeSpan + minSize,
           positive = TRUE, fixed = list(add = ac))
    checkmate::assertChoice(qTry, c("heuristic", "all"))
    checkmate::assertNumeric(potentialAdducts, lower = 0, upper = 1, any.missing = FALSE, min.len = 2,
                             names = "unique", add = ac)
    aapply(checkmate::assertNumber, . ~ retWindow + mzWindow + relMinAdductAbundance, finite = TRUE, lower = 0,
           fixed = list(add = ac))
    checkmate::assertNumber(minRTOverlap, lower = 0, upper = 1, add = ac)
    checkmate::assertList(extraOpts, any.missing = FALSE, names = "unique", null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!numEQ(sum(potentialAdducts), 1))
        stop("The sum of all adduct probabilities should be one.")
    
    anaInfo <- analysisInfo(fGroups)
    featTable <- featureTable(fGroups)
    
    padds <- sapply(names(potentialAdducts), as.adduct, simplify = FALSE)
    
    if (any(sapply(padds, slot, "molMult") > 1))
        stop("Molecule multiplier of potential adducts should not be >1")
    
    oadds <- sapply(padds, as.character, format = "openms")
    charges <- sapply(padds, slot, "charge")
    chsign <- if (ionization == "positive") "+" else "-"
    pa <- sprintf("%s:%s:%f", oadds, strrep(chsign, abs(charges)), potentialAdducts)
    
    params <- list(ionization = ionization, chargeMin = chargeMin, chargeMax = chargeMax, chargeSpan = chargeSpan,
                   qTry = qTry, potentialAdducts = pa, minRTOverlap = minRTOverlap,
                   retWindow = retWindow, mzWindow = mzWindow, extraOpts = extraOpts)
    baseHash <- makeHash(params, minSize, relMinAdductAbundance)
    
    printf("Annotating all features with OpenMS for %d analyses ...\n", nrow(anaInfo))
    
    cmdQueue <- lapply(seq_len(nrow(anaInfo)), function(anai)
    {
        hash <- makeHash(featTable[[anai]], baseHash)
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
        fcmp <- Filter(function(cmp) nrow(cmp) > 1 || nzchar(cmp$adduct_ion), fcmp)
        # convert to data.tables and fix adducts
        fcmp <- lapply(fcmp, function(cmp)
        {
            setDT(cmp)
            cmp[, adduct_ion := mapply(adduct_ion, charge,
                                       FUN = function(a, c) as.character(as.adduct(a, format = "openms", charge = c)))]
            cmp[, charge := NULL]
        })
        unlink(cmd[c("inFile", "outFile")]) # remove temporary files, as their size may be considerable
        return(fcmp)
    }, prepareHandler = function(cmd)
    {
        inFile <- tempfile(fileext = ".featureXML")
        writeFeatureXML(cmd$fTable, inFile, TRUE)
        outFile <- tempfile(fileext = ".consensusXML")
        cmdMAD <- do.call(patRoon:::getOpenMSMADCommand, c(list(inFile = inFile, outFile = outFile), params))
        return(c(cmd, list(outFile = outFile), cmdMAD))
    }, logSubDir = "openms", cacheName = "componentsOpenMS")
    
    return(componentsOpenMS(fGroups = fGroups, mzWindow = mzWindow, minSize = minSize,
                            relMinAdductAbundance = relMinAdductAbundance, featureComponents = featComponents))
})

#' @rdname component-generation
#' @export
setMethod("generateComponentsOpenMS", "featureGroupsSet", function(fGroups, ...)
{
    generateComponentsSet(fGroups, generateComponentsOpenMS, ...)
})


getOpenMSMADCommand <- function(inFile, outFile, ionization, chargeMin, chargeMax, chargeSpan, qTry,
                                potentialAdducts, minRTOverlap, retWindow, mzWindow, extraOpts)
{
    boolToChr <- function(b) if (b) "true" else "false"
    
    if (ionization == "negative")
    {
        # ensure they are negative
        chargeMin <- -abs(chargeMin); chargeMax <- -abs(chargeMax)
    }
    
    settings <- list("-algorithm:MetaboliteFeatureDeconvolution:negative_mode" = boolToChr(ionization == "negative"),
                     "-algorithm:MetaboliteFeatureDeconvolution:charge_min" = chargeMin,
                     "-algorithm:MetaboliteFeatureDeconvolution:charge_max" = chargeMax,
                     "-algorithm:MetaboliteFeatureDeconvolution:charge_span_max" = chargeSpan,
                     "-algorithm:MetaboliteFeatureDeconvolution:q_try" = qTry,
                     "-algorithm:MetaboliteFeatureDeconvolution:retention_max_diff" = retWindow,
                     "-algorithm:MetaboliteFeatureDeconvolution:mass_max_diff" = mzWindow,
                     "-algorithm:MetaboliteFeatureDeconvolution:min_rt_overlap" = minRTOverlap)
    
    if (!is.null(extraOpts))
        settings <- modifyList(settings, extraOpts)
    
    settingsArgs <- OpenMSArgListToOpts(settings)
    # add potential adducts later as OpenMSArgListToOpts() doesn't handle this currently...
    settingsArgs <- c(settingsArgs, "-algorithm:MetaboliteFeatureDeconvolution:potential_adducts", potentialAdducts)
    
    return(list(command = getCommandWithOptPath("MetaboliteAdductDecharger", "OpenMS"),
                args = c(settingsArgs, "-in", inFile, "-out_cm", outFile)))
}

#' @rdname component-generation
#' @export
defaultOpenMSAdducts <- function(ionization)
{
    checkmate::assertChoice(ionization, c("positive", "negative"))
    if (ionization == "positive")
        return(c("[M+H]+" = 0.5,
                 "[M+Na]+" = 0.2,
                 "[M+NH4]+" = 0.2,
                 "[M+K]+" = 0.1))
    # UNDONE: more for neg?
    return(c("[M-H]-" = 0.8,
             "[M-H2O-H]-" = 0.2))
}
