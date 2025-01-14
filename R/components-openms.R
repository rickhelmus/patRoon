# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include components-features.R
NULL

#' @rdname components-class
#' @export
componentsOpenMS <- setClass("componentsOpenMS", contains = "componentsFeatures")

setMethod("initialize", "componentsOpenMS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "openms", ...))

#' Componentization of adducts, isotopes etc. with OpenMS
#'
#' Uses the
#' \href{https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_MetaboliteAdductDecharger.html}{MetaboliteAdductDecharger}
#' utility (see \url{http://www.openms.de}) to generate components.
#'
#' @templateVar algo OpenMS
#' @templateVar do generate components
#' @templateVar generic generateComponents
#' @templateVar algoParam openms
#' @template algo_generator
#'
#' @details Features that show highly similar chromatographic elution profiles are grouped, and subsequently annotated
#'   with their adducts.
#'
#' @param chargeMin,chargeMax The minimum/maximum charge to consider. Corresponds to the
#'   \command{algorithm:MetaboliteFeatureDeconvolution:charge_min}/\command{algorithm:MetaboliteFeatureDeconvolution:charge_min}
#'    options.
#' @param chargeSpan The maximum charge span for a single analyte. Corresponds to
#'   \command{algorithm:MetaboliteFeatureDeconvolution:charge_span_max}.
#' @param qTry Sets how charges are determined. Corresponds to \command{algorithm:MetaboliteFeatureDeconvolution:q_try}.
#'   Valid options are \code{"heuristic"} and \code{"all"} (the \code{"feature"} option from \command{OpenMS} is
#'   currently not supported).
#' @param potentialAdducts The adducts to consider. Should be a \code{numeric} vector with probabilities for each
#'   adduct, \emph{e.g.} \code{potentialAdducts=c("[M+H]+" = 0.8, "[M+Na]+" = 0.2)}. Note that the sum of probabilities
#'   should always be \samp{1}. Furthermore, note that additions of multiple adducts should be controlled by the
#'   \code{chargeMin}/\code{chargeMax} arguments (and \emph{not} with \code{potentialAdducts}), \emph{e.g.} if
#'   \code{chargeMax=2} then both \code{[M+H]+} and \code{[2M+H]2+} may be considered. Please see the
#'   \command{algorithm:MetaboliteFeatureDeconvolution:potential_adducts} option of
#'   \href{https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_MetaboliteAdductDecharger.html}{MetaboliteAdductDecharger}
#'    for more details. If \code{NULL} then the a default is chosen with \code{\link{defaultOpenMSAdducts}} (which is
#'   \emph{not} the same as \command{OpenMS}).
#'
#'   \setsWF Should be a \code{list} where each entry specifies the potential adducts for a set. Should either be named
#'   with the sets names or follow the same order as \code{sets(fGroups)}. Example:
#'   \code{potentialAdducts=list(positive=c("[M+H]+" = 0.8, "[M+Na]+" = 0.2), negative=c("[M-H]-" = 0.8, "[M-H2O-H]-" =
#'   0.2))}
#' @param minRTOverlap,retWindow Sets feature retention tolerances when grouping features. Sets the
#'   \command{"algorithm:MetaboliteFeatureDeconvolution:retention_max_diff"} and
#'   \command{algorithm:MetaboliteFeatureDeconvolution:min_rt_overlap} options.
#' @param extraOpts Named character vector with extra command line parameters directly passed to
#'   \command{MetaboliteAdductDecharger}. Set to \code{NULL} to ignore.
#'
#' @templateVar ion TRUE
#' @templateVar minSize TRUE
#' @templateVar absMzDev \command{algorithm:MetaboliteFeatureDeconvolution:mass_max_diff} option
#' @template compon_algo-args
#'
#' @inheritParams generateComponents
#'
#' @return A \code{\link{componentsFeatures}} derived object.
#'
#' @template compon_gen-feat
#'
#' @templateVar class componentsSet
#' @template compon_gen-sets-merged
#'
#' @templateVar what \code{generateComponentsOpenMS}
#' @template uses-multiProc
#'
#' @references \insertRef{Bielow2010}{patRoon}
#'
#' @templateVar what generateComponentsOpenMS
#' @template main-rd-method
#' @export
setMethod("generateComponentsOpenMS", "featureGroups", function(fGroups, ionization = NULL, chargeMin = 1,
                                                                chargeMax = 1, chargeSpan = 3,
                                                                qTry = "heuristic",
                                                                potentialAdducts = NULL,
                                                                minRTOverlap = 0.66, retWindow = 1,
                                                                absMzDev = 0.005, minSize = 2,
                                                                relMinAdductAbundance = 0.75,
                                                                adductConflictsUsePref = TRUE,
                                                                NMConflicts = c("preferential", "mostAbundant", "mostIntense"),
                                                                prefAdducts = c("[M+H]+", "[M-H]-"),
                                                                extraOpts = NULL)
{
    ac <- checkmate::makeAssertCollection()
    ionization <- checkAndGetIonization(ionization, fGroups, add = ac)
    aapply(checkmate::assertCount, . ~ chargeMin + chargeMax + chargeSpan + minSize,
           positive = TRUE, fixed = list(add = ac))
    checkmate::assertChoice(qTry, c("heuristic", "all"))
    checkmate::assertNumeric(potentialAdducts, lower = 0, upper = 1, any.missing = FALSE, min.len = 2,
                             names = "unique", null.ok = TRUE, add = ac)
    aapply(checkmate::assertNumber, . ~ retWindow + absMzDev + relMinAdductAbundance, finite = TRUE, lower = 0,
           fixed = list(add = ac))
    checkmate::assertNumber(minRTOverlap, lower = 0, upper = 1, add = ac)
    checkmate::assertFlag(adductConflictsUsePref, add = ac)
    checkmate::assertSubset(NMConflicts, c("preferential", "mostAbundant", "mostIntense"), empty.ok = FALSE, add = ac)
    checkmate::assertCharacter(prefAdducts, min.chars = 1, any.missing = FALSE, unique = TRUE, add = ac)
    checkmate::assertList(extraOpts, any.missing = FALSE, names = "unique", null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (is.null(potentialAdducts))
        potentialAdducts <- defaultOpenMSAdducts(ionization)
    
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
                   retWindow = retWindow, absMzDev = absMzDev, extraOpts = extraOpts)
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
        writeFeatureXML(cmd$fTable, "analysis.mzML", inFile, TRUE) # NOTE: dummy input mzML file name --> not needed
        outFile <- tempfile(fileext = ".consensusXML")
        cmdMAD <- do.call(patRoon:::getOpenMSMADCommand, c(list(inFile = inFile, outFile = outFile), params))
        return(c(cmd, list(outFile = outFile), cmdMAD))
    }, logSubDir = "openms", cacheName = "componentsOpenMS")
    
    return(componentsOpenMS(fGroups = fGroups, absMzDev = absMzDev, minSize = minSize,
                            relMinAdductAbundance = relMinAdductAbundance,
                            adductConflictsUsePref = adductConflictsUsePref, NMConflicts = NMConflicts,
                            prefAdducts = prefAdducts, featureComponents = featComponents))
})

#' @rdname generateComponentsOpenMS
#' @export
setMethod("generateComponentsOpenMS", "featureGroupsSet", function(fGroups, ionization = NULL, chargeMin = 1,
                                                                   chargeMax = 1, chargeSpan = 3, qTry = "heuristic",
                                                                   potentialAdducts = NULL, ...)
{
    if (is.null(potentialAdducts))
        potentialAdducts <- rep(list(NULL), length(sets(fGroups)))
    else
    {
        checkmate::assertList(potentialAdducts, "numeric", any.missing = FALSE, all.missing = FALSE,
                              len = length(sets(fGroups)))
        checkmate::assert(
            checkmate::checkNames(names(potentialAdducts), "unnamed"),
            checkmate::checkNames(names(potentialAdducts), "unique", must.include = sets(fGroups)),
            .var.name = "potentialAdducts"
        )
        if (checkmate::testNames(names(potentialAdducts), "unnamed"))
            names(potentialAdducts) <- sets(fGroups)
        
        potentialAdducts <- potentialAdducts[sets(fGroups)]
    }
    
    generateComponentsSet(fGroups, ionization, generateComponentsOpenMS, setIonization = TRUE, chargeMin = chargeMin,
                          chargeMax = chargeMax, chargeSpan = chargeSpan, qTry = qTry, ...,
                          setArgs = lapply(potentialAdducts, function(pa) list(potentialAdducts = pa)))
})


getOpenMSMADCommand <- function(inFile, outFile, ionization, chargeMin, chargeMax, chargeSpan, qTry,
                                potentialAdducts, minRTOverlap, retWindow, absMzDev, extraOpts)
{
    boolToChr <- function(b) if (b) "true" else "false"
    
    if (ionization == "negative")
    {
        # ensure they are negative
        chargeMin <- -abs(chargeMin); chargeMax <- -abs(chargeMax)
    }
    
    settings <- list("-algorithm:MetaboliteFeatureDeconvolution:charge_min" = chargeMin,
                     "-algorithm:MetaboliteFeatureDeconvolution:charge_max" = chargeMax,
                     "-algorithm:MetaboliteFeatureDeconvolution:charge_span_max" = chargeSpan,
                     "-algorithm:MetaboliteFeatureDeconvolution:q_try" = qTry,
                     "-algorithm:MetaboliteFeatureDeconvolution:retention_max_diff" = retWindow,
                     "-algorithm:MetaboliteFeatureDeconvolution:mass_max_diff" = absMzDev,
                     "-algorithm:MetaboliteFeatureDeconvolution:min_rt_overlap" = minRTOverlap)
    
    recentMAD <- OpenMSVersionAtLeast("MetaboliteAdductDecharger", "3.0")
    
    if (!recentMAD) # negative_mode param usage was changed in OpenMS 3.0
        settings <- c(settings, list("-algorithm:MetaboliteFeatureDeconvolution:negative_mode" = boolToChr(ionization == "negative")))
    
    if (!is.null(extraOpts))
        settings <- modifyList(settings, extraOpts)
    
    settingsArgs <- OpenMSArgListToOpts(settings)
    
    if (recentMAD && ionization == "negative")
        settingsArgs <- c(settingsArgs, "-algorithm:MetaboliteFeatureDeconvolution:negative_mode")
    
    # add potential adducts later as OpenMSArgListToOpts() doesn't handle this currently...
    settingsArgs <- c(settingsArgs, "-algorithm:MetaboliteFeatureDeconvolution:potential_adducts", potentialAdducts)
    
    return(list(command = getExtDepPath("openms", "MetaboliteAdductDecharger"),
                args = c(settingsArgs, "-in", inFile, "-out_cm", outFile)))
}

#' Default adducts for OpenMS componentization
#'
#' Returns the default adducts and their probabilities when the OpenMS algorithm is used for componentization.
#' 
#' See the \code{potentialAdducts} argument of \code{\link{generateComponentsOpenMS}} for more details.
#' 
#' @param ionization The ionization polarity: either \code{"positive"} or \code{"negative"}.
#' 
#' @export
defaultOpenMSAdducts <- function(ionization)
{
    checkmate::assertChoice(ionization, c("positive", "negative"))
    if (ionization == "positive")
        return(c("[M+H]+" = 0.4,
                 "[M+Na]+" = 0.2,
                 "[M+NH4]+" = 0.2,
                 "[M+K]+" = 0.1,
                 "[M-H2O]+" = 0.1))
    # UNDONE: more for neg?
    return(c("[M-H]-" = 0.8,
             "[M-H2O-H]-" = 0.2))
}
