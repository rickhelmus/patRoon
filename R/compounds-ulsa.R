#' @include main.R
#' @include compounds.R
#' @include feature_groups-set.R
NULL

makeULSAInput <- function(fGroups, MSPeakLists, adduct, absMzDev, out)
{
    # NOTE: ULSA calculates the mass tolerance as mass +/- (MaxMass-MinMass)
    # To emulate our own tolerance, simply make sure MaxMass-MinMass == absMzDev/2
    
    MSPLTab <- as.data.table(MSPeakLists)
    MSPLTab <- MSPLTab[group %in% MSPLTab[type == "MSMS"]$group]
    gNames <- unique(MSPLTab$group)
    fGroupsSub <- fGroups[, gNames]
    gInfo <- groupInfo(fGroupsSub)
    
    neutralMasses <- if (is.null(adduct))
        annotations(fGroupsSub)$neutralMass
    else
        calculateMasses(gInfo$mzs, adduct, "neutral")
    
    MS2Fields <- MSPLTab[type == "MSMS",
                         .(paste0("[", paste0(mz, collapse = ", "), "]"),
                           paste0("[", paste0(intensity, collapse = ", "), "]")),
                         by = "group"]
    
    inputTab <- data.table(Nr = seq_along(gNames), ScanNum = seq_along(gNames), ScanInPeak = 10,
                           Rt = gInfo$rts, RtStart = gInfo$rts - 1, RtEnd = gInfo$rts + 1, SecInPeak = 10,
                           MeasMass = gInfo$mzs, MinMass = gInfo$mzs - (absMzDev/2), MaxMass = gInfo$mzs + (absMzDev/2),
                           Area = 1E6, Intensity = 2.5E5, FeatPurity = 1, MediRes = 25000, Parent = 1,
                           AccuMass = neutralMasses, MS1Comp = 0, MS1CompInt = 0,
                           MS2Comp = MS2Fields[["V1"]], MS2CompInt = MS2Fields[["V2"]])
    fwrite(inputTab, out)
}

#' @templateVar what generateCompoundsULSA
#' @template main-rd-method
#' @export
setMethod("generateCompoundsULSA", "featureGroups", function(fGroups, MSPeakLists, adduct = NULL, source = "ESI",
                                                             weightF = 1, absMzDev = 0.002)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::reportAssertions(ac)

    if (length(fGroups) == 0)
        return(compounds(algorithm = "ulsa"))
    
    gCount <- length(fGroups)
    gNames <- names(fGroups)
    
    checkAdductForPos <- if (!is.null(adduct)) adduct else as.adduct(annotations(fGroups)$adduct[1])
    ionization <- if (checkAdductForPos@charge > 0) "POSITIVE" else "NEGATIVE"
    weightF <- rep(weightF, length.out = 7)
    
    baseHash <- makeHash(adduct, ionization, source, weightF, absMzDev)
    
    printf("Processing %d feature groups with ULSA...\n", gCount)
    
    inputFile <- tempfile("ulsa_input", fileext = ".csv")
    makeULSAInput(fGroups, MSPeakLists, adduct, absMzDev, inputFile)
    executeCommand("julia", c(system.file("misc", "runULSA.jl", package = "patRoon"), ionization,
                              source, as.character(weightF), inputFile))
    browser()
    cmdQueue <- pruneList(sapply(gNames, function(grp)
    {
        if (is.null(MSPeakLists[[grp]]) || is.null(MSPeakLists[[grp]][["MS"]]) || is.null(MSPeakLists[[grp]][["MSMS"]]))
            return(NULL)
        
        pmz <- MSPeakLists[[grp]]$MS[precursor == TRUE]$mz
        pl <- MSPeakLists[[grp]]$MSMS
        
        hash <- makeHash(pmz, pl, baseHash)
        
        
        return(list(command = "julia", args = c(system.file("misc", "runULSA.jl", package = "patRoon"), ionization,
                                                source, as.character(weightF)),
                    hash = hash, logFile = paste0(ana, ".txt"), precursorMZ = pmz, MSMSPL = pl))
    }, simplify = FALSE))

    results <- executeMultiProcess(cmdQueue, patRoon:::ULSAMPFinishHandler,
                                   prepareHandler = patRoon:::ULSAMPPrepareHandler,
                                   logSubDir = "ulsa", cacheName = "compoundsULSA")
    
    ngrp <- length(results)
    printf("Assigned %d compounds to %d feature groups (%.2f%%).\n", sum(unlist(lapply(results, function(r) nrow(r$comptab)))),
           ngrp, if (gCount == 0) 0 else ngrp * 100 / gCount)

    return(compounds(groupAnnotations = lapply(results, "[[", "comptab"), scoreTypes = "score",
                     scoreRanges = lapply(results, "[[", "scRanges"),
                     algorithm = "ulsa"))
})

#' @template featAnnSets-gen_args
#' @rdname generateCompoundsULSA
#' @export
setMethod("generateCompoundsULSA", "featureGroupsSet", function(fGroups, MSPeakLists, ...,
                                                                setThreshold = 0, setThresholdAnn = 0)
{
    generateCompoundsSet(fGroups, MSPeakLists, adduct = NULL, generateCompoundsULSA, ...,
                         setThreshold = setThreshold, setThresholdAnn = setThresholdAnn)
})
