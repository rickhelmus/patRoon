#' @include main.R
#' @include components-features.R
NULL

makeCompCreateCommand <- function(inPath, fileName, onlyMS, mzRange, massWinPer, retWinPer, RThreshold, deltaMass,
                                  minInt)
{
    # UNDONE: check if julia exists? allow to configure path?
    return(list(command = "julia", args = c(system.file("misc", "runCompCreate.jl", package = "patRoon"), inPath,
                                            fileName, mzRange[1], mzRange[2], onlyMS, massWinPer, retWinPer,
                                            RThreshold, deltaMass, minInt),
                onlyMS = onlyMS, fileName = fileName))
}

compCreateMPPrepareHandler <- function(cmd)
{
    # export features to SAFD format, with some dummy values
    # NOTE: column order is important...
    featInput <- data.table(Nr = seq_len(nrow(cmd$fTable)), ScanNum = seq_len(nrow(cmd$fTable)), ScanInPeak = 10,
                            Rt = cmd$fTable$ret, RtStart = cmd$fTable$retmin, RtEnd = cmd$fTable$retmax,
                            MinInPeak = cmd$fTable$retmax - cmd$fTable$retmin, MeasMass = cmd$fTable$mz,
                            MinMass = cmd$fTable$mzmin, MaxMass = cmd$fTable$mzmax, Area = cmd$fTable$area,
                            Int = cmd$fTable$intensity, FeatPurity = 1, MediRes = 25000)
    
    featInputFile <- tempfile("compcreate_features", fileext = ".csv")
    fwrite(featInput, featInputFile)
    
    cmd$args <- c(cmd$args, featInputFile)
    return(c(cmd, list(featInputFile = featInputFile)))
}

compCreateMPFinishHandler <- function(cmd)
{
    fExt <- if (cmd$onlyMS) "_compMS1.csv" else "_comp.csv"
    results <- fread(paste0(tools::file_path_sans_ext(cmd$featInputFile), fExt))
    return(results)
}

#' @rdname components-class
#' @export
componentsCompCreate <- setClass("componentsCompCreate", contains = "componentsFeatures")

setMethod("initialize", "componentsCompCreate",
          function(.Object, ...) callNextMethod(.Object, algorithm = "openms", ...))

#'
#' @templateVar what generateComponentsCompCreate
#' @template main-rd-method
#' @export
setMethod("generateComponentsCompCreate", "featureGroups", function(fGroups, profPath = NULL, onlyMS,
                                                                    mzRange = c(0, 400), massWinPer = 0.8,
                                                                    retWinPer = 0.5, RThreshold = 0.8,
                                                                    deltaMass = 0.004, minInt = 300)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::reportAssertions(ac)
    
    if (mzRange[1] > mzRange[2])
        stop("First element of mzRange should be smaller than second.")
    
    anaInfo <- analysisInfo(fGroups)
    featTable <- featureTable(fGroups)
    
    if (!is.null(profPath))
        profPath <- rep(profPath, length.out = nrow(anaInfo))
    
    params <- list(onlyMS, mzRange, massWinPer, retWinPer, RThreshold, deltaMass, minInt)
    baseHash <- makeHash(params)
    
    printf("Annotating all features with CompCreate for %d analyses ...\n", nrow(anaInfo))
    
    anaPaths <- if (is.null(profPath)) anaInfo$path else profPath
    cmdQueue <- Map(anaInfo$analysis, anaInfo$path, f = function(ana, path)
    {
        fp <- getMzXMLAnalysisPath(ana, path)
        # UNDONE: really want to hash big files?
        hash <- makeHash(makeFileHash(fp), featTable[[ana]], params)
        
        cmd <- do.call(makeCompCreateCommand, c(list(dirname(fp), basename(fp)), params))
        
        return(c(cmd, list(inputPath = fp, hash = hash, logFile = paste0(ana, ".txt"), fTable = featTable[[ana]])))
    })
    names(cmdQueue) <- anaInfo$analysis
    
    if (length(cmdQueue) == 0)
        return(componentsCompCreate())
    
    featComponents <- executeMultiProcess(cmdQueue, patRoon:::compCreateMPFinishHandler,
                                          prepareHandler = patRoon:::compCreateMPPrepareHandler,
                                          logSubDir = "compcreate", cacheName = "componentsCompCreate")
    
    browser()
    return(componentsCompCreate(fGroups = fGroups, absMzDev = absMzDev, minSize = minSize,
                            relMinAdductAbundance = relMinAdductAbundance,
                            adductConflictsUsePref = adductConflictsUsePref, NMConflicts = NMConflicts,
                            prefAdducts = prefAdducts, featureComponents = featComponents))
})

#' @rdname generateComponentsCompCreate
#' @export
setMethod("generateComponentsCompCreate", "featureGroupsSet", function(fGroups, ...)
{
    # UNDONE
})


