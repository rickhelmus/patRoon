param <- setClass("param", slots = c(name = "character", baseName = "ANY", description = "character",
                                     version = "character", date = "POSIXct", data = "list"),
                  contains = "VIRTUAL")

setMethod("initialize", "param", function(.Object, name, baseName, description, version, ...)
{
    obj <- callNextMethod(.Object, name = name, baseName = baseName, description = description, version = version,
                          date = Sys.time(), ...)
})

setAs("param", "list", function(from) from@data)
setMethod("as.list", "param", function(x, ...) x@data)

setMethod("show", "param", function(object)
{
    printf("Parameter class: %s\n", object@name)
    printf("Base name: %s\n", object@baseName)
    printf("Description: %s\n", object@description)
    printf("Version: %s\n", object@version)
    printf("Date created: %s\n", format(object@date))
})

setMethod("$", "param", function(x, name) x@data[[name]])
setReplaceMethod("$", "param", function(x, name, value) { x@data[[name]] <- value; return(x) })

setMethod("export", "param", function(obj, type, out)
{
    checkmate::assertChoice(type, c("R", "json", "yaml"))
    data <- as.list(obj)
    if (type == "R")
        constructive::construct_dump(data, out)
    else if (type == "json")
        jsonlite::write_json(data, out, pretty = TRUE)
    else if (type == "yaml")
        yaml::write_yaml(data, out)
})

FeatureOpenMSParam <- setClass("FeatureOpenMSParam", contains = "param")
setMethod("initialize", "FeatureOpenMSParam", function(.Object, data = list(), ...)
{
    .Object@data <- list(
        noiseThrInt = 1000,
        chromSNR = 3,
        chromFWHM = 5,
        mzPPM = defaultLim("mz", "medium_rel"),
        reEstimateMTSD = TRUE,
        traceTermCriterion = "sample_rate",
        traceTermOutliers = 5,
        minSampleRate = 0.5,
        minTraceLength = 3,
        maxTraceLength = -1,
        widthFiltering = "fixed",
        minFWHM = 1,
        maxFWHM = 30,
        traceSNRFiltering = FALSE,
        localRTRange = 10,
        localMZRange = 6.5,
        isotopeFilteringModel = "metabolites (5% RMS)",
        MZScoring13C = FALSE,
        useSmoothedInts = TRUE,
        extraOpts = NULL,
        useFFMIntensities = FALSE
    )
    .Object@data <- modifyList(.Object@data, data, keep.null = TRUE)
    
    callNextMethod(.Object, name = "FeatureOpenMSParam", baseName = "FeatureOpenMSParam",
                   description = "Parameters for OpenMS feature detection", version = "1.0", ...)
})

setValidity("FeatureOpenMSParam", function(object)
{
    # UNDONE: check if no other fields are added?
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumber(object$noiseThrInt, lower = 0, finite = TRUE, add = ac)
    checkmate::assertNumber(object$chromSNR, lower = 0, finite = TRUE, add = ac)
    checkmate::assertNumber(object$chromFWHM, lower = 0, finite = TRUE, add = ac)
    checkmate::assertNumber(object$mzPPM, lower = 0, finite = TRUE, add = ac)
    checkmate::assertFlag(object$reEstimateMTSD, add = ac)
    checkmate::assertChoice(object$traceTermCriterion, c("sample_rate", "outliers"), add = ac)
    checkmate::assertCount(object$traceTermOutliers, positive = TRUE, add = ac)
    checkmate::assertNumber(object$minSampleRate, lower = 0, upper = 1, finite = TRUE, add = ac)
    checkmate::assertCount(object$minTraceLength, positive = TRUE, add = ac)
    checkmate::assertNumber(object$maxTraceLength, lower = -1, finite = TRUE, add = ac)
    checkmate::assertChoice(object$widthFiltering, c("fixed", "adaptive"), add = ac)
    checkmate::assertNumber(object$minFWHM, lower = 0, finite = TRUE, add = ac)
    checkmate::assertNumber(object$maxFWHM, lower = 0, finite = TRUE, add = ac)
    checkmate::assertFlag(object$traceSNRFiltering, add = ac)
    checkmate::assertNumber(object$localRTRange, lower = 0, finite = TRUE, add = ac)
    checkmate::assertNumber(object$localMZRange, lower = 0, finite = TRUE, add = ac)
    checkmate::assertChoice(object$isotopeFilteringModel, c("metabolites (5% RMS)", "metabolites (10% RMS)", "lipids (5% RMS)", "lipids (10% RMS)"), add = ac)
    checkmate::assertFlag(object$MZScoring13C, add = ac)
    checkmate::assertFlag(object$useSmoothedInts, add = ac)
    checkmate::assertList(object$extraOpts, types = "character", names = "unique", null.ok = TRUE, add = ac)
    checkmate::assertFlag(object$useFFMIntensities, add = ac)
    OK <- tryCatch(checkmate::reportAssertions(ac), error = function(e) e)
    
    return(OK)
})
