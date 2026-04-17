defaultParam <- setClass("defaultParam", slots = c(dummy = "ANY"))
DEFAULT <- defaultParam()
isDefault <- function(x) inherits(x, "defaultParam")

paramListFillDefaults <- function(paramList, definitions)
{
    fillDefaults <- function(li, defs)
    {
        Map(names(li), li, f = function(name, value)
        {
            if (isDefault(value))
                return(defs[[name]])
            else if (is.list(value) && !is.data.frame(value))
                return(fillDefaults(value, defs[[name]]))
            else
                return(value)
        })
    }
    return(fillDefaults(paramList, sapply(definitions, "[[", "default", simplify = FALSE)))
}

setMethod("show", "defaultParam", function(object)
{
    printf("A %s object.\n", class(object))
})

param <- setClass("param", slots = c(name = "character", baseName = "ANY", description = "character",
                                     version = "character", date = "POSIXct", definitions = "list", data = "list"),
                  contains = "VIRTUAL")

setMethod("initialize", "param", function(.Object, name, baseName, description, version,
          definitions, ...)
{
    # construct template data with defaults
    makeDefsFromList <- function(li)
    {
        sapply(li, function(e)
        {
            return(if (is.list(e) && !is.data.frame(e)) makeDefsFromList(e) else DEFAULT)
        }, simplify = FALSE)
    }
    data <- sapply(definitions, "[[", "default", simplify = FALSE)
    data <- makeDefsFromList(data)
    data <- modifyList(data, list(...), keep.null = TRUE)
 
    callNextMethod(.Object, name = name, baseName = baseName, description = description, version = version,
                   date = Sys.time(), definitions = definitions, data = data)
})

setValidity("param", function(object)
{
    ac <- checkmate::makeAssertCollection()
    
    # Automated validation using definitions
    for (param_name in names(object@definitions))
    {
        def <- object@definitions[[param_name]]
        value <- object@data[[param_name]]
        
        if (isDefault(value))
            next
        if (is.null(value) && def$type == "list" && def$null.ok)
            next
        
        # Validate based on parameter type
        if (def$type == "number")
        {
            args <- list(x = value, finite = TRUE)
            if (!is.null(def$positive) && def$positive)
                args$lower <- 0
            if (!is.null(def$lower))
                args$lower <- def$lower
            if (!is.null(def$upper))
                args$upper <- def$upper
            if (!is.null(def$finite))
                args$finite <- def$finite
            do.call(checkmate::assertNumber, c(args, list(add = ac)))
        }
        else if (def$type == "flag")
            checkmate::assertFlag(value, add = ac)
        else if (def$type == "choice")
            checkmate::assertChoice(value, def$choices, add = ac)
        else if (def$type == "list")
        {
            list_args <- list(x = value, types = def$types, names = def$names, null.ok = def$null.ok)
            do.call(checkmate::assertList, c(list_args, list(add = ac)))
        }    
        else if (def$type == "count")
        {
            count_args <- list(x = value, finite = TRUE)
            if (!is.null(def$positive) && def$positive) count_args$positive <- TRUE
            do.call(checkmate::assertCount, c(count_args, list(add = ac)))
        }
    }
    
    OK <- tryCatch(checkmate::reportAssertions(ac), error = function(e) e)
    
    return(OK)
})

setMethod("as.list", "param", function(x, fill = TRUE) if (fill) paramListFillDefaults(x@data, x@definitions) else x@data)
setAs("param", "list", function(from) as.list(from))

setMethod("show", "param", function(object)
{
    printf("Parameter class: %s\n", object@name)
    printf("Base name: %s\n", object@baseName)
    printf("Description: %s\n", object@description)
    printf("Version: %s\n", object@version)
    printf("Date created: %s\n", format(object@date))
    printf("Parameters (%d in total)\n", length(object@definitions))
    # UNDONE: mark defaults
    print(data.table(parameter = names(object@definitions), value = sapply(as.list(object), as.character),
                     description = sapply(object@definitions, "[[", "description")))
})

setMethod("[[", c("param", "ANY", "missing"), function(x, i, j)
{
    return(x@data[[i]])
})

setReplaceMethod("[[", c("param", "ANY", "missing"), function(x, i, j, value)
{
    checkmate::assertChoice(i, names(x@definitions))
    x@data[[i]] <- value
    validObject(x)
    return(x)
})

setMethod("$", "param", function(x, name)
{
    eval(substitute(x@data$NAME_ARG, list(NAME_ARG = name)))
})

setReplaceMethod("$", "param", function(x, name, value)
{
    checkmate::assertChoice(name, names(x@definitions))
    eval(substitute(x@data$NAME_ARG <- VALUE_ARG, list(NAME_ARG = name, VALUE_ARG = value)))
    validObject(x)
    return(x)
})

setMethod("names", "param", function(x) names(x@definitions))
setMethod("length", "param", function(x) length(x@definitions))

setMethod("export", "param", function(obj, type, out)
{
    checkmate::assertChoice(type, c("R", "json", "yaml"))
    exp <- list(
        name = obj@name,
        classType = class(obj),
        version = obj@version,
        date = obj@date,
        data = obj@data
    )
    
    if (type != "R")
    {
        # prep defaults
        exp$data <- rapply(exp$data, classes = "defaultParam", how = "replace", function(x)
        {
            if (type == "json") jsonlite::unbox("__default") else "__default"
        })
    }
    
    if (type == "R")
        constructive::construct_dump(list(param = exp), out, data = list(DEFAULT = DEFAULT))
    else if (type == "json")
    {
        # unbox scalars
        n <- c("name", "classType", "version", "date")
        exp[n] <- lapply(exp[n], jsonlite::unbox)
        exp$data <- Map(names(exp$data), exp$data, f = function(name, value)
        {
            if (obj@definitions[[name]]$type %in% c("number", "string", "flag", "choice", "count"))
                value <- jsonlite::unbox(value)
            return(value)
        })
        jsonlite::write_json(exp, out, pretty = TRUE)
    }
    else if (type == "yaml")
        yaml::write_yaml(exp, out) # UNDONE: keep this? Format is probably only suitable for simple params
})

importParam <- function(type, inFile)
{
    checkmate::assertChoice(type, c("R", "json", "yaml"))
    imp <- if (type == "R")
        eval(parse(inFile))
    else if (type == "json")
        jsonlite::read_json(inFile, simplifyVector = TRUE)
    else # if (type == "yaml")
        yaml::read_yaml(inFile)
    
    imp <- rapply(imp, classes = "character", how = "replace", f = function(x)
    {
        if (x == "__default")
            x <- DEFAULT
        return(x)
    })
    
    ret <- new(imp$class)
    ret@version <- imp$version
    ret@date <- as.POSIXct(imp$date)
    ret@data <- imp$data
    
    # UNDONE: handle version updates
    
    validObject(ret)
    return(ret)
}

FeaturesOpenMSParam <- setClass("FeaturesOpenMSParam", contains = "param")
setMethod("initialize", "FeaturesOpenMSParam", function(.Object, ...)
{
    # UNDONE: move defs to global list, so it is RO and can be shared between objects
    defs <- list(
        noiseThrInt = list(
            default = 1000,
            description = "Intensity threshold",
            type = "number",
            positive = TRUE,
            finite = TRUE
        ),
        chromSNR = list(
            default = 3,
            description = "Chromatographic signal-to-noise ratio threshold",
            type = "number",
            positive = TRUE,
            finite = TRUE
        ),
        chromFWHM = list(
            default = 5,
            description = "Chromatographic full width at half maximum",
            type = "number",
            positive = TRUE,
            finite = TRUE
        ),
        mzPPM = list(
            default = defaultLim("mz", "medium_rel"),
            description = "Mass-to-charge ratio tolerance in parts per million",
            type = "number",
            positive = TRUE,
            finite = TRUE
        ),
        reEstimateMTSD = list(
            default = TRUE,
            description = "Re-estimate mass trace signal-to-noise ratio",
            type = "flag"
        ),
        traceTermCriterion = list(
            default = "sample_rate",
            description = "Criterion for trace termination",
            type = "choice",
            choices = c("sample_rate", "outliers")
        ),
        traceTermOutliers = list(
            default = 5,
            description = "Number of outliers for trace termination",
            type = "number",
            positive = TRUE,
            finite = TRUE
        ),
        minSampleRate = list(
            default = 0.5,
            description = "Minimum sample rate for feature detection",
            type = "number",
            lower = 0,
            upper = 1,
            finite = TRUE
        ),
        minTraceLength = list(
            default = 3,
            description = "Minimum trace length",
            type = "number",
            positive = TRUE,
            finite = TRUE
        ),
        maxTraceLength = list(
            default = -1,
            description = "Maximum trace length (-1 for unlimited)",
            type = "number",
            lower = -1,
            finite = TRUE
        ),
        widthFiltering = list(
            default = "fixed",
            description = "Width filtering method",
            type = "choice",
            choices = c("fixed", "off", "auto")
        ),
        minFWHM = list(
            default = 1,
            description = "Minimum full width at half maximum",
            type = "number",
            positive = TRUE,
            finite = TRUE
        ),
        maxFWHM = list(
            default = 30,
            description = "Maximum full width at half maximum",
            type = "number",
            positive = TRUE,
            finite = TRUE
        ),
        traceSNRFiltering = list(
            default = FALSE,
            description = "Enable trace signal-to-noise ratio filtering",
            type = "flag"
        ),
        localRTRange = list(
            default = 10,
            description = "Local retention time range for feature detection",
            type = "number",
            positive = TRUE,
            finite = TRUE
        ),
        localMZRange = list(
            default = 6.5,
            description = "Local mass-to-charge range for feature detection",
            type = "number",
            positive = TRUE,
            finite = TRUE
        ),
        isotopeFilteringModel = list(
            default = "metabolites (5% RMS)",
            description = "Isotope filtering model",
            type = "choice",
            choices = c("metabolites (5% RMS)", "metabolites (10% RMS)", "lipids (5% RMS)", "lipids (10% RMS)")
        ),
        MZScoring13C = list(
            default = FALSE,
            description = "Enable 13C scoring for mass-to-charge",
            type = "flag"
        ),
        useSmoothedInts = list(
            default = TRUE,
            description = "Use smoothed intensities",
            type = "flag"
        ),
        extraOpts = list(
            default = NULL,
            description = "Extra options as character list",
            type = "list",
            types = "character",
            names = "unique",
            null.ok = TRUE
        ),
        useFFMIntensities = list(
            default = FALSE,
            description = "Use feature finding MS intensities",
            type = "flag"
        )
    )

    callNextMethod(.Object, name = "FeaturesOpenMSParam", baseName = "FeaturesOpenMSParam",
                   description = "Parameters for OpenMS feature detection", version = "1.0",
                   definitions = defs, ...)
})


FeaturesPiekParam <- setClass("FeaturesPiekParam", contains = "param")
setMethod("initialize", "FeaturesPiekParam", function(.Object, ...)
{
    defs <- list(
        genEICParams = list(
            default = getPiekEICParams(),
            description = "Parameters for EIC generation (use getPiekEICParams())",
            type = "list",
            null.ok = FALSE
        ),
        peakParams = list(
            default = getDefPeakParams("chrom", "piek"),
            description = "Parameters for peak detection (use getDefPeakParams())",
            type = "list",
            null.ok = FALSE
        ),
        IMS = list(
            default = FALSE,
            description = "Use ion mobility separation (IMS)",
            type = "flag"
        ),
        assignMethod = list(
            default = "basepeak",
            description = "Method to assign m/z/mobility across EIC datapoints",
            type = "choice",
            choices = c("basepeak", "weighted.mean")
        ),
        assignRTWindow = list(
            default = defaultLim("retention", "very_narrow"),
            description = "Retention time window (+/- s) for assignment",
            type = "number",
            lower = 0,
            finite = TRUE
        ),
        rtWindowDup = list(
            default = defaultLim("retention", "narrow"),
            description = "RT window for duplicate feature detection",
            type = "number",
            lower = 0,
            finite = TRUE
        ),
        mzWindowDup = list(
            default = defaultLim("mz", "medium"),
            description = "m/z window for duplicate feature detection",
            type = "number",
            lower = 0,
            finite = TRUE
        ),
        mobWindowDup = list(
            default = defaultLim("mobility", "medium"),
            description = "Mobility window for duplicate feature detection",
            type = "number",
            lower = 0,
            finite = TRUE
        ),
        minPeakOverlapDup = list(
            default = 0.25,
            description = "Minimum retention time overlap (fraction) to consider duplicates",
            type = "number",
            lower = 0,
            upper = 1,
            finite = TRUE
        ),
        minIntensityIMS = list(
            default = 25,
            description = "Minimum intensity for IMS datapoints",
            type = "number",
            positive = TRUE,
            finite = TRUE
        ),
        EICBatchSize = list(
            default = Inf,
            description = "Number of EICs processed per batch",
            type = "number",
            positive = TRUE,
            finite = FALSE
        ),
        keepDups = list(
            default = FALSE,
            description = "Keep duplicate / non-centered features",
            type = "flag"
        ),
        verbose = list(
            default = TRUE,
            description = "Verbose output",
            type = "flag"
        )
    )

    callNextMethod(.Object, name = "FeaturesPiekParam", baseName = "FeaturesPiekParam",
                   description = "Parameters for piek feature detection", version = "1.0",
                   definitions = defs, ...)
})

setValidity("FeaturesPiekParam", function(object)
{
    parsFilled <- paramListFillDefaults(object@data, object@definitions)
    
    ac <- checkmate::makeAssertCollection()
    assertPiekGenEICParams(parsFilled$genEICParams, .var.name = "genEICParams", add = ac)
    assertFindPeakParams(parsFilled$peakParams, .var.name = "peakParams", add = ac)
    
    OK <- tryCatch(checkmate::reportAssertions(ac), error = function(e) e)
    return(OK)
})
