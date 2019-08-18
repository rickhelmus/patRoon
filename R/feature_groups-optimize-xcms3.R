#' @include main.R
#' @include feature_groups-optimize-xcms.R
NULL

# UNDONE: do the subset/subsetAdjust/minFraction parameters currently work for ret correction?

convertXCMS3GroupListToParam <- function(args, sGroups)
{
    method <- args$method
    args <- args[names(args) != "method"]
    args$sampleGroups <- sGroups
    
    if (method == "density")
    {
        # update to newstyle naming
        args <- setListNamesIfPresent(args,
                                      c("minfrac", "minsamp", "mzwid", "max"),
                                      c("minFraction", "minSamples", "binSize", "maxFeatures"))
        return(do.call(xcms::PeakDensityParam, args))
    }
    
    if (method == "nearest")
    {
        args <- setListNamesIfPresent(args,
                                      c("mzCheck", "rtCheck"),
                                      c("absMz", "absRt"))
        return(do.call(xcms::NearestPeaksParam, args))
    }
    
    if (method == "mzClust") # UNDONE: not really supported in XCMS base class
    {
        args <- setListNamesIfPresent(args,
                                      c("mzppm", "mzabs", "minsamp", "minfrac"),
                                      c("ppm", "absMz", "minSamples", "minFraction"))
        return(do.call(xcms::MzClustParam, args))
    }
}

convertXCMS3RetAlignListToParam <- function(args)
{
    method <- args$method
    args <- args[names(args) != "method"]
    
    if (method == "obiwarp")
    {
        args <- setListNamesIfPresent(args,
                                      c("profStep", "center", "distFunc"),
                                      c("binSize", "centerSample", "distFun"))
        return(do.call(xcms::ObiwarpParam, args))
    }
    
    if (method == "peakgroups")
    {
        # UNDONE: old style 'missing' param doesn't translate
        args <- setListNamesIfPresent(args,
                                      c("extra"),
                                      c("extraPeaks"))
        return(do.call(xcms::PeakGroupsParam, args))
    }
}

convertXCMS3GroupParamToList <- function(param)
{
    if (isXCMSClass(param, "PeakDensityParam"))
        return(list(minfrac = xcms::minFraction(param),
                    minsamp = xcms::minSamples(param),
                    bw = xcms::bw(param),
                    mzwid = xcms::binSize(param),
                    max = xcms::maxFeatures(param)))
    
    if (isXCMSClass(param, "NearestPeaksParam"))
        return(list(mzVsRTbalance = xcms::mzVsRtBalance(param),
                    mzCheck = xcms::absMz(param),
                    rtCheck = xcms::absRt(param),
                    kNN = xcms::kNN(param)))
    
    if (isXCMSClass(param, "MzClustParam"))
        return(list(mzppm = xcms::ppm(param),
                    mzabs = xcms::absMz(param),
                    minsamp = xcms::minSamples(param),
                    minfrac = xcms::minFraction(param)))
}

convertXCMS3RetAlignParamToList <- function(param)
{
    if (isXCMSClass(param, "ObiwarpParam"))
        return(list(profStep = xcms::binSize(param),
                    center = xcms::centerSample(param),
                    response = xcms::response(param),
                    distFunc = xcms::distFun(param),
                    gapInit = xcms::gapInit(param),
                    gapExtend = xcms::gapExtend(param),
                    factorDiag = xcms::factorDiag(param),
                    factorGap = xcms::factorGap(param),
                    localAlignment = xcms::localAlignment(param),
                    initPenalty = xcms::initPenalty(param)))
    
    if (isXCMSClass(param, "peakgroups"))
        return(list(extra = xcms::extraPeaks(param),
                    smooth = xcms::smooth(param),
                    span = xcms::span(param),
                    family = xcms::family(param)))
}

# Split the ranges for each entry over two lists, so we can convert it to
# lower/upper param objects.
splitXCMSOptRangeArgs <- function(args)
{
    rangeArgs <- lengths(args) == 2
    lower <- lapply(args[rangeArgs], "[[", 1)
    upper <- lapply(args[rangeArgs], "[[", 2)
    other <- args[!rangeArgs]
    return(list(c(other, lower), c(other, upper)))
}

featureGroupsOptimizerXCMS3 <- setRefClass("featureGroupsOptimizerXCMS3", contains = "featureGroupsOptimizerXCMS")

featureGroupsOptimizerXCMS3$methods(
    
    convertOptToCallParams = function(params)
    {
        params <- callSuper(params)
        ret <- list()
        
        if (!is.null(params$groupArgs[["method"]]))
            ret$groupParam <- convertXCMS3GroupListToParam(params$groupArgs,
                                                           analysisInfo(features)$group)

        if (!is.null(params$retcorArgs[["method"]]))
            ret$retAlignParam <- convertXCMS3RetAlignListToParam(params$retcorArgs)

        return(ret)
    },
    
    checkInitialParams = function(params)
    {
        # params should be a list containing:
        # - groupParam: a *list* with two xcms group param objects, specifying lower/upper bounds for grouping
        # - retAlignParam: as above, but for alignment instead of grouping
        #
        # These objects are converted to a list format the optimizer can use (or
        # actually the old-style XCMS optimizer can use). This is done by
        # merging the lower/upper values after conversion to a list object.

        mergeArgList <- function(lower, upper)
        {
            # combine lists first: this ensures unique values in either are kept.
            ret <- modifyList(lower, upper)
            
            inBoth <- intersect(names(lower), names(upper))
            ret[inBoth] <- mapply(lower[inBoth], upper[inBoth], SIMPLIFY = FALSE, FUN = c)

            return(ret)
        }
                
        ret <- list(groupArgs = mergeArgList(convertXCMS3GroupParamToList(params$groupParam[[1]]),
                                             convertXCMS3GroupParamToList(params$groupParam[[2]])),
                    retcorArgs = mergeArgList(convertXCMS3RetAlignParamToList(params$retAlignParam[[1]]),
                                              convertXCMS3RetAlignParamToList(params$retAlignParam[[2]])))
        
        return(callSuper(ret))
    }
)

generateFGroupsOptPSetXCMS3 <- function(..., features)
{
    givenArgs <- list(...)
    
    if (!is.null(givenArgs[["groupParam"]]) || !is.null(givenArgs[["retAlignParam"]]))
        stop(paste0("Overriding default settings for grouping/retention alignment currently not supported!",
                    "Please modify these values after the parameter set has been created."))
    
    # the defaults are based on the old style XCMS base class
    args <- generateFGroupsOptPSetXCMS(...)
    
    return(list(groupParam = lapply(splitXCMSOptRangeArgs(args$groupArgs), convertXCMS3GroupListToParam,
                                    sGroups = analysisInfo(features)$group),
                retAlignParam = lapply(splitXCMSOptRangeArgs(args$retcorArgs), convertXCMS3RetAlignListToParam)))
}

getDefFGroupsOptParamRangesXCMS3 <- function(features)
{
    args <- getDefFGroupsOptParamRangesXCMS()
    
    # necessary for conversion to param objects
    args$groupArgs$method <- "density"; args$retcorArgs$method <- "obiwarp"
    
    return(list(groupParam = lapply(splitXCMSOptRangeArgs(args$groupArgs), convertXCMS3GroupListToParam,
                                    sGroups = analysisInfo(features)$group),
                retAlignParam = lapply(splitXCMSOptRangeArgs(args$retcorArgs), convertXCMS3RetAlignListToParam)))
}
