# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include feature_groups-optimize.R
NULL

convertXCMS3GroupListToParam <- function(args, method, sGroups)
{
    args$sampleGroups <- sGroups
    
    if (method == "density")
        return(do.call(xcms::PeakDensityParam, args))
    if (method == "nearest")
        return(do.call(xcms::NearestPeaksParam, args))
    if (method == "mzClust") # UNDONE: support this?
        return(do.call(xcms::MzClustParam, args))
    
    stop(paste("Unknown/Unsupported grouping method:", method))
}

convertXCMS3RetAlignListToParam <- function(args, method)
{
    if (method == "obiwarp")
        return(do.call(xcms::ObiwarpParam, args))
    if (method == "peakgroups")
        return(do.call(xcms::PeakGroupsParam, args))
    
    stop(paste("Unknown/Unsupported retention alignment method:", method))
}

featureGroupsOptimizerXCMS3 <- setRefClass("featureGroupsOptimizerXCMS3", contains = "featureGroupsOptimizer",
                                           fields = list(groupParamNames = "character", retAlignParamNames = "character"))

featureGroupsOptimizerXCMS3$methods(

    flattenParams = function(params, setNames)
    {
        # combine groupParams and retAlignParams in one list for easier optimization
        
        # filter out invalid params (e.g. when user forgets to use
        # groupParams/retAlignParams for specifying params to optimize)
        params <- params[names(params) %in% c("rtalign", "loadRawData", "groupMethod", "groupParams",
                                              "retAlignMethod", "retAlignParams")]

        if (is.null(params[["groupParams"]]))
            params$groupParams <- list()
        if (is.null(params[["retAlignParams"]]))
            params$retAlignParams <- list()
        
        if (setNames)
        {
            # save which parameters belong to where so we can re-construct the param objects later
            groupParamNames <<- if (length(params$groupParams) > 0) names(params$groupParams) else character()
            retAlignParamNames <<- if (length(params$retAlignParams) > 0) names(params$retAlignParams) else character()
        }
        
        # flatten
        params <- c(params, params[["groupParams"]], params[["retAlignParams"]])
        params[c("groupParams", "retAlignParams")] <- NULL # and remove original sublists
        
        return(params)
    },
    
    checkInitialParams = function(params)
    {
        if (is.null(params[["groupMethod"]]))
            params$groupMethod <- "density"
        if (is.null(params[["retAlignMethod"]]))
            params$retAlignMethod <- "obiwarp"
        return(flattenParams(params, TRUE))
    },
    
    getOptSettingRange = function(settingName, params, paramRanges) callSuper(settingName, params, flattenParams(paramRanges, FALSE)),

    fixDesignParam = function(param, value) if (param == "extraPeaks") round(value) else value,
    
    fixOptParamBounds = function(param, bounds)
    {
        if (param == "extraPeaks")
            return(round(bounds, 0))
        if (param %in% c("binSize", "minFraction"))
        {
            if (bounds[2] > 1)
            {
                # 1 is max value
                bounds <- round(c(1-(diff(bounds)*0.8), 1), 2)
                printf("binSize or minFraction greater 1, decreasing to %s\n", bounds)
            }
        }
        
        return(bounds)
    },
    
    convertOptToCallParams = function(params)
    {
        ret <- params[names(params) %in% c("rtalign", "loadRawData")]
        
        # convert flattened group/ret align lists back to XCMS usable param objects
        ret$groupParam <- convertXCMS3GroupListToParam(params[groupParamNames], params$groupMethod,
                                                       analysisInfo(features)$group)
        ret$retAlignParam <- convertXCMS3RetAlignListToParam(params[retAlignParamNames],
                                                             params$retAlignMethod)
        
        return(ret)
    }
)

generateFGroupsOptPSetXCMS3 <- function(...)
{
    givenArgs <- list(...)
    
    groupMethod <- givenArgs[["groupMethod"]]
    if (is.null(groupMethod))
        groupMethod <- "density"
    
    retAlignMethod <- givenArgs[["retAlignMethod"]]
    if (is.null(retAlignMethod))
        retAlignMethod <- "obiwarp"
    
    groupParams <- retAlignParams <- list()
    
    if (groupMethod == "density")
    {
        groupParams$bw <- c(22, 28)
        groupParams$binSize <- c(0.015, 0.035)
    }
    else if (groupMethod == "nearest")
    {
        # UNDONE: check if this makes sense (does this grouping method?)
        groupParams$absMz <- c(0.015, 0.035)
        groupParams$absRt <- c(5, 15)
    }
    else
        stop("Only density and nearest group methods supported.")
    
    if (is.null(givenArgs[["rtalign"]]) || givenArgs[["rtalign"]])
    {
        if (retAlignMethod == "obiwarp")
        {
            distFun <- givenArgs[["retAlignParams"]][["distFun"]]
            if (is.null(distFun))
            {
                distFun <- "cor_opt"
                # only set when not specified, otherwise it will be doubly set in generateFGroupsOptPSet()
                retAlignParams$distFun <- distFun
            }
            
            retAlignParams[c("gapInit", "gapExtend")] <- switch(distFun,
                                                                cor = list(c(0, 0.4), c(2.1, 2.7)),
                                                                cor_opt = list(c(0, 0.4), c(2.1, 2.7)),
                                                                cov = list(c(0, 0.4), c(11.4, 12)),
                                                                prd = list(c(0, 0.4), c(7.5, 8.1)),
                                                                euc = list(c(0.7, 1.1), c(1.5, 2.1)))
            retAlignParams$binSize <- c(0.7, 1)
        }
        else if (retAlignMethod == "peakgroups")
        {
            retAlignParams$extraPeaks <- c(1, 3)
            retAlignParams$span <- c(0.1, 0.3)
        }
        else
            stop("Only obiwarp and peakgroups alignment methods are supported.")
    }
    
    return(list(groupParams = groupParams, retAlignParams = retAlignParams))
}

getDefFGroupsOptParamRangesXCMS3 <- function()
{
    return(list(groupParams = list(binSize = c(0.0001, Inf), bw = c(0.25, Inf)),
                retAlignParams = list(binSize = c(0.3, Inf), span = c(0.001, Inf))))
}
