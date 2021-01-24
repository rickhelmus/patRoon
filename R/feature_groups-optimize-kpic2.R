#' @include main.R
#' @include feature_groups-optimize.R
NULL

featureGroupsOptimizerKPIC2 <- setRefClass("featureGroupsOptimizerKPIC2", contains = "featureGroupsOptimizer",
                                           fields = list(groupArgs = "list", alignArgs = "list"))

featureGroupsOptimizerKPIC2$methods(
    
    flattenParams = function(params, setFields)
    {
        # combine groupArgs and alignArgs in one list for easier optimization
        
        # filter out invalid params (e.g. when user forgets to use
        # groupArgs/alignArgs for specifying params to optimize)
        params <- params[names(params) %in% c("rtalign", "exportedData", "groupArgs", "alignArgs")]
        
        for (p in c("groupArgs", "alignArgs"))
        {
            if (!is.null(params[[p]]))
            {
                if (setFields)
                    .self[[p]] <- params[[p]]
                
                # don't save method here: both groupArgs and alignArgs have method parameter
                params <- c(params, params[[p]][names(params[[p]]) != "method"])
                params[[p]] <- NULL
            }
        }
        
        return(params)
    },
    
    checkInitialParams = function(params)
    {
        # set defaults
        defs <- list(mz_tolerance = 0.01, rt_tolerance = 10, mz_weight = 0.8, rt_weight = 0.2)
        if (is.null(params[["groupArgs"]]))
            params$groupArgs <- list()
        params$groupArgs[names(defs)] <- ifelse(sapply(names(defs), function(x) is.null(params$groupArgs[[x]])), defs,
                                                params$groupArgs[names(defs)])
        
        return(flattenParams(params, TRUE))
    },
    getOptSettingRange = function(settingName, params, paramRanges) callSuper(settingName, params,
                                                                              flattenParams(paramRanges, FALSE)),
    
    fixOptParamBounds = function(param, bounds)
    {
        if (param == "frac")
        {
            if (bounds[2] > 1)
            {
                bounds <- round(c(1-(diff(bounds)*0.8), 1), 2)
                printf("frac greater 1, decreasing to %s\n", bounds)
            }
        }
        return(bounds)
    },
    
    convertOptToCallParams = function(params)
    {
        # general params
        ret <- params[names(params) %in% c("rtalign", "exportedData")]
        
        for (p in c("groupArgs", "alignArgs"))
        {
            if (!is.null(.self[[p]]))
            {
                pn <- names(.self[[p]])
                pn <- pn[pn != "method"]
                ret[[p]] <- params[pn]
                
                # re-add method
                method <- .self[[p]][["method"]]
                if (!is.null(method))
                    ret[[p]] <- c(ret[[p]], list(method = method))
            }
        }
        
        if (!is.null(ret$groupArgs[["mz_tolerance"]])) # also implies rt_tolerance
        {
            ret$groupArgs$tolerance <- c(ret$groupArgs$mz_tolerance, ret$groupArgs$rt_tolerance)
            ret$groupArgs[c("mz_tolerance", "rt_tolerance")] <- NULL
        }
        if (!is.null(ret$groupArgs[["mz_weight"]])) # also implies rt_weight
        {
            ret$groupArgs$weight <- c(ret$groupArgs$mz_weight, ret$groupArgs$rt_weight)
            ret$groupArgs[c("mz_weight", "rt_weight")] <- NULL
        }
        
        return(ret)
    }
)

generateFGroupsOptPSetKPIC2 <- function(...)
{
    givenArgs <- list(...)
    
    groupMethod <- givenArgs[["groupArgs"]][["method"]]
    if (is.null(groupMethod))
        groupMethod <- "score"
    
    alignMethod <- givenArgs[["alignArgs"]][["method"]]
    if (is.null(alignMethod))
        alignMethod <- "fftcc"
    
    groupArgs <- list(method = groupMethod)
    alignArgs <- list(method = alignMethod)
    
    groupArgs$mz_tolerance <- c(0.001, 0.01)
    groupArgs$rt_tolerance <- c(4, 12)
    groupArgs$mz_weight <- c(0.2, 1)
    # groupArgs$rt_weight <- c(0.05, 0.4) # UNDONE: optimizing both simultaneously setting both makes less sense...?
    
    # UNDONE: set align args? (only span and only if move=="loess")
    
    if (is.null(givenArgs[["alignArgs"]][["move"]]))
        alignArgs$move <- "direct"
    
    return(list(groupArgs = groupArgs, alignArgs = alignArgs))
}

getDefFGroupsOptParamRangesKPIC2 <- function()
{
    return(list(groupArgs = list(mz_tolerance = c(0.0001, Inf), rt_tolerance = c(1, Inf))))
}
