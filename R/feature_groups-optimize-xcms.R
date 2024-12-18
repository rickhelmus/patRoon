# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include feature_groups-optimize.R
NULL

featureGroupsOptimizerXCMS <- setRefClass("featureGroupsOptimizerXCMS", contains = "featureGroupsOptimizer",
                                          fields = list(groupArgs = "list", retcorArgs = "list"))

featureGroupsOptimizerXCMS$methods(

    flattenParams = function(params, setFields)
    {
        # combine groupArgs and retcorArgs in one list for easier optimization

        # filter out invalid params (e.g. when user forgets to use
        # groupArgs/retcorArgs for specifying params to optimize)
        params <- params[names(params) %in% c("rtalign", "loadRawData", "groupArgs", "retcorArgs")]
        
        for (p in c("groupArgs", "retcorArgs"))
        {
            if (!is.null(params[[p]]))
            {
                if (setFields)
                    .self[[p]] <- params[[p]]
                
                # don't save method here: both groupArgs and retcorArgs have method parameter
                params <- c(params, params[[p]][names(params[[p]]) != "method"])
                params[[p]] <- NULL
            }
        }
        
        return(params)
    },
    
    checkInitialParams = function(params) flattenParams(params, TRUE),
    getOptSettingRange = function(settingName, params, paramRanges) callSuper(settingName, params, flattenParams(paramRanges, FALSE)),
    
    fixDesignParam = function(param, value) if (param %in% c("extra", "missing")) round(value) else value,
    
    fixOptParamBounds = function(param, bounds)
    {
        if (param %in% c("extra", "missing"))
            return(round(bounds, 0))
        if (param %in% c("profStep", "minfrac"))
        {
            if (bounds[2] > 1)
            {
                # 1 is max value for profStep
                bounds <- round(c(1-(diff(bounds)*0.8), 1), 2)
                printf("profStep or minfrac greater 1, decreasing to %s\n", bounds)
            }
        }

        return(bounds)
    },

    convertOptToCallParams = function(params)
    {
        # general params
        ret <- params[names(params) %in% c("rtalign", "loadRawData")]

        for (p in c("groupArgs", "retcorArgs"))
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

        return(ret)
    }
)

generateFGroupsOptPSetXCMS <- function(...)
{
    givenArgs <- list(...)

    groupMethod <- givenArgs[["groupArgs"]][["method"]]
    if (is.null(groupMethod))
        groupMethod <- "density"

    retcorMethod <- givenArgs[["retcorArgs"]][["method"]]
    if (is.null(retcorMethod))
        retcorMethod <- "obiwarp"

    groupArgs <- list(method = groupMethod)
    retcorArgs <- list(method = retcorMethod)

    if (groupMethod == "density")
    {
        groupArgs$bw <- c(22, 28)
        groupArgs$mzwid <- c(0.015, 0.035)
    }
    else if (groupMethod == "nearest")
    {
        # UNDONE: check if this makes sense (does this grouping method?)
        groupArgs$mzCheck <- c(0.015, 0.035)
        groupArgs$rtCheck <- c(5, 15)
    }
    else
        stop("Only density and nearest group methods supported.")

    if (is.null(givenArgs[["rtalign"]]) || givenArgs[["rtalign"]])
    {
        if (retcorMethod == "obiwarp")
        {
            distFunc <- givenArgs[["retcorArgs"]][["distFunc"]]
            if (is.null(distFunc))
            {
                distFunc <- "cor_opt"
                # only set when not specified, otherwise it will be doubly set in generateFGroupsOptPSet()
                retcorArgs$distFunc <- distFunc
            }

            retcorArgs[c("gapInit", "gapExtend")] <- switch(distFunc,
                                                            cor = list(c(0, 0.4), c(2.1, 2.7)),
                                                            cor_opt = list(c(0, 0.4), c(2.1, 2.7)),
                                                            cov = list(c(0, 0.4), c(11.4, 12)),
                                                            prd = list(c(0, 0.4), c(7.5, 8.1)),
                                                            euc = list(c(0.7, 1.1), c(1.5, 2.1)))
            retcorArgs$profStep <- c(0.7, 1)
        }
        else if (retcorMethod == "peakgroups")
        {
            retcorArgs$missing <- c(1, 3)
            retcorArgs$extra <- c(1, 3)
            retcorArgs$span <- c(0.1, 0.3)
        }
        else
            stop("Only obiwarp and peakgroups alignment methods are supported.")
    }

    return(list(groupArgs = groupArgs, retcorArgs = retcorArgs))
}

getDefFGroupsOptParamRangesXCMS <- function()
{
    return(list(groupArgs = list(mzwid = c(0.0001, Inf), bw = c(0.25, Inf)),
                retcorArgs = list(profStep = c(0.3, Inf), span = c(0.001, Inf))))
}
