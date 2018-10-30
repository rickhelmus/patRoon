#' @include main.R
#' @include feature_groups-optimize.R
NULL

featureGroupsOptimizerXCMS <- setRefClass("featureGroupsOptimizerXCMS", contains = "featureGroupsOptimizer",
                                          fields = list(groupArgs = "list", retcorArgs = "list"))

featureGroupsOptimizerXCMS$methods(

    checkInitialParams = function(params)
    {
        for (p in c("groupArgs", "retcorArgs"))
        {
            if (!is.null(params[[p]]))
            {
                .self[[p]] <- params[[p]]
                # don't save method here: both groupArgs and retcorArgs have method parameter
                params <- c(params, params[[p]][names(params[[p]]) != "method"])
                params[[p]] <- NULL
            }
        }

        return(params)
    },

    defaultParamRanges = function(params)
    {
        return(list(profStep = c(0.3, Inf),
                    mzwid = c(0.0001, Inf),
                    bw = c(0.25, Inf),
                    span = c(0.001, Inf)))
    },

    fixOptParamBounds = function(param, bounds)
    {
        if (param %in% c("extra", "missing"))
            return(round(bounds, 0))
        if (param %in% c("profStep", "minfrac"))
        {
            if (bounds[2] > 1)
            {  # 1 is max value for profStep
                bounds <- round(c(1-(diff(bounds)*0.8), 1), 2)
                printf("profStep or minfrac greater 1, decreasing to %s\n", bounds)
            }
        }

        return(bounds)
    },

    convertOptToCallParams = function(params)
    {
        # general params
        ret <- params[names(params) %in% c("rtalign", "exportedData")]

        for (p in c("groupArgs", "retcorArgs"))
        {
            if (!is.null(.self[[p]]))
            {
                pn <- names(.self[[p]])
                pn <- pn[names(pn) != "method"]
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
