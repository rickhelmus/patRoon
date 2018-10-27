#' @include main.R
#' @include feature_groups-optimize.R
NULL

featureGroupsOptimizerXCMS <- setRefClass("featureGroupsOptimizerXCMS", contains = "featureGroupsOptimizer",
                                          fields = list(groupArgs = "character", retcorArgs = "character"))

featureGroupsOptimizerXCMS$methods(
    
    checkInitialParams = function(params)
    {
        for (p in c("groupArgs", "retcorArgs"))
        {
            if (!is.null(params[[p]]))
            {
                .self[[p]] <- names(params[[p]])
                params <- c(params, params[[p]])
                params[[p]] <- NULL
            }
        }
        return(params)
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
    
    getMinOptSetting = function(settingName, params)
    {
        if (settingName == "profStep")
            return(0.3)
        if (settingName == "mzwid")
            return(0.0001)
        if (settingName == "bw")
            return(0.25)
        if (settingName == "span")
            return(0.001)
        
        return(0)
    },
    
    convertOptToCallParams = function(params)
    {
        # general params
        ret <- params[names(params) %in% c("rtalign", "exportedData")]
        
        for (p in c("groupArgs", "retcorArgs"))
        {
            if (!is.null(.self[[p]]))
                ret[[p]] <- params[.self[[p]]]
        }
        
        return(ret)
    }
)
