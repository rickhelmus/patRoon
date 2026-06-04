# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

defaultParam <- setClass("defaultParam", slots = c(dummy = "ANY"))
DEFAULT <- defaultParam()
isDefault <- function(x) inherits(x, "defaultParam")

# The param config defs are an environment, so they can be shared across objects and be read only
ParamConfigDefs <- setClass("ParamConfigDefs", contains = "environment")
setMethod("initialize", "ParamConfigDefs", function(.Object, li, ...)
{
    .Object <- callNextMethod(.Object, ...)
    .Object@.xData <- list2env(li, envir = .Object@.xData)
    lockEnvironment(.Object, bindings = TRUE)
    return(.Object)
})

# The param definitions should be instantiated at run time, eg when defaults are based on limits or external packages.
paramConfigDefsFact <- function(li)
{
    defs <- NULL
    function(...)
    {
        if (is.null(defs))
            defs <<- ParamConfigDefs(li)
        return(defs)
    }
}

prepAndVerifyParamForCall <- function(param, classN, ..., exOptsToDots = FALSE)
{
    checkmate::assertClass(param, classN)
    
    overrideArgs <- list(...)
    for (arg in names(overrideArgs))
        param[[arg]] <- overrideArgs[[arg]]
    
    ret <- as.list(param)
    
    if (exOptsToDots)
    {
        ret <- c(ret, "..." = ret$extraOpts)
        ret$extraOpts <- NULL
    }    
    
    return(ret)
}
