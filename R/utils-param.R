# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

.default <- setClass(".default", slots = c(def = "ANY"))
isDefault <- function(x) inherits(x, ".default")
hasDefaultSet <- function(x) length(x@def) > 1 || !is.na(x@def) || is.null(attr(x@def, "missing")) || !attr(x@def, "missing")

setMethod("initialize", ".default", function(.Object, d, ...)
{
    if (missing(d))
    {
        d <- NA
        attr(d, "missing") <- TRUE
    }
    callNextMethod(.Object, def = d, ...)
})

setMethod("show", ".default", function(object)
{
    printf("A %s object.\n", class(object))
})

#' @exportS3Method constructive::.cstr_construct
.cstr_construct..default <- function(x, ...)
{
    args <- if (hasDefaultSet(x)) list(x@def) else list()
    code <- constructive::.cstr_apply(args, fun = ".default", ...)
    constructive::.cstr_repair_attributes(x, code, ..., ignore = "def", idiomatic_class = class(x), flag_s4 = FALSE)
}

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
        ret <- c(ret, ret$extraOpts)
        ret$extraOpts <- NULL
    }

    return(ret)
}
