# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include generics.R
NULL

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

param <- setClass("param", slots = c(name = "character", baseName = "character", description = "character",
                                     version = "character", date = "POSIXct", definitions = "ParamConfigDefs",
                                     data = "list"),
                  contains = "VIRTUAL")

setMethod("initialize", "param", function(.Object, name, baseName, description, version,
          definitions, ..., template = character())
{
    # construct template data with defaults
    makeDefsFromList <- function(li)
    {
        sapply(li, function(e)
        {
            return(if (is.list(e) && !is.data.frame(e) && length(e) > 0) makeDefsFromList(e) else DEFAULT)
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
    
    parsFilled <- paramListFillDefaults(object@data, object@definitions)

    # Automated validation using definitions
    for (paramName in names(object@definitions))
    {
        def <- object@definitions[[paramName]]
        value <- parsFilled[[paramName]]

        if (is.null(value) && (!is.null(def[["null.ok"]]) && def$null.ok))
            next
        if (is.null(def[["type"]]))
            next # checked elsewhere
        
        cmfunc <- switch(def$type,
            number = checkmate::assertNumber,
            character = checkmate::assertCharacter,
            numeric = checkmate::assertNumeric,
            string = checkmate::assertString,
            flag = checkmate::assertFlag,
            choice = checkmate::assertChoice,
            subset = checkmate::assertSubset,
            list = checkmate::assertList,
            data.frame = checkmate::assertDataFrame,
            count = checkmate::assertCount,
            specSimParams = assertSpecSimParams,
            IMS = assertIMSArg,
            TPStructParams = assertTPStructParams
        )
        
        if (is.null(cmfunc))
            stop(sprintf("Unknown type '%s' for parameter '%s'", def$type, paramName), call. = FALSE)
        
        do.call(cmfunc, c(list(x = value), def$typeCheckArgs, list(add = ac)))
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
