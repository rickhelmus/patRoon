# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include workflow-step.R
NULL

#' Base transformation products (TP) class
#'
#' Holds information for all TPs for a set of parents.
#'
#' This class holds all generated data for transformation products for a set of parents. The class is \code{virtual} and
#' derived objects are created by \link[=generateTPs]{TP generators}.
#'
#' The TP data in objects from this class include a \code{retDir} column. These are \code{numeric} values that hint what
#' the the chromatographic retention order of a TP might be compared to its parent: a value of \samp{-1} means it will
#' elute earlier, \samp{1} it will elute later and \samp{0} that there is no significant difference or the direction is
#' unknown. These values are based on a typical reversed phase separation. When structural information is available
#' (\emph{e.g.} when \code{\link{generateTPsBioTransformer}} or \code{\link{generateTPsLibrary}} was used to generate
#' the data), the \code{retDir} values are based on calculated \code{log P} values of the parent and its TPs.
#'
#' @param TPs,x,obj,object \code{transformationProducts} object to be accessed
#'
#' @seealso The derived \code{\link{transformationProductsStructure}} class for more methods and
#'   \code{\link{generateTPs}}
#'
#' @slot parents A \code{\link{data.table}} with metadata for all parents that have TPs in this object. Use the
#'   \code{parents} method for access.
#' @slot products A \code{list} with \code{\link{data.table}} entries with TP information for each parent. Use the
#'   \code{products} method for access.
#'
#' @templateVar seli parents
#' @templateVar selOrderi names()
#' @templateVar del TRUE
#' @templateVar deli parents
#' @templateVar delj transformation products
#' @templateVar deljtype numeric index (row) or name of the TP
#' @templateVar delfwhat parent
#' @templateVar delfa the TP info table (a \code{data.table}), the parent name
#' @templateVar delfr the TP indices (rows) (specified as an \code{integer} or \code{logical} vector) or names to be removed
#' @templateVar dollarOpName parent
#' @template sub_sel_del-args
#' 
#' @param \dots For \code{delete}: passed to the function specified as \code{j}. Otherwise ignored.
#'
#' @templateVar class transformationProducts
#' @template class-hierarchy
#'
#' @export
transformationProducts <- setClass("transformationProducts",
                          slots = c(parents = "data.table", products = "list"),
                          contains = c("VIRTUAL", "workflowStep"))

#' @describeIn transformationProducts Accessor method for the \code{parents} slot of a
#'   \code{transformationProducts} class.
#' @aliases parents
#' @export
setMethod("parents", "transformationProducts", function(TPs) TPs@parents)

#' @describeIn transformationProducts Accessor method for the \code{products} slot.
#' @aliases products
#' @export
setMethod("products", "transformationProducts", function(TPs) TPs@products)

#' @describeIn transformationProducts Obtain total number of transformation products.
#' @export
setMethod("length", "transformationProducts", function(x) if (length(products(x)) == 0) 0 else sum(sapply(products(x), nrow)))

#' @describeIn transformationProducts Obtain the names of all parents in this object.
#' @export
setMethod("names", "transformationProducts", function(x) parents(x)$name)

#' @describeIn transformationProducts Show summary information for this object.
#' @export
setMethod("show", "transformationProducts", function(object)
{
    callNextMethod()
    
    printf("Parents: %s (%d total)\n", getStrListWithMax(names(object), 6, ", "), nrow(parents(object)))
    printf("Total TPs: %d\n", length(object))
})

#' @describeIn transformationProducts Subset on parents.
#' @export
setMethod("[", c("transformationProducts", "ANY", "missing", "missing"), function(x, i, ...)
{
    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, names(x))
        x <- delete(x, setdiff(names(x), i))
    }
    
    return(x)
})

#' @describeIn transformationProducts Extracts a table with TPs for a parent.
#' @export
setMethod("[[", c("transformationProducts", "ANY", "missing"), function(x, i, j)
{
    assertExtractArg(i)
    return(x@products[[i]])
})

#' @describeIn transformationProducts Extracts a table with TPs for a parent.
#' @export
setMethod("$", "transformationProducts", function(x, name)
{
    eval(substitute(x@products$NAME_ARG, list(NAME_ARG = name)))
})

# UNDONE: add more parent info?
#' @describeIn transformationProducts Returns all TP data in a table.
#' @export
setMethod("as.data.table", "transformationProducts", function(x) rbindlist(products(x), idcol = "parent"))

#' @describeIn transformationProducts Converts this object to a suspect list that can be used as input for
#'   \code{\link{screenSuspects}}.
#' @param includeParents If \code{TRUE} then parents are also included in the returned suspect list.
#' @export
setMethod("convertToSuspects", "transformationProducts", function(obj, includeParents = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(includeParents, add = ac)
    checkmate::reportAssertions(ac)

    if (length(obj) == 0)
        stop("Cannot create suspect list: no data", call. = FALSE)

    prodAll <- rbindlist(products(obj))
    
    # remove TPs that are equal for a parent that were predicted through different routes
    prodAll <- unique(prodAll, by = "name")
    
    keepCols <- c("name", "SMILES", "InChI", "InChIKey", "formula", "neutralMass", "molNeutralized")
    prodAll <- prodAll[, intersect(keepCols, names(prodAll)), with = FALSE]
    prodAll <- prepareSuspectList(prodAll, NULL, FALSE, FALSE, FALSE, FALSE)
    
    if (includeParents)
        prodAll <- rbind(parents(obj), prodAll, fill = TRUE)

    return(prodAll[])
})

setMethod("needsScreening", "transformationProducts", function(TPs) TRUE)

setMethod("linkTPsToFGroups", "transformationProducts", function(TPs, fGroups)
{
    TPNames <- as.data.table(TPs)$name
    ret <- screenInfo(fGroups)[name %in% TPNames, c("group", "name"), with = FALSE]
    setnames(ret, "name", "TP_name")
    return(ret)
})

#' @templateVar where transformationProducts
#' @templateVar what transformation product data
#' @template delete
#' @export
setMethod("delete", "transformationProducts", function(obj, i = NULL, j = NULL, ...)
{
    ac <- checkmate::makeAssertCollection()
    i <- assertDeleteArgAndToChr(i, names(obj), add = ac)
    checkmate::assert(
        checkmate::checkIntegerish(j, any.missing = FALSE, null.ok = TRUE),
        checkmate::checkString(j, min.chars = 1, null.ok = TRUE),
        checkmate::checkFunction(j, null.ok = TRUE),
        .var.name = "j"
    )
    checkmate::reportAssertions(ac)
    
    if (length(i) == 0 || (!is.null(j) && length(j) == 0))
        return(obj) # nothing to remove...
    
    # UNDONE: NULL for i and j will remove all?
    
    # i = NULL; j = vector: remove from all parents
    # i = vector; j = NULL: remove all data for specified parents
    # j = function: remove specific TPs from given parents (or all if i=NULL)
    
    if (!is.function(j))
    {
        if (is.null(j))
            obj@products <- obj@products[setdiff(names(obj), i)]
        else
        {
            obj@products[i] <- lapply(obj@products[i], function(tab)
            {
                if (is.character(j))
                    return(tab[!name %chin% j])
                inds <- j[j <= nrow(tab)]
                return(if (length(inds) > 0) tab[-inds] else tab)
            })
        }
    }
    else
    {
        obj@products[i] <- Map(obj@products[i], i, f = function(tab, par)
        {
            rm <- j(tab, par, ...)
            if (is.logical(rm))
                return(tab[!rm])
            if (is.character(rm))
                return(tab[!name %chin% rm])
            return(tab[setdiff(seq_len(nrow(tab)), rm)])
        })
    }
    
    obj@products <- pruneList(obj@products, checkZeroRows = TRUE)
    obj@parents <- obj@parents[name %in% names(obj@products), ]
    
    return(obj)
})

#' @describeIn transformationProducts Performs rule-based filtering. Useful to simplify and clean-up the data.
#'
#' @param verbose If set to \code{FALSE} then no text output is shown.
#' @param negate If \code{TRUE} then filters are performed in opposite manner.
#'
#' @templateVar getProps names(TPs)[[1]]
#' @templateVar ex likelihood=c("LIKELY","PROBABLE")
#' @template filter-properties
#' 
#' @return \code{filter} returns a filtered \code{transformationProducts} object.
#'
#' @export
setMethod("filter", "transformationProducts", function(obj, properties = NULL, verbose = TRUE, negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(properties, any.missing = FALSE, null.ok = TRUE, add = ac)
    if (!is.null(properties))
    {
        props <- if (length(obj) > 0) names(as.data.table(obj)[, -"parent"]) else NULL
        checkmate::assertNames(names(properties), type = "unique", subset.of = props, add = ac)
        checkmate::qassertr(properties, "V")
    }
    aapply(checkmate::assertFlag, . ~ verbose + negate, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0 || is.null(properties) || length(properties) == 0)
        return(obj)
    
    oldn <- length(obj)
    
    hash <- makeHash(obj, properties, negate)
    cache <- loadCacheData("filterTPs", hash)
    if (!is.null(cache))
        obj <- cache
    else
    {
        pred <- if (negate)
            function(v, p) !v %in% p
        else
            function(v, p) v %in% p
        
        obj <- delete(obj, j = function(tab, ...)
        {
            tab <- copy(tab)
            tab[, keep := TRUE]
            for (prop in names(properties))
            {
                if (is.null(tab[[prop]]))
                    stop(sprintf("Property %s not present.", prop), call. = FALSE)
                tab[keep == TRUE, keep := pred(get(prop), properties[[prop]])]
            }
            return(!tab$keep)
        })

        saveCacheData("filterTPs", obj, hash)
    }
    
    if (verbose)
    {
        newn <- length(obj)
        printf("Done! Filtered %d (%.2f%%) TPs. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)
    }
    
    return(obj)
})


#' Generation of transformation products (TPs)
#'
#' Functionality to automatically obtain transformation products for a given set of parent compounds.
#'
#' @templateVar func generateTPs
#' @templateVar what generate transformation products
#' @templateVar ex1 generateTPsBioTransformer
#' @templateVar ex2 generateTPsLogic
#' @templateVar algos biotransformer,logic,library,library_formula,cts
#' @templateVar algosSuffix BioTransformer,Logic,Library,LibraryFormula,CTS
#' @templateVar ret transformationProducts
#' @template generic-algo
#'
#' @param \dots Any parameters to be passed to the selected TP generation algorithm.
#'
#' @return A \code{\link{transformationProducts}} (derived) object containing all generated TPs.
#'
#' @seealso The derived class \code{\link{transformationProductsStructure}} for more specific methods to post-process TP
#'   data.
#'
#' @export
generateTPs <- function(algorithm, ...)
{
    checkmate::assertChoice(algorithm, c("biotransformer", "logic", "library", "library_formula", "cts"))
    
    f <- switch(algorithm,
                biotransformer = generateTPsBioTransformer,
                logic = generateTPsLogic,
                library = generateTPsLibrary,
                library_formula = generateTPsLibraryFormula,
                cts = generateTPsCTS)
    
    f(...)
}
