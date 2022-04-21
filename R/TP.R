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
#' unknown. These values are based on a typical reversed phase separation. When \code{\link{generateTPsBioTransformer}}
#' or \code{\link{generateTPsLibrary}} was used to generate the data, the \code{retDir} values are based on calculated
#' \code{log P} values of the parent and its TPs.
#'
#' @param TPs,x,object \code{transformationProducts} object to be accessed
#'
#' @seealso Derived classes \code{\link{transformationProductsBT}} and \code{\link{transformationProductsLibrary}} for
#'   specific algorithm methods and \code{\link{generateTPs}}
#'
#' @slot parents A \code{\link{data.table}} with metadata for all parents that have TPs in this object. Use the
#'   \code{parents} method for access.
#' @slot products A \code{list} with \code{\link{data.table}} entries with TP information for each parent. Use the
#'   \code{products} method for access.
#'
#' @templateVar seli parents
#' @templateVar selOrderi names()
#' @templateVar dollarOpName parent
#' @template sub_sel_del-args
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
#' @param \dots Unused.
#' @export
setMethod("[", c("transformationProducts", "ANY", "missing", "missing"), function(x, i, ...)
{
    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, names(x))
        x@products <- x@products[i]
        x@parents <- x@parents[name %in% i]
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
setMethod("convertToSuspects", "transformationProducts", function(TPs, includeParents = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(includeParents, add = ac)
    checkmate::reportAssertions(ac)

    if (length(TPs) == 0)
        stop("Cannot create suspect list: no data", call. = FALSE)

    prodAll <- rbindlist(products(TPs))
    
    # remove TPs that are equal for a parent that were predicted through different routes
    prodAll <- unique(prodAll, by = "name")
    
    keepCols <- c("name", "SMILES", "InChI", "InChIKey", "formula", "neutralMass")
    prodAll <- prodAll[, intersect(keepCols, names(prodAll)), with = FALSE]
    prodAll <- prepareSuspectList(prodAll, NULL, FALSE, FALSE)
    
    if (includeParents)
        prodAll <- rbind(parents(TPs), prodAll, fill = TRUE)

    return(prodAll)
})

#' @export
setMethod("plotGraph", "transformationProducts", function(obj, components = NULL, structuresMax = 25, prune = TRUE,
                                                          onlyCompletePaths = FALSE)
{
    # UNDONE: move to structure class
    # UNDONE: don't make name unique, but use IDs?
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(components, "componentsTPs", null.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ prune + onlyCompletePaths, fixed = list(add = ac))
    checkmate::assertCount(structuresMax, add = ac)
    checkmate::reportAssertions(ac)
    
    TPTab <- copy(as.data.table(obj))
    TPTab[, c("name_orig", "name") := .(name, make.unique(name))]
    TPTab[, parent_name := fifelse(is.na(parent_ID), parent, name[match(parent_ID, ID)]), by = "parent"]
    
    if (!is.null(components))
    {
        cmpTab <- as.data.table(components)
        TPTab <- TPTab[parent %chin% cmpTab$parent_name] # omit missing root parents
        TPTab[, present := name_orig %chin% cmpTab$TP_name]
        
        TPTab[, childPresent := FALSE]
        markChildPresent <- function(TPNames)
        {
            if (length(TPNames) == 0)
                return()
            TPTab[name %chin% TPNames, childPresent := TRUE]
            pars <- TPTab[name %chin% TPNames]$parent_name
            markChildPresent(pars[TPTab[name %chin% pars]$childPresent == FALSE])
        }
        markChildPresent(TPTab[present == TRUE]$parent_name)
        
        if (prune)
            TPTab <- TPTab[present == TRUE | childPresent == TRUE]
        if (onlyCompletePaths)
        {
            TPTab <- TPTab[present == TRUE]
            # keep removing TPs without parent until no change
            oldn <- nrow(TPTab)
            repeat
            {
                TPTab <- TPTab[parent_name == parent | parent_name %chin% name]
                newn <- nrow(TPTab)
                if (oldn == newn)
                    break
                oldn <- newn
            }
        }
    }
        
    pars <- parents(obj)
    TPTab[, parent_formula := fifelse(is.na(parent_ID),
                                      pars$formula[match(parent_name, pars$name)],
                                      formula[match(parent_name, name)])]
    TPTab[, formulaDiff := mapply(formula, parent_formula, FUN = function(f, pf)
    {
        sfl <- splitFormulaToList(subtractFormula(f, pf))
        ret <- ""
        subfl <- sfl[sfl < 0]
        if (length(subfl) > 0)
            ret <- paste0("-", formulaListToString(abs(subfl)))
        addfl <- sfl[sfl > 0]
        if (length(addfl) > 0)
            ret <- if (nzchar(ret)) paste0(ret, " +", formulaListToString(addfl)) else paste0("+", formulaListToString(addfl))
        return(ret)
    })]
    
    nodes <- data.table(id = union(TPTab$parent, TPTab$name))
    nodes[, isTP := id %chin% TPTab$name]
    nodes[isTP == TRUE, label := paste0("TP", TPTab$chem_ID[match(id, TPTab$name)])]
    nodes[isTP == FALSE, label := id]
    nodes[, group := if (.N > 1) label else "unique", by = "label"]
    nodes[, present := isTP == FALSE | TPTab$present[match(id, TPTab$name)]]
    nodes[present == TRUE, shapeProperties := list(list(list(useBorderWithImage = TRUE)))]
    nodes[present == FALSE, shapeProperties := list(list(list(useBorderWithImage = FALSE)))]
    nodes[, present := NULL]
    nodes[isTP == FALSE, level := 0]
    nodes[isTP == TRUE, level := TPTab$generation[match(id, TPTab$name)]]
    
    if (nrow(nodes) <= structuresMax)
    {
        # UNDONE: make util?
        imgf <- tempfile(fileext = ".png") # temp file is re-used
        getURIFromSMILES <- function(SMILES)
        {
            mol <- getMoleculesFromSMILES(SMILES, emptyIfFails = TRUE)[[1]]
            withr::with_png(imgf, withr::with_par(list(mar = rep(0, 4)), plot(getRCDKStructurePlot(mol, 150, 150))))
            return(knitr::image_uri(imgf))
        }
        nodes[, shape := "image"]
        nodes[, SMILES := fifelse(isTP, TPTab$SMILES[match(id, TPTab$name)], pars$SMILES[match(id, pars$name)])]
        nodes[, image := getURIFromSMILES(SMILES[1]), by = "SMILES"]
        nodes[, SMILES := NULL]
    }
    else
        nodes[, shape := "ellipse"]
    
    TPCols <- intersect(c("name", "name_lib", "SMILES", "formula", "routes", "generation", "accumulation", "production",
                          "globalAccumulation", "likelihood", "Lipinski_Violations", "Insecticide_Likeness_Violations",
                          "Post_Em_Herbicide_Likeness_Violations", "transformation", "transformation_ID", "enzyme",
                          "biosystem", "evidencedoi", "evidencedref", "sourcecomment", "datasetref"), names(TPTab))
    nodes[isTP == TRUE, title := sapply(id, function(TP)
    {
        TPTabSub <- TPTab[name == TP, TPCols, with = FALSE]
        return(paste0(names(TPTabSub), ": ", TPTabSub, collapse = "<br>"))
    })]
    
    edges <- data.table(from = TPTab$parent_name, to = TPTab$name, label = TPTab$formulaDiff)

    visNetwork::visNetwork(nodes = nodes, edges = edges) %>%
        visNetwork::visNodes(shapeProperties = list(useBorderWithImage = FALSE)) %>%
        visNetwork::visEdges(arrows = "to", font = list(align = "top", size = 12)) %>%
        visNetwork::visOptions(highlightNearest = list(enabled = TRUE, hover = TRUE, algorithm = "hierarchical"),
                               selectedBy = list(variable = "group", main = "Select duplicate TPs",
                                                 values = unique(nodes$group[nodes$group != "unique"]))) %>%
        visNetwork::visHierarchicalLayout(enabled = TRUE, sortMethod = "directed")
})

setMethod("needsScreening", "transformationProducts", function(TPs) TRUE)

setMethod("linkTPsToFGroups", "transformationProducts", function(TPs, fGroups)
{
    TPNames <- as.data.table(TPs)$name
    ret <- screenInfo(fGroups)[name %in% TPNames, c("group", "name"), with = FALSE]
    setnames(ret, "name", "TP_name")
    return(ret)
})

#' Generation of transformation products (TPs)
#'
#' Functionality to automatically obtain transformation products for a given set of parent compounds.
#'
#' @templateVar func generateTPs
#' @templateVar what generate transformation products
#' @templateVar ex1 generateTPsBioTransformer
#' @templateVar ex2 generateTPsLogic
#' @templateVar algos biotransformer,logic,library,cts
#' @templateVar algosSuffix BioTransformer,Logic,Library,CTS
#' @templateVar ret transformationProducts
#' @template generic-algo
#'
#' @param \dots Any parameters to be passed to the selected TP generation algorithm.
#'
#' @return A \code{\link{transformationProducts}} (derived) object containing all generated TPs.
#'
#' @seealso In addition, the derived classes \code{\link{transformationProductsBT}},
#'   \code{\link{transformationProductsLibrary}}, \code{\link{transformationProductsCTS}}, which are derived from
#'   \code{\link{transformationProductsStructure}}, for algorithm specific methods to post-process TP data.
#'
#' @export
generateTPs <- function(algorithm, ...)
{
    checkmate::assertChoice(algorithm, c("biotransformer", "logic", "library", "cts"))
    
    f <- switch(algorithm,
                biotransformer = generateTPsBioTransformer,
                logic = generateTPsLogic,
                library = generateTPsLibrary,
                cts = generateTPsCTS)
    
    f(...)
}
