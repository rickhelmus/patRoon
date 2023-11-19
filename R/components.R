#' @include main.R
#' @include feature_groups.R
#' @include workflow-step.R
#' @include utils-components.R
NULL

#' Component class
#'
#' Contains data for feature groups that are related in some way. These
#' \emph{components} commonly include adducts, isotopes and homologues.
#'
#' \code{components} objects are obtained from \code{\link{generateComponents}}.
#'
#' @slot components List of all components in this object. Use the
#'   \code{componentTable} method for access.
#' @slot componentInfo A \code{\link{data.table}} containing general information
#'   for each component. Use the \code{componentInfo} method for access.
#'
#' @param obj,object,x The \code{component} object.
#' @param index The index of the component. Can be a numeric index or a
#'   character with its name.
#' @param \dots For \code{delete}: passed to the function specified as \code{j}.
#' 
#'   For \code{plotChroms}: Further (optional) arguments passed to the
#'   \code{plotChroms} method for the \code{\link{featureGroups}} class. Note that
#'   the \code{colourBy}, \code{showPeakArea}, \code{showFGroupRect} and
#'   \code{topMost} arguments cannot be set as these are set by this method.
#'
#'   For \code{plotSpectrum}: Further arguments passed to
#'   \code{\link[graphics]{plot}}.
#'
#'   For \code{consensus}: \code{components} objects that should be used to
#'   generate the consensus.
#'   
#'   \setsPassedArgs1{components}
#'
#' @templateVar seli components
#' @templateVar selOrderi names()
#' @templateVar selj feature groups
#' @templateVar selOrderj groupNames()
#' @templateVar optionalj TRUE
#' @templateVar del TRUE
#' @templateVar deli components
#' @templateVar delj feature groups
#' @templateVar deljtype numeric index/logical (relative to component) or character
#' @templateVar delfwhat component
#' @templateVar delfa the component (a \code{data.table}), the component name
#' @templateVar delfr the feature groups to be removed (same format as \code{j})
#' @templateVar dollarOpName component
#' @template sub_sel_del-args
#'
#' @templateVar class components
#' @template class-hierarchy
#'
#' @seealso \code{\link{generateComponents}}
#'
#' @export
components <- setClass("components",
                       slots = c(components = "list", componentInfo = "data.table"),
                       contains = "workflowStep")

setMethod("initialize", "components", function(.Object, ...)
{
    .Object <- callNextMethod(.Object, ...)
    .Object@components <- makeEmptyListNamed(.Object@components)
    return(.Object)
})

setMethod("collapseComponents", "components", function(obj) obj)

setMethod("groupNamesResults", "components", function(obj) groupNames(obj))

#' @describeIn components Accessor method for the \code{components} slot of a
#'   \code{components} class. Each component is stored as a
#'   \code{\link{data.table}}.
#' @aliases componentTable
#' @export
setMethod("componentTable", "components", function(obj) obj@components)

#' @describeIn components Accessor method for the \code{componentInfo} slot of a
#'   \code{components} class.
#' @aliases componentInfo
#' @export
setMethod("componentInfo", "components", function(obj) obj@componentInfo)

#' @templateVar class components
#' @templateVar what feature groups
#' @template strmethod
#' @export
setMethod("groupNames", "components", function(obj) unique(unlist(sapply(obj@components, "[[", "group", simplify = FALSE), use.names = FALSE)))

#' @describeIn components Obtain total number of components.
#' @export
setMethod("length", "components", function(x) length(x@components))

#' @describeIn components Obtain the names of all components.
#' @export
setMethod("names", "components", function(x) names(x@components))

#' @describeIn components Show summary information for this object.
#' @export
setMethod("show", "components", function(object)
{
    callNextMethod()

    printf("Components: %s (%d total)\n", getStrListWithMax(names(object), 6, ", "), length(object))

    if (length(object@components) == 0)
        rowCounts <- 0
    else
        rowCounts <- sapply(object@components, nrow)

    printf("Number of unique feature groups in this object: %d\n", length(groupNames(object)))
    printf("Size of components: %d (total), %.1f (mean), %d - %d (min - max)\n", sum(rowCounts), mean(rowCounts),
           min(rowCounts), max(rowCounts))
})

#' @describeIn components Subset on components/feature groups.
#' @export
setMethod("[", c("components", "ANY", "ANY", "missing"), function(x, i, j, ...)
{
    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, names(x))
        x <- delete(x, setdiff(names(x), i))
    }
    
    if (!missing(j))
    {
        j <- assertSubsetArgAndToChr(j, groupNames(x))
        x <- delete(x, j = setdiff(groupNames(x), j))
    }
    
    return(x)
})

#' @describeIn components Extracts a component table, optionally filtered by a feature group.
#' @export
setMethod("[[", c("components", "ANY", "ANY"), function(x, i, j)
{
    assertExtractArg(i)

    if (missing(j))
        return(x@components[[i]])

    assertExtractArg(j)
    if (!is.character(j))
        j <- groupNames(x)[j]
    return(x@components[[i]][group == j])
})

#' @describeIn components Extracts a component table by component name.
#' @export
setMethod("$", "components", function(x, name)
{
    eval(substitute(x@components$NAME_ARG, list(NAME_ARG = name)))
})

#' @templateVar where components
#' @templateVar what (parts of) components
#' @template delete
#' @export
setMethod("delete", "components", function(obj, i = NULL, j = NULL, ...)
{
    gNames <- groupNames(obj)
    
    ac <- checkmate::makeAssertCollection()
    i <- assertDeleteArgAndToChr(i, names(obj), add = ac)
    checkmate::assert(
        checkmate::checkFunction(j, null.ok = TRUE),
        checkmate::checkIntegerish(j, any.missing = FALSE, null.ok = TRUE),
        checkmate::checkCharacter(j, any.missing = FALSE),
        .var.name = "j"
    )
    checkmate::reportAssertions(ac)
    
    if (length(i) > 0 && (is.null(j) || length(j) > 0))
    {
        # i = NULL; j = vector: remove from all components
        # i = vector; j = NULL: remove specified components
        # j = function: remove specific feature groups from given components (or all components if i=NULL)
        
        if (!is.function(j))
        {
            if (is.null(j))
                obj@components <- obj@components[setdiff(names(obj), i)]
            else
            {
                obj@components[i] <- lapply(obj@components[i], function(ct)
                {
                    if (is.character(j))
                        return(ct[!group %chin% j])
                    inds <- j[j <= nrow(ct)]
                    return(if (length(inds) > 0) ct[-inds] else ct)
                })
            }
        }
        else
        {
            obj@components[i] <- Map(obj@components[i], i, f = function(ct, cmp)
            {
                rm <- j(ct, cmp, ...)
                if (is.logical(rm))
                    return(ct[!rm])
                else if (is.character(rm))
                    return(ct[!group %chin% rm])
                return(ct[setdiff(seq_len(nrow(ct)), rm)])
            })
        }
        
        obj@components <- pruneList(obj@components, checkZeroRows = TRUE) # remove empty components
        # update infos
        obj@componentInfo <- obj@componentInfo[name %in% names(obj@components)]
        
        if (!is.null(j))
            obj@componentInfo[, size := sapply(obj@components, nrow)][] # update
    }
    
    return(obj)
})

#' @describeIn components Returns all component data in a table.
#' @export
setMethod("as.data.table", "components", function(x)
{
    if (length(x) == 0)
        return(data.table())
    ret <- rbindlist(componentTable(x), idcol = "name")
    ret <- merge(componentInfo(x), ret, by = "name")
    ret <- ret[orderComponentsNames(name)]
    return(ret)
})

#' @describeIn components Provides rule based filtering for components.
#'
#' @param size Should be a two sized vector with the minimum/maximum size of a component. Set to \code{NULL} to ignore.
#' @param adducts Remove any feature groups within components that do not match given adduct rules. If \code{adducts} is
#'   a logical then only results are kept when an adduct is assigned (\code{adducts=TRUE}) or not assigned
#'   (\code{adducts=FALSE}). Otherwise, if \code{adducts} contains one or more \code{\link{adduct}} objects (or
#'   something that can be converted to it with \code{\link{as.adduct}}) then only results are kept that match the given
#'   adducts. Set to \code{NULL} to ignore this filter.
#' @param isotopes Only keep results that match a given isotope rule. If \code{isotopes} is a logical then only results
#'   are kept with (\code{isotopes=TRUE}) or without (\code{isotopes=FALSE}) isotope assignment. Otherwise
#'   \code{isotopes} should be a numeric vector with isotope identifiers to keep (\emph{e.g.} \samp{0} for monoisotopic
#'   results, \samp{1} for \samp{M+1} results etc.). Set to \code{NULL} to ignore this filter.
#' @param rtIncrement,mzIncrement Should be a two sized vector with the minimum/maximum retention or mz increment of a
#'   homologous series. Set to \code{NULL} to ignore.
#' @param checkComponentsSession If set then components and/or feature groups are removed that were selected for removal
#'   (see \link{check-GUI} and the \code{\link{checkComponents}} function). The value of \code{checkComponentsSession}
#'   should either by a path to the session file or \code{TRUE}, in which case the default session file name is used. If
#'   \code{negate=TRUE} then all non-selected data is removed instead.
#' @param negate If \code{TRUE} then filters are applied in opposite manner.
#' @param verbose If set to \code{FALSE} then no text output is shown.
#'
#' @note \code{filter} Applies only those filters for which a component has data available. For instance, filtering by
#'   adduct will only filter any results within a component if that component contains adduct information.
#'
#' @export
setMethod("filter", "components", function(obj, size = NULL, adducts = NULL, isotopes = NULL, rtIncrement = NULL,
                                           mzIncrement = NULL, checkComponentsSession = NULL, negate = FALSE,
                                           verbose = TRUE)
{
    if (isTRUE(checkComponentsSession))
        checkComponentsSession <- "checked-components.yml"
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertIntegerish(size, lower = 0, any.missing = FALSE, len = 2, null.ok = TRUE, add = ac)
    checkmate::assert(checkmate::checkFlag(isotopes, null.ok = TRUE),
                      checkmate::checkIntegerish(isotopes, lower = 0),
                      .var.name = "isotopes")
    checkmate::assertNumeric(rtIncrement, any.missing = FALSE, len = 2, null.ok = TRUE, add = ac)
    checkmate::assertNumeric(mzIncrement, lower = 0, any.missing = FALSE, len = 2, null.ok = TRUE, add = ac)
    if (!is.logical(checkComponentsSession))
        assertCheckSession(checkComponentsSession, mustExist = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertFlag(negate, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)

    if (!is.null(adducts) && !is.logical(adducts))
        adducts <- sapply(adducts, function(a) as.character(checkAndToAdduct(a)))

    old <- obj
    
    if (verbose)
        cat("Filtering components... ")

    if (!is.null(checkComponentsSession))
    {
        session <- readCheckSession(checkComponentsSession, "components")
        
        if (negate)
            obj <- obj[union(session$removeFully, names(session$removePartially))]
        else
            obj <- delete(obj, i = session$removeFully)
        
        if (length(session$removePartially) > 0)
        {
            obj <- delete(obj, i = names(session$removePartially), j = function(cmp, cName)
            {
                gRm <- cmp$group %chin% session$removePartially[[cName]]
                return(if (negate) !gRm else gRm)
            })
        }
    }
    
    inRange <- function(x, range) x >= range[1] & x <= range[2]
    mark <- if (negate) function(x) !x else function(x) x

    cInfo <- componentInfo(obj)
    
    obj <- delete(obj, j = function(cmp, cname, ...)
    {
        cmp <- copy(cmp)
        cmp[, keep := TRUE]
        if (!is.null(adducts) && !is.null(cmp[["adduct_ion"]]))
        {
            # NOTE: RAMClustR and CAMERA follow the generic adduct format
            # applied here, hence we can simply call as.adduct w/out specific
            # formatting.
            
            if (is.logical(adducts) && adducts)
                cmp[, keep := mark(!is.na(adduct_ion))]
            else if (is.logical(adducts) && !adducts)
                cmp[, keep := mark(is.na(adduct_ion))]
            else
                cmp[, keep := mark(adduct_ion %in% adducts)]
            
        }
        
        if (!is.null(isotopes) && !is.null(cmp[["isonr"]]))
        {
            if (is.logical(isotopes) && isotopes)
                cmp[keep == TRUE, keep := mark(!is.na(isonr))]
            else if (is.logical(isotopes) && !isotopes)
                cmp[keep == TRUE, keep := mark(is.na(isonr))]
            else
                cmp[keep == TRUE, keep := mark(isonr %in% isotopes)]
        }
        
        ci <- cInfo[name == cname]
        if (!is.null(size))
            cmp[keep == TRUE, keep := mark(inRange(ci$size, ..size))]
        if (!is.null(rtIncrement) && !is.null(ci[["ret_increment"]]))
            cmp[keep == TRUE, keep := mark(inRange(ci$ret_increment, rtIncrement))]
        if (!is.null(mzIncrement) && !is.null(ci[["mz_increment"]]))
            cmp[keep == TRUE, keep := mark(inRange(ci$mz_increment, mzIncrement))]
        
        return(!cmp$keep)
    })
    
    if (verbose)
        printComponentsFiltered(old, obj)

    return(obj)
})

#' @describeIn components Returns the component id(s) to which a feature group
#'   belongs.
#' @param fGroup The name (thus a character) of the feature group that should be
#'   searched for.
#' @aliases findFGroup
#' @export
setMethod("findFGroup", "components", function(obj, fGroup)
{
    checkmate::assertString(fGroup, min.chars = 1)
    if (length(obj) == 0)
        return(numeric())
    which(sapply(componentTable(obj), function(ct) fGroup %in% ct$group))
})

#' @describeIn components Plot a \emph{pseudo} mass spectrum for a single
#'   component.
#'
#' @param markFGroup If specified (\emph{i.e.} not \code{NULL}) this argument
#'   can be used to mark a feature group in the plotted spectrum. The value
#'   should be a character with the name of the feature group. Setting this to
#'   \code{NULL} will not mark any peak.
#'
#' @template plot-lim
#'
#' @export
setMethod("plotSpectrum", "components", function(obj, index, markFGroup = NULL, xlim = NULL, ylim = NULL, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assert(
        checkmate::checkInt(index, lower = 1, upper = length(componentTable(obj))),
        checkChoiceSilent(index, names(obj))
    , .var.name = "index")
    checkmate::assertString(markFGroup, min.chars = 1, null.ok = TRUE, add = ac)
    assertXYLim(xlim, ylim, add = ac)
    checkmate::reportAssertions(ac)

    plotData <- copy(componentTable(obj)[[index]]) # UNDONE: allow multiple selection?

    haveIso <- !is.null(plotData[["isogroup"]])
    haveAdd <- !is.null(plotData[["adduct_ion"]])
    haveHom <- !is.null(plotData[["hsnr"]])

    if (haveAdd && any(!is.na(plotData[["adduct_ion"]])))
    {
        # merge any adduct rows
        plotData[!is.na(adduct_ion), adduct_ion := paste0(adduct_ion, collapse = "/"), by = "group"]
        plotData <- plotData[!duplicated(plotData, by = c("group", "adduct_ion"))]
    }

    if (haveHom)
    {
        # merge merged homologue entries
        plotData[!is.na(hsnr), c("rGroup", "intensity") :=
                     .(paste0(rGroup, collapse = "/"),
                       max(intensity)), by = "hsnr"]
        plotData <- plotData[!duplicated(plotData, by = c("hsnr", "rGroup"))]
    }

    plotData[, label := sapply(seq_len(nrow(plotData)), function(r)
    {
        if (haveIso && !is.na(isogroup[r]) && isonr[r] != 0)
            sprintf("iso %d-[M+%d]", as.integer(isogroup[r]), as.integer(isonr[r]))
        # adduct?
        else if (haveAdd && !is.na(adduct_ion[r]))
            adduct_ion[r]
        # homologue?
        else if (haveHom && !is.na(hsnr[r]))
            paste0("HS #", hsnr[r])
        else # unknown
            ""
    })]

    plotData[, lwd := ifelse(nzchar(label), 1, 0.5)]

    if (haveHom)
        plotData[, categ := rGroup]
    else
        plotData[, categ := ifelse(nzchar(label), "assigned", "unassigned")]

    if (!is.null(markFGroup) && markFGroup %in% plotData$group)
    {
        plotData[group == markFGroup, c("categ", "lwd") := .(markFGroup, 2)]

        # convert to factor for sorting
        plotData[, categ := factor(categ, levels = c(markFGroup, unique(categ[categ != markFGroup])))]
    }
    else if (!haveHom)
        plotData[, categ := factor(categ, levels = c("assigned", "unassigned"))]

    plotData[nzchar(label), label := paste(label, round(mz, 4), sep = "\n")]
    plotData[!nzchar(label), label := as.character(round(mz, 4))]

    allCateg <- unique(plotData$categ)
    specCols <- getBrewerPal(length(allCateg), "Dark2")
    names(specCols) <- allCateg
    
    intMax <- max(plotData$intensity)
    mzRange <- range(plotData$mz)
    
    makeLegend <- function(x, y, ...) legend(x, y, allCateg, col = specCols[as.character(allCateg)],
                                             text.col = specCols[as.character(allCateg)], lty = 1,
                                             xpd = NA, ncol = 1, cex = 0.75, bty = "n", ...)
    
    oldp <- par(no.readonly = TRUE)
    plot.new()
    
    leg <- makeLegend(0, 0, plot = FALSE)
    lw <- (grconvertX(leg$rect$w, to = "ndc") - grconvertX(0, to = "ndc"))
    par(omd = c(0, 1 - lw, 0, 1), new = TRUE)
    
    if (is.null(xlim))
        xlim <- mzRange
    if (is.null(ylim))
        ylim <- c(0, intMax * 1.25)
    
    plot(0, xlab = "m/z", ylab = "Intensity", xlim = xlim, ylim = ylim,
         type = "n", bty = "l", ...)
    
    segments(plotData$mz, 0, plotData$mz, plotData$intensity,
             col = specCols[as.character(plotData$categ)], lwd = plotData$lwd * 2)
    
    tyOffset <- max(plotData$intensity * 0.02)
    tx <- plotData$mz
    ty <- plotData$intensity + tyOffset
    
    text(tx, ty, plotData$label, srt = 90, adj = 0, cex = 0.75)
    
    makeLegend(par("usr")[2], par("usr")[4])
    
    par(oldp)
})

setMethod("plotSpectrumHash", "components", function(obj, index, markFGroup = NULL, xlim = NULL, ylim = NULL, ...)
{
    return(makeHash(obj[[index]], markFGroup, xlim, ylim, ...))
})

#' @describeIn components Plot an extracted ion chromatogram (EIC) for all feature groups within a single component.
#' @param fGroups The \code{\link{featureGroups}} object that was used to generate the components.
#' @template EICParams-arg
#' @note For \code{plotChroms}: The \code{topMost} and \code{topMostByRGroup} EIC parameters are ignored unless the
#'   components are from homologous series.
#' @export
setMethod("plotChroms", "components", function(obj, index, fGroups, EICParams = getDefEICParams(rtWindow = 5), ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assert(
        checkmate::checkInt(index, lower = 1, upper = length(componentTable(obj))),
        checkChoiceSilent(index, names(obj))
    , .var.name = "index")
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    assertEICParams(EICParams)
    checkmate::reportAssertions(ac)

    comp <- componentTable(obj)[[index]]

    isHom <- !is.null(comp[["hsnr"]]) # homologues?

    showPeakArea <- isHom
    showFGroupRect <- !isHom
    colourBy = if (!isHom) "fGroup" else "rGroup"

    if (isHom)
    {
        rGroups <- unique(comp$rGroup)
        fGroups <- replicateGroupFilter(fGroups, rGroups, verbose = FALSE)
    }
    else
    {
        EICParams$topMost <- 1
        EICParams$topMostByRGroup <- FALSE
    }

    fGroups <- fGroups[, unique(comp$group)]

    if (length(fGroups) > 0)
        plotChroms(fGroups, EICParams = EICParams, colourBy = colourBy, showPeakArea = showPeakArea,
                   showFGroupRect = showFGroupRect, ...)
})

setMethod("plotChromsHash", "components", function(obj, index, fGroups, EICParams = getDefEICParams(rtWindow = 5), ...)
{
    comp <- componentTable(obj)[[index]]
    anas <- analyses(fGroups)
    if (!is.null(comp[["hsnr"]])) # homologues?
        anas <- analysisInfo(fGroups)[group %chin% comp$rGroup]$analysis
    makeHash(comp, plotChromsHash(fGroups, EICParams = EICParams, analyses = anas, groupName = comp$group, ...))
})

#' @describeIn components Generates a consensus from multiple \code{components}
#'   objects. At this point results are simply combined and no attempt is made to
#'   merge similar components.
#' @return \code{consensus} returns a \code{components} object that is produced
#'   by merging multiple specified \code{components} objects.
#' @export
setMethod("consensus", "components", function(obj, ...)
{
    allComponents <- c(list(obj), list(...))

    checkmate::assertList(allComponents, types = "components", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...")

    compNames <- make.unique(sapply(allComponents, algorithm))
    mcmp <- mergeComponents(allComponents, compNames, "algorithm")
    
    return(components(components = mcmp$components, componentInfo = mcmp$componentInfo,
                      algorithm = paste0(compNames, collapse = ",")))
})

#' Grouping feature groups in components
#'
#' Functionality to automatically group related feature groups (\emph{e.g.} isotopes, adducts and homologues) to assist
#' and simplify annotation.
#'
#' Several algorithms are provided to group feature groups that are related in some (chemical) way to each other. How
#' feature groups are related depends on the algorithm: examples include adducts, statistics and parents/transformation
#' products. The linking of this data is generally useful for annotation purposes and reducing data complexity.
#'
#' @templateVar func generateComponents
#' @templateVar what generate components
#' @templateVar ex1 generateComponentsRAMClustR
#' @templateVar ex2 generateComponentsNontarget
#' @templateVar algos ramclustr,camera,nontarget,intclust,openms,cliquems,specclust,tp
#' @templateVar algosSuffix RAMClustR,CAMERA,Nontarget,IntClust,OpenMS,CliqueMS,SpecClust,TPs
#' @templateVar ret components
#' @template generic-algo
#'
#' @param fGroups \code{\link{featureGroups}} object for which components should be generated.
#' @param \dots Any parameters to be passed to the selected component generation algorithm.
#'
#' @return A \code{\link{components}} (derived) object containing all generated components.
#'
#' @section Sets workflows: In a \link[=sets-workflow]{sets workflow} the componentization data is generated differently
#'   depending on the used algorithm. Please see the details in the algorithm specific functions linked in the \verb{See Also} section.
#'
#' @templateVar what generateComponents
#' @template main-rd-method
#' @export
setMethod("generateComponents", "featureGroups", function(fGroups, algorithm, ...)
{
    checkmate::assertChoice(algorithm, c("ramclustr", "camera", "nontarget", "intclust",
                                         "openms", "cliquems", "specclust", "tp"))
    
    f <- switch(algorithm,
                ramclustr = generateComponentsRAMClustR,
                camera = generateComponentsCAMERA,
                nontarget = generateComponentsNontarget,
                intclust = generateComponentsIntClust,
                openms = generateComponentsOpenMS,
                cliquems = generateComponentsCliqueMS,
                specclust = generateComponentsSpecClust,
                tp = generateComponentsTPs)

    f(fGroups, ...)
})
