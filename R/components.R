#' @include main.R
#' @include workflow-step.R
#' @include utils-components.R
NULL

#' Component class
#'
#' Contains data for feature groups that are related in some way. These
#' \emph{components} commonly include adducts, isotopes and homologues.
#'
#' \code{components} objects are obtained from
#' \link[=component-generation]{component generators}.
#'
#' @slot components List of all components in this object. Use the
#'   \code{componentTable} method for access.
#' @slot componentInfo A \code{\link{data.table}} containing general information
#'   for each component. Use the \code{componentInfo} method for access.
#'
#' @param obj,object,x The \code{component} object.
#' @param index The index of the component. Can be a numeric index or a
#'   character with its name.
#' @param \dots For \code{plotChroms}: Further (optional) arguments passed to the
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
#' @templateVar seli components
#' @templateVar selOrderi names()
#' @templateVar selj feature groups
#' @templateVar selOrderj groupNames()
#' @templateVar optionalj TRUE
#' @templateVar dollarOpName component
#' @template sub_op-args
#'
#' @return The subset operator (\code{"["}) and \code{filter} method return the
#'   data subset in an object from the \code{componentsReduced} class. This
#'   object does not contain any algorithm specific data and as such, algorithm
#'   specific methods (\emph{e.g.} \code{treeCut}) will not work on this object.
#'   The reason for this is that it is often very difficult or impossible to
#'   subset the algorithmic data.
#'
#' @templateVar class components
#' @template class-hierarchy
#'
#' @seealso \link{component-generation},  \code{\link{componentsNT}} and
#'   \code{\link{componentsIntClust}}
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

#' @rdname components-class
componentsReduced <- setClass("componentsReduced",
                              contains = "components")
setMethod("initialize", "componentsReduced",
          function(.Object, ...) callNextMethod(.Object, algorithm = "reduced", ...))


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
        gCounts <- 0
    else
        gCounts <- sapply(object@components, nrow)

    printf("Number of feature groups in components: %d (total), %.1f (mean), %d - %d (min - max)\n",
           sum(gCounts), mean(gCounts), min(gCounts), max(gCounts))
})

#' @describeIn components Subset on components/feature groups.
#' @export
setMethod("[", c("components", "ANY", "ANY", "missing"), function(x, i, j, ...)
{
    if (!missing(i))
        i <- assertSubsetArgAndToChr(i, names(x))
    if (!missing(j))
        j <- assertSubsetArgAndToChr(j, groupNames(x))

    # non-existing indices result in NULL values --> prune

    # for some reason we need to explicitly give a data.table() as otherwise it complains componentInfo was assigned a list()
    if (length(x) == 0)
        return(componentsReduced(components = list(), componentInfo = data.table()))

    if (!missing(i))
    {
        x@components <- pruneList(x@components[i])
        x@componentInfo <- x@componentInfo[name %in% i]
    }

    if (!missing(j))
    {
        x@components <- sapply(x@components, function(cmp) cmp[group %in% j],
                               simplify = FALSE)
        x@components <- x@components[sapply(x@components, nrow) > 0]
        x@componentInfo <- x@componentInfo[name %in% names(x@components)]
        x@componentInfo[, size := sapply(x@components, nrow)] # update in case groups were filtered away
    }

    return(componentsReduced(components = x@components, componentInfo = x@componentInfo))
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
#' @param size Should be a two sized vector with the minimum/maximum size of a
#'   component. Set to \code{NULL} to ignore.
#' @param adducts Remove any feature groups within components that do not match
#'   given adduct rules. If \code{adducts} is a logical then only results are
#'   kept when an adduct is assigned (\code{adducts=TRUE}) or not assigned
#'   (\code{adducts=FALSE}). Otherwise, if \code{adducts} contains one or more
#'   \code{\link{adduct}} objects (or something that can be converted to it with
#'   \code{\link{as.adduct}}) then only results are kept that match the given
#'   adducts. Set to \code{NULL} to ignore this filter.
#' @param isotopes Only keep results that match a given isotope rule. If
#'   \code{isotopes} is a logical then only results are kept with
#'   (\code{isotopes=TRUE}) or without (\code{isotopes=FALSE}) isotope
#'   assignment. Otherwise \code{isotopes} should be a numeric vector with
#'   isotope identifiers to keep (\emph{e.g.} \samp{0} for monoisotopic results,
#'   \samp{1} for \samp{M+1} results etc.). Set to \code{NULL} to ignore this
#'   filter.
#' @param rtIncrement,mzIncrement Should be a two sized vector with the
#'   minimum/maximum retention or mz increment of a homologous series. Set to
#'   \code{NULL} to ignore.
#' @param negate If \code{TRUE} then filters are applied in opposite manner.
#' @param verbose If set to \code{FALSE} then no text output is shown.
#'
#' @note \code{filter} Applies only those filters for which a component has data
#'   available. For instance, filtering by adduct will only filter any results
#'   within a component if that component contains adduct information.
#'
#' @export
setMethod("filter", "components", function(obj, size = NULL, adducts = NULL, isotopes = NULL, rtIncrement = NULL,
                                           mzIncrement = NULL, checkComponentsSession = NULL, negate = FALSE,
                                           verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertIntegerish(size, lower = 0, any.missing = FALSE, len = 2, null.ok = TRUE, add = ac)
    checkmate::assert(checkmate::checkFlag(isotopes, null.ok = TRUE),
                      checkmate::checkIntegerish(isotopes, lower = 0),
                      .var.name = isotopes)
    checkmate::assertNumeric(rtIncrement, lower = 0, any.missing = FALSE, len = 2, null.ok = TRUE, add = ac)
    checkmate::assertNumeric(mzIncrement, lower = 0, any.missing = FALSE, len = 2, null.ok = TRUE, add = ac)
    assertCheckComponentsSession(checkComponentsSession, obj, mustExist = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertFlag(negate, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)

    if (!is.null(adducts) && !is.logical(adducts))
        adducts <- sapply(adducts, function(a) as.character(checkAndToAdduct(a)))

    oldn <- length(obj); oldresn <- if (oldn > 0) sum(sapply(obj@components, nrow)) else 0
    
    if (verbose)
        cat("Filtering components... ")

    if (!is.null(checkComponentsSession))
    {
        session <- readRDS(getCheckSessionPath(checkComponentsSession, "components"))
        if (negate)
            obj <- obj[, setdiff(names(obj), session$primarySelections)]
        else
            obj <- obj[session$primarySelections]

        cNames <- names(obj)
        obj@components <- pruneList(Map(obj@components, session$secondarySelections[cNames], f = function(cmp, ec)
        {
            wh <- ec[match(cmp$group, session$secondarySelections$name)]
            if (negate)
                return(cmp[!wh])
            return(cmp[wh])
        }))
    }
    
    obj@components <- pruneList(lapply(obj@components, function(cmp)
    {
        if (!is.null(adducts) && !is.null(cmp[["adduct_ion"]]))
        {
            # NOTE: RAMClustR and CAMERA follow the generic adduct format
            # applied here, hence we can simply call as.adduct w/out specific
            # formatting.

            if (is.logical(adducts) && adducts)
                keep <- !is.na(cmp$adduct_ion)
            else if (is.logical(adducts) && !adducts)
                keep <- is.na(cmp$adduct_ion)
            else
                keep <- cmp$adduct_ion %in% adducts
            
            if (negate)
                keep <- !keep
            cmp <- cmp[keep]
        }

        if (!is.null(isotopes) && !is.null(cmp[["isonr"]]))
        {
            if (is.logical(isotopes) && isotopes)
                keep <- !is.na(cmp$isonr)
            else if (is.logical(isotopes) && !isotopes)
                keep <- is.na(cmp$isonr)
            else
                keep <- cmp$isonr %in% isotopes
            
            if (negate)
                keep <- !keep
            cmp <- cmp[keep]
        }

        return(cmp)
    }), checkZeroRows = TRUE)

    if (length(obj) == 0)
        obj@componentInfo <- data.table()
    else
    {
        obj@componentInfo <- obj@componentInfo[name %in% names(obj@components)]
        obj@componentInfo[, size := sapply(obj@components, nrow)] # update in case groups were filtered away

        inRange <- function(x, range) x >= range[1] & x <= range[2]
        if (negate)
            inRange <- Negate(inRange)
        
        if (!is.null(size))
        {
            csize <- size # rename for DT
            obj@componentInfo <- obj@componentInfo[inRange(size, csize)]
        }
        if (!is.null(rtIncrement) && !is.null(obj@componentInfo[["ret_increment"]]))
            obj@componentInfo <- obj@componentInfo[inRange(ret_increment, rtIncrement)]
        if (!is.null(mzIncrement) && !is.null(obj@componentInfo[["mz_increment"]]))
            obj@componentInfo <- obj@componentInfo[inRange(mz_increment, mzIncrement)]

        obj@components <- obj@components[names(obj@components) %in% obj@componentInfo$name]
    }

    if (verbose)
    {
        newn <- length(obj); newresn <- if (newn > 0) sum(sapply(obj@components, nrow)) else 0
        printf("Done! Filtered %d (%.2f%%) components and %d (%.2f%%) feature groups. Remaining: %d components with %d feature groups\n",
               oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100,
               oldresn - newresn, if (oldresn == 0) 0 else (1-(newresn/oldresn))*100,
               newn, newresn)
    }

    return(componentsReduced(components = obj@components, componentInfo = obj@componentInfo))
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
#' @template useGGplot2
#'
#' @template plot-lim
#'
#' @export
setMethod("plotSpectrum", "components", function(obj, index, markFGroup = NULL, useGGPlot2 = FALSE,
                                                 xlim = NULL, ylim = NULL, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assert(
        checkmate::checkInt(index, lower = 1, upper = length(componentTable(obj))),
        checkChoiceSilent(index, names(obj))
    , .var.name = index)
    checkmate::assertString(markFGroup, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertFlag(useGGPlot2, add = ac)
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

    if (!useGGPlot2)
    {
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
    }
    else
    {
        ret <- ggplot(plotData, aes_string(x = "mz", y = 0, label = "label")) +
            geom_segment(aes_string(xend = "mz", yend = "intensity", colour = "categ",
                                    size = "lwd")) + scale_size(range = range(plotData$lwd), guide = FALSE) +
            ggrepel::geom_text_repel(aes_string(y = "intensity", angle = 0), min.segment.length = 0.1,
                                     nudge_y = grid::convertUnit(grid::unit(5, "mm"), "npc", valueOnly = TRUE), size = 3.2) +

            xlab("m/z") + ylab("Intensity") +
            cowplot::theme_cowplot(font_size = 12) + theme(legend.position = "bottom", legend.title = element_blank())

        return(ret)
    }
})

setMethod("plotSpectrumHash", "components", function(obj, index, markFGroup = NULL, useGGPlot2 = FALSE,
                                                     xlim = NULL, ylim = NULL, ...)
{
    return(makeHash(obj[[index]], markFGroup, useGGPlot2, xlim, ylim, ...))
})

#' @describeIn components Plot an extracted ion chromatogram (EIC) for all
#'   feature groups within a single component.
#' @param fGroups The \code{\link{featureGroups}} object that was used to
#'   generate the components.
#' @param rtWindow Retention window: see the \code{plotChroms} method for the
#'   \code{\link{featureGroups}} class.
#' @export
setMethod("plotChroms", "components", function(obj, index, fGroups, rtWindow = 5, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assert(
        checkmate::checkInt(index, lower = 1, upper = length(componentTable(obj))),
        checkChoiceSilent(index, names(obj))
    , .var.name = index)
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertNumber(rtWindow, lower = 0, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    comp <- componentTable(obj)[[index]]

    isHom <- !is.null(comp[["hsnr"]]) # homologues?

    topMost <- if (!isHom) 1 else NULL
    showPeakArea <- isHom
    showFGroupRect <- !isHom
    colourBy = if (!isHom) "fGroup" else "rGroup"

    if (isHom)
    {
        rGroups <- unique(comp$rGroup)
        fGroups <- replicateGroupFilter(fGroups, rGroups, verbose = FALSE)
    }

    fGroups <- fGroups[, unique(comp$group)]

    if (length(fGroups) > 0)
        plotChroms(fGroups, rtWindow = rtWindow, colourBy = colourBy, showPeakArea = showPeakArea,
                   showFGroupRect = showFGroupRect, topMost = topMost, ...)
})

setMethod("plotChromsHash", "components", function(obj, index, fGroups, rtWindow = 5, ...)
{
    comp <- componentTable(obj)[[index]]
    if (!is.null(comp[["hsnr"]])) # homologues?
    {
        rGroups <- unique(comp$rGroup)
        fGroups <- replicateGroupFilter(fGroups, rGroups, verbose = FALSE)
    }
    fGroups <- fGroups[, unique(comp$group)]
    makeHash(comp, fGroups, rtWindow, ...)
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

    allComponents <- allComponents[lengths(allComponents) > 0]
    if (length(allComponents) < 2)
        stop("Need at least two non-empty components objects")

    compNames <- make.unique(sapply(allComponents, algorithm))
    mcmp <- mergeComponents(allComponents, compNames, "algorithm")
    
    return(components(components = mcmp$components, componentInfo = mcmp$componentInfo,
                      algorithm = paste0(compNames, collapse = ",")))
})

#' @templateVar func generateComponents
#' @templateVar what generate components
#' @templateVar ex1 generateComponentsRAMClustR
#' @templateVar ex2 generateComponentsNontarget
#' @templateVar algos ramclustr,camera,nontarget,intclust,openms,cliquems
#' @template generic-algo
#'
#' @param ... Any parameters to be passed to the selected component generation
#'   algorithm.
#'
#' @rdname component-generation
#' @aliases generateComponents
#' @export
setMethod("generateComponents", "featureGroups", function(fGroups, algorithm, ...)
{
    checkmate::assertChoice(algorithm, c("ramclustr", "camera", "nontarget", "intclust",
                                         "openms", "cliquems"))
    
    f <- switch(algorithm,
                ramclustr = generateComponentsRAMClustR,
                camera = generateComponentsCAMERA,
                nontarget = generateComponentsNontarget,
                intclust = generateComponentsIntClust,
                openms = generateComponentsOpenMS,
                cliquems = generateComponentsCliqueMS)

    f(fGroups, ...)
})
