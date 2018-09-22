#' @include main.R
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
#' @slot algorithm The algorithm that was used to generate the components. Use
#'   the \code{algorithm} method for access.
#'
#' @param obj,object,x The \code{component} object.
#' @param index The index of the component. Can be a numeric index or a
#'   character with its name.
#' @param \dots For \code{plotEIC}: Further (optional) arguments passed to the
#'   \code{plotEIC} method for the \code{\link{featureGroups}} class. Note that
#'   the \code{colourBy}, \code{showPeakArea}, \code{showFGroupRect} and
#'   \code{topMost} arguments cannot be set as these are set by this method. For
#'   \code{consensus}: \code{components} objects that should be used to generate
#'   the consensus.
#'
#' @return The subset/extract operator (\code{"["}) and \code{filter} method
#'   return the data subset in an object from the \code{componentsReduced}
#'   class. This object does not contain any algorithm specific data and as
#'   such, algorithm specific methods (\emph{e.g.} \code{treeCut}) will not work
#'   on this object. The reason for this is that it is often very difficult or
#'   impossible to subset the algorithmic data.
#'
#' @seealso \link{component-generation} and \code{\link{componentsIntClust}}
#'
#' @export
components <- setClass("components",
                       slots = c(components = "list", componentInfo = "data.table",
                                 algorithm = "character"),
                       prototype = c(components = list(), componentInfo = data.table(),
                                     algorithm = "none"))

componentsReduced <- setClass("componentsReduced",
                              contains = "components")

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

#' @describeIn components Retrieve the algorithm (a character string) used to
#'   generate components.
#' @export
setMethod("algorithm", "components", function(obj) obj@algorithm)

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
    printf("A components object (%s)\n", class(object))
    printf("Algorithm: %s\n", algorithm(object))
    printf("Components: %s (%d total)\n", getStrListWithMax(names(object), 6, ", "), length(object))

    if (length(object@components) == 0)
        gCounts <- 0
    else
        gCounts <- sapply(object@components, nrow)

    printf("Number of feature groups in components: %d (total), %.1f (mean), %d - %d (min - max)\n",
           sum(gCounts), mean(gCounts), min(gCounts), max(gCounts))

    showObjectSize(object)
})

#' @templateVar class components
#' @templateVar whati components
#' @templateVar orderi names()
#' @templateVar whatj feature groups
#' @templateVar orderj groupNames()
#' @template extr_op-args
#'
#' @export
setMethod("[", c("components", "ANY", "ANY", "ANY"), function(x, i, j, ..., drop = TRUE)
{
    if (!missing(i))
        assertExtractArg(i)
    if (!missing(j))
        assertExtractArg(j)
    
    # non-existing indices result in NULL values --> prune
    
    if (!missing(i))
    {
        if (!is.character(i))
            i <- names(x)[i]
        x@components <- pruneList(x@components[i])
        x@componentInfo <- x@componentInfo[name %in% i]
    }
    
    if (!missing(j))
    {
        if (!is.character(j))
            j <- groupNames(x)[j]
        x@components <- sapply(x@components, function(cmp) cmp[group == j],
                               simplify = FALSE)
        x@components <- x@components[sapply(x@components, nrow) > 0]
        x@componentInfo <- x@componentInfo[name %in% names(x@components)]
        x@componentInfo[, size := sapply(x@components, nrow)] # update in case groups were filtered away
    }
    
    return(componentsReduced(components = x@components, componentInfo = x@componentInfo, algorithm = "reduced"))
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

#' @templateVar class components
#' @templateVar withoutFGroups TRUE
#' @template filterby
#' @export
setMethod("filterBy", "components", function(obj, fGroups, negate, index = NULL)
{
    # UNDONE: filter by presence of isotopes, adducts, ...
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertFlag(negate, add = ac)
    checkmate::assert(
        checkmate::checkInt(index, lower = 1, upper = length(componentTable(obj))),
        checkChoiceSilent(index, names(obj)),
        checkmate::checkNull(index)
        , .var.name = index)
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
        grps <- character()
    else
    {
        cTable <- componentTable(obj)
        
        if (!is.null(index))
            cTable <- cTable[index]
        grps <- unique(unlist(sapply(cTable, "[[", "group")))
    }
    
    return(groupNamesFilter(fGroups, "components", grps, negate, hashParam = c(grps, negate, index)))
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
#' @export
setMethod("plotSpec", "components", function(obj, index, markFGroup = NULL, useGGPlot2 = FALSE, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assert(
        checkmate::checkInt(index, lower = 1, upper = length(componentTable(obj))),
        checkChoiceSilent(index, names(obj))
    , .var.name = index)
    checkmate::assertString(markFGroup, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertFlag(useGGPlot2, add = ac)
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

        plot(0, xlab = "m/z", ylab = "Intensity", xlim = mzRange, ylim = c(0, intMax * 1.25),
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

#' @describeIn components Plot an extracted ion chromatogram (EIC) for all
#'   feature groups within a single component.
#' @param fGroups The \code{\link{featureGroups}} object that was used to
#'   generate the components.
#' @param rtWindow Retention window: see the \code{plotEIC} method for the
#'   \code{\link{featureGroups}} class.
#' @export
setMethod("plotEIC", "components", function(obj, index, fGroups, rtWindow = 5, ...)
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
        plotEIC(fGroups, rtWindow = rtWindow, colourBy = colourBy, showPeakArea = showPeakArea,
                showFGroupRect = showFGroupRect, topMost = topMost, ...)
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
    retCInfo <- copy(componentInfo(allComponents[[1]]))
    retCInfo[, algorithm := compNames[[1]]]
    retCTable <- componentTable(allComponents[[1]])
    
    for (mi in seq(2, length(allComponents)))
    {
        if (length(allComponents[[mi]]) == 0)
            next
        rci <- copy(componentInfo(allComponents[[mi]]))
        rci[, algorithm := compNames[[mi]]]
        retCInfo <- rbind(retCInfo, rci, fill = TRUE)
        retCTable <- c(retCTable, componentTable(allComponents[[mi]]))
    }
    
    retCInfo[, name := paste0(name, "-", algorithm)]
    names(retCTable) <- retCInfo[["name"]]
    
    return(components(components = retCTable, componentInfo = retCInfo, algorithm = paste0(compNames, collapse = ",")))
})

#' @templateVar func generateComponents
#' @templateVar what generate components
#' @templateVar ex1 generateComponentsRAMClustR
#' @templateVar ex2 generateComponentsNontarget
#' @templateVar algos ramclustr,camera,nontarget,intclust
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
    f <- switch(algorithm,
                ramclustr = generateComponentsRAMClustR,
                camera = generateComponentsCAMERA,
                nontarget = generateComponentsNontarget,
                intclust = generateComponentsIntClust,
                stop("Invalid algorithm! Should be: ramclustr, camera, nontarget or intclust"))

    f(fGroups, ...)
})
