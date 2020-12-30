#' @include main.R
#' @include features.R
#' @include workflow-step.R
NULL

#' Base class for grouped features.
#'
#' This class holds all the information for grouped features.
#'
#' The \code{featureGroup} class is the workhorse of \pkg{patRoon}: almost all
#' functionality operate on its instantiated objects. The class holds all
#' information from grouped features (obtained from \code{\link{features}}).
#' This class itself is \code{virtual}, hence, objects are not created directly
#' from it. Instead, 'feature groupers' such as \code{\link{groupFeaturesXCMS}}
#' return a \code{featureGroups} derived object after performing the actual
#' grouping of features across analyses.
#'
#' @param fGroups,obj,x,object \code{featureGroups} object to be accessed.
#' @param retMin Plot retention time in minutes (instead of seconds).
#' @param \dots Ignored for \code{"["} operator or passed to
#'   \code{\link[graphics]{plot}} (\code{plot} and \code{plotChroms}),
#'   \code{\link[graphics]{lines}} (\code{plotInt}), \pkg{\link{VennDiagram}}
#'   plotting functions (\code{plotVenn}), \code{\link{chordDiagram}}
#'   (\code{plotChord}) or \code{\link[UpSetR]{upset}} (\code{plotUpSet}).
#' @param average Average data within replicate groups.
#' @param areas If set to \code{TRUE} then areas are considered instead of peak
#'   intensities.
#' @param pch,type,lty Common plotting parameters passed to \emph{e.g.}
#'   \code{\link[graphics]{plot}}. For \code{plot}: if \code{pch=NULL} then
#'   values are automatically assigned.
#' @param col Colour(s) used. If \code{col=NULL} then colours are automatically
#'   generated.
#' @param which A character vector with replicate groups used for comparison.
#'   
#'   For plotting functions: set to \code{NULL} for all replicate groups.
#'   
#'   For \code{plotVenn}: alternatively a named \code{list} containing elements
#'   of \code{character} vectors with replicate groups to compare. For instance,
#'   \code{which=list(infl = c("influent-A", "influent-B"), effl =
#'   c("effluent-A", "effluent-B"))}, will compare the features in replicate
#'   groups \samp{"influent-A/B"} against those in \samp{"effluent-A/B"}.
#'   The names of the list are used for labelling in the plot, and will be made
#'   automatically if not specified.
#' @param colourBy Sets the automatic colour selection: \code{"none"} for a
#'   single colour or \code{"rGroups"}/\code{"fGroups"} for a distinct colour
#'   per replicate/feature group.
#' @param showLegend If \code{TRUE} a legend will be shown with either replicate
#'   groups (\code{colourBy == "rGroups"}) or feature groups (\code{colourBy ==
#'   "fGroups"}, only for \code{plotChroms}). If \code{colourBy} is \code{"none"}
#'   no legend will be shown.
#'
#' @templateVar seli analyses
#' @templateVar selOrderi analyses()
#' @templateVar selj feature groups
#' @templateVar selOrderj names()
#' @templateVar optionalji TRUE
#' @templateVar dollarOpName feature group
#' @template sub_op-args
#'
#' @slot groups Matrix (\code{\link{data.table}}) with intensities for each
#'   feature group (columns) per analysis (rows). Access with \code{groups}
#'   method.
#' @slot analysisInfo,features \link[=analysis-information]{Analysis info} and
#'   \code{\link{features}} class associated with this object. Access with
#'   \code{analysisInfo} and \code{featureTable} methods, respectively.
#' @slot groupInfo \code{data.frame} with retention time (\code{rts} column, in
#'   seconds) and \emph{m/z} (\code{mzs} column) for each feature group. Access
#'   with \code{groupInfo} method.
#' @slot ftindex Matrix (\code{\link{data.table}}) with feature indices for each
#'   feature group (columns) per analysis (rows). Each index corresponds to the
#'   row within the feature table of the analysis (see
#'   \code{\link{featureTable}}).
#'
#' @templateVar class featureGroups
#' @template class-hierarchy
#'
#' @export
featureGroups <- setClass("featureGroups",
                          slots = c(groups = "data.table", analysisInfo = "data.frame", groupInfo = "data.frame",
                                    features = "features", ftindex = "data.table"),
                          contains = c("VIRTUAL", "workflowStep"))

setMethod("initialize", "featureGroups", function(.Object, ...)
{
    args <- list(...)

    # data.table's don't seem to initialize well (gives error that slot is init as list)
    if (is.null(args[["groups"]]))
        args$groups <- data.table()
    if (is.null(args[["ftindex"]]))
        args$ftindex <- data.table()

    do.call(callNextMethod, c(list(.Object), args))
})

#' @describeIn featureGroups Obtain feature group names.
#' @export
setMethod("names", "featureGroups", function(x) names(x@groups))

#' @templateVar class featureGroups
#' @templateVar what analyses
#' @template strmethod
#' @export
setMethod("analyses", "featureGroups", function(obj) analysisInfo(obj)$analysis)

#' @templateVar class featureGroups
#' @templateVar what replicate groups
#' @template strmethod
#' @export
setMethod("replicateGroups", "featureGroups", function(obj) unique(analysisInfo(obj)$group))

#' @describeIn featureGroups Same as \code{names}. Provided for consistency to other classes.
#' @export
setMethod("groupNames", "featureGroups", function(obj) names(obj))

#' @describeIn featureGroups Obtain number of feature groups.
#' @export
setMethod("length", "featureGroups", function(x) ncol(x@groups))

#' @describeIn featureGroups Shows summary information for this object.
#' @export
setMethod("show", "featureGroups", function(object)
{
    callNextMethod(object)
    anaInfo <- analysisInfo(object)
    printf("Feature groups: %s (%d total)\n", getStrListWithMax(names(object), 6, ", "),
           ncol(groupTable(object)))
    showAnaInfo(analysisInfo(object))
})

#' @describeIn featureGroups Accessor for \code{groups} slot.
#' @aliases groupTable
#' @export
setMethod("groupTable", "featureGroups", function(object, areas = FALSE)
{
    checkmate::assertFlag(areas)

    if (areas)
    {
        anaInfo <- analysisInfo(object)
        ret <- copy(object@groups)
        ftindex <- object@ftindex
        fTable <- featureTable(object)
        for (cl in seq_along(ret))
        {
            ftinds <- ftindex[[cl]]
            anainds <- seq_len(nrow(ret))[ftinds != 0]
            ftinds <- ftinds[ftinds != 0]
            as <- mapply(anainds, ftinds, SIMPLIFY = TRUE, FUN = function(a, i)
            {
                fTable[[anaInfo$analysis[a]]][["area"]][i]
            })
            set(ret, anainds, cl, as)
        }
        return(ret)
    }
    return(object@groups)
})

#' @describeIn featureGroups Obtain analysisInfo (see analysisInfo slot in \code{\link{features}}).
#' @export
setMethod("analysisInfo", "featureGroups", function(obj) obj@analysisInfo)

#' @describeIn featureGroups Accessor for \code{groupInfo} slot.
#' @aliases groupInfo
#' @export
setMethod("groupInfo", "featureGroups", function(fGroups) fGroups@groupInfo)

#' @describeIn featureGroups Obtain feature information (see \code{\link{features}}).
#' @export
setMethod("featureTable", "featureGroups", function(obj) featureTable(obj@features))

#' @describeIn featureGroups Accessor for \code{features} slot.
#' @export
setMethod("getFeatures", "featureGroups", function(obj) obj@features)

#' @describeIn featureGroups Accessor for \code{ftindex} slot.
#' @aliases groupFeatIndex
#' @export
setMethod("groupFeatIndex", "featureGroups", function(fGroups) fGroups@ftindex)

setMethod("removeAnalyses", "featureGroups", function(fGroups, indices)
{
    if (length(indices) > 0)
    {
        if (length(fGroups@groups) > 0)
        {
            fGroups@groups <- fGroups@groups[-indices]
            fGroups@ftindex <- fGroups@ftindex[-indices]
        }
        fGroups@analysisInfo <- fGroups@analysisInfo[-indices, ]
        fGroups@features <- fGroups@features[-indices]
    }
    return(fGroups)
})

setMethod("removeGroups", "featureGroups", function(fGroups, indices)
{
    if (length(indices) > 0)
    {
        if (length(fGroups@groups) > 0)
        {
            fGroups@groups <- fGroups@groups[, -indices, with = FALSE]
            fGroups@ftindex <- fGroups@ftindex[, -indices, with = FALSE]
        }
        fGroups@groupInfo <- fGroups@groupInfo[-indices, ]
    }
    return(fGroups)
})

#' @describeIn featureGroups Subset on analyses/feature groups.
#' @param rGroups An optional \code{character} vector: if specified only keep
#'   results for the given replicate groups (equivalent to the \code{rGroups}
#'   argument to \code{\link[=filter,featureGroups-method]{filter}}).
#' @export
setMethod("[", c("featureGroups", "ANY", "ANY", "missing"), function(x, i, j, ..., rGroups, drop = TRUE)
{
    if (!missing(rGroups))
        x <- filter(x, rGroups = rGroups)

    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, analyses(x))
        x <- removeAnalyses(x, which(!analyses(x) %in% i))
    }

    if (!missing(j))
    {
        j <- assertSubsetArgAndToChr(j, names(x))
        x <- removeGroups(x, which(!names(x) %in% j))
    }

    return(removeEmptyGroups(x))
})

#' @describeIn featureGroups Extract intensity values.
#' @export
setMethod("[[", c("featureGroups", "ANY", "ANY"), function(x, i, j)
{
    assertExtractArg(i)

    if (missing(j))
        return(x@groups[[i]])

    assertExtractArg(j)
    if (is.character(i))
        i <- match(i, analyses(x))
    return(x@groups[[i, j]])
})

#' @describeIn featureGroups Extract intensity values for a feature group.
#' @export
setMethod("$", "featureGroups", function(x, name)
{
    eval(substitute(x@groups$NAME_ARG, list(NAME_ARG = name)))
})

setMethod("removeEmptyGroups", "featureGroups", function(fGroups)
{
    if (length(fGroups) > 0)
    {
        empty <- unlist(fGroups@groups[, lapply(.SD, function(x) sum(x) == 0)])

        if (any(empty))
            fGroups <- removeGroups(fGroups, which(empty))
    }
    return(fGroups)
})

# UNDONE: make this public?
setMethod("removeEmptyAnalyses", "featureGroups", function(fGroups)
{
    if (length(fGroups) > 0)
    {
        trGT <- transpose(groupTable(fGroups))

        empty <- trGT[, sapply(.SD, sum) == 0]
        if (any(empty))
            fGroups <- removeAnalyses(fGroups, which(empty))
    }
    return(fGroups)
})

setMethod("averageGroups", "featureGroups", function(fGroups, areas)
{
    gTable <- copy(groupTable(fGroups, areas))
    if (nrow(gTable) == 0)
        return()

    gNames <- names(fGroups)
    anaInfo <- analysisInfo(fGroups)

    gTable[, sgroup := anaInfo$group]

    gTable[, (gNames) := lapply(.SD, function(v) { if (any(v > 0)) mean(v[v>0]) else 0 }), by = sgroup, .SDcols = gNames]
    gTable <- unique(gTable, by = "sgroup")
    gTable[, sgroup := NULL]

    return(gTable[])
})

setMethod("updateFeatIndex", "featureGroups", function(fGroups)
{
    # remove feature indices from feature groups that were removed later (e.g. filtered out)

    gTable <- groupTable(fGroups)
    fGroups@ftindex <- copy(groupFeatIndex(fGroups))
    for (i in seq_along(gTable))
        set(fGroups@ftindex, which(gTable[[i]] == 0), i, 0)
    return(fGroups)
})

#' @describeIn featureGroups Exports feature groups to a \file{.csv} file that
#'   is readable to Bruker ProfileAnalysis (a 'bucket table'), Bruker TASQ (an
#'   analyte database) or that is suitable as input for the \verb{Targeted peak
#'   detection} functionality of \href{http://mzmine.github.io/}{MZmine}.
#' @param out The destination file for the exported data.
#' @aliases export
#' @export
setMethod("export", "featureGroups", function(fGroups, type, out)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertChoice(type, c("brukerpa", "brukertasq", "mzmine"), add = ac)
    checkmate::assertPathForOutput(out, overwrite = TRUE, add = ac) # NOTE: assert doesn't work on Windows...
    checkmate::reportAssertions(ac)

    if (length(fGroups) == 0)
        stop("Cannot export empty feature groups object")

    if (type == "brukerpa")
    {
        # UNDONE: do we need this?
        #files <- sapply(bucketInfo$fInfo$analysis, function(f) file.path(bucketInfo$dataPath, paste0(f, ".d")), USE.NAMES = F)
        files <- fGroups@analysisInfo$analysis

        # col.names: if NA an empty initial column is added
        write.table(fGroups@groups, out, na = "", sep = "\t", quote = FALSE, row.names = files, col.names = NA)
    }
    else if (type == "brukertasq")
    {
        hdr <- c("name", "formula", "m/z", "rt", "rn", "CAS", "Qual1", "Qual2", "Qual3", "Qual4",
                 "Qual5", "Qual6", "quantifier ion 1 low", "quantifier ion 1 up",
                 "quantifier ion 1 formula", "quantifier ion 2 low", "quantifier ion 2 up",
                 "quantifier ion 2 formula", "quantifier ion 3 low", "quantifier ion 3 up",
                 "quantifier ion 3 formula", "precursor ion MS2", "precursor ion MS2 formula",
                 "precursor ion MS3", "precursor ion MS3 formula", "precursor ion MS4",
                 "precursor ion MS4 formula")

        df <- data.frame(matrix(ncol=length(hdr), nrow=ncol(fGroups@groups)))
        colnames(df) <- hdr

        df["name"] <- colnames(fGroups@groups)
        df["m/z"] <- fGroups@groupInfo$mzs
        df["rt"] <- fGroups@groupInfo$rts / 60

        write.csv(df, out, row.names = FALSE, na = "")
    }
    else if (type == "mzmine")
    {
        df <- fGroups@groupInfo
        df$name <- rownames(df)
        df <- df[, c("mzs", "rts", "name")]
        df$rts <- df$rts / 60
        write.table(df, out, row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
    }
})

#' @describeIn featureGroups Obtain a summary table (a \code{\link{data.table}})
#'   with retention, \emph{m/z}, intensity and optionally other feature data.
#' @param features If \code{TRUE} then feature specific data will be added. If
#'   \code{average=TRUE} this data will be averaged for each feature group.
#' @param regression Set to \code{TRUE} to add regression data for each feature
#'   group. For this a linear model is created (intensity/area [depending on
#'   \code{areas} argument] \emph{vs} concentration). The model concentrations
#'   (e.g. of a set of standards) is derived from the \code{conc} column of the
#'   \link[=analysis-information]{analysis information}. From this model the
#'   intercept, slope and R2 is added to the output. In addition, when
#'   \code{features=TRUE}, concentrations for each feature are added. Note that
#'   no regression information is added when no \code{conc} column is present in
#'   the analysis information or when less than two concentrations are specified
#'   (\emph{i.e.} the minimum amount).
#' @export
setMethod("as.data.table", "featureGroups", function(x, average = FALSE, areas = FALSE, features = FALSE, regression = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(average, add = ac)
    checkmate::assertFlag(areas, add = ac)
    checkmate::assertFlag(features, add = ac)
    checkmate::assertFlag(regression, add = ac)
    checkmate::reportAssertions(ac)

    if (length(x) == 0)
        return(data.table(mz = numeric(), ret = numeric(), group = character()))

    if (features && average && regression)
        stop("Cannot add regression data for averaged features.")

    anaInfo <- analysisInfo(x)
    gNames <- names(x)
    gInfo <- groupInfo(x)
    doConc <- regression && !is.null(anaInfo[["conc"]]) && sum(!is.na(anaInfo[["conc"]]) > 1)

    if (regression && is.null(anaInfo[["conc"]]))
        warning("No concentration information specified in the analysis information (i.e. conc column, see ?`analysis-information`)")

    if (features)
    {
        ftindex <- groupFeatIndex(x)
        fTable <- featureTable(x)
        snames <- anaInfo$analysis

        ret <- rbindlist(lapply(gNames, function(grp)
        {
            rbindlist(lapply(snames, function(s)
            {
                fTable[[s]][ftindex[[grp]][match(s, snames)]]
            }), idcol = "analysis")
        }), idcol = "group")

        if (doConc)
            ret[, conc := anaInfo$conc[analysis]]

        ret[, analysis := snames[analysis]]
        ret[, group := gNames[group]]

        if (average)
        {
            ret <- ret[, -c("isocount", "analysis", "ID")]
            numCols <- setdiff(names(ret), c("group"))
            ret[, (numCols) := lapply(.SD, mean), .SDcols = numCols, by = "group"]
            ret <- unique(ret, by = "group")
        }
        else
        {
            doConc <- doConc && length(snames) > 1
            if (doConc)
            {
                ret[, c("RSQ", "intercept", "slope") :=
                        {
                            notna <- !is.na(conc)
                            if (sum(notna) < 2)
                                NA_real_
                            else
                            {
                                reg <- summary(lm(intensity[notna] ~ conc[notna]))
                                list(reg[["r.squared"]], reg[["coefficients"]][1, 1], reg[["coefficients"]][2, 1])
                            }
                        }, by = group]
                ret[, conc_reg := (intensity - intercept) / slope] # y = ax+b
            }
        }

        ret[, c("group_ret", "group_mz") := gInfo[group, c("rts", "mzs")]]
        setcolorder(ret, c("group", "group_ret", "group_mz"))
    }
    else
    {
        if (average)
        {
            gTable <- averageGroups(x, areas)
            snames <- unique(anaInfo$group)
            if (doConc)
                concs <- anaInfo[!duplicated(anaInfo$group), "conc"] # conc should be same for all replicates
        }
        else
        {
            gTable <- groupTable(x, areas)
            snames <- anaInfo$analysis
            if (doConc)
                concs <- anaInfo$conc
        }

        ret <- transpose(gTable)
        setnames(ret, snames)

        doConc <- doConc && length(snames) > 1 && sum(!is.na(concs)) > 1
        if (doConc)
        {
            notna <- !is.na(concs)
            notnaconcs <- concs[notna]
            regr <- lapply(gTable, function(grp) summary(lm(grp[notna] ~ notnaconcs)))

            ret[!sapply(regr, is.null), c("RSQ", "intercept", "slope") :=
                    .(sapply(regr, "[[", "r.squared"),
                      sapply(regr, function(r) r$coefficients[1, 1]),
                      sapply(regr, function(r) r$coefficients[2, 1]))]
        }

        ret[, c("group", "ret", "mz") := .(gNames, gInfo$rts, gInfo$mzs)]
        setcolorder(ret, c("group", "ret", "mz"))
    }

    return(ret[])
})

#' @describeIn featureGroups Generates an \emph{m/z} \emph{vs} retention time
#'   plot for all featue groups. Optionally highlights unique/overlapping
#'   presence amongst replicate groups.
#' @param onlyUnique If \code{TRUE} and \code{colourBy="rGroups"} then only
#'   feature groups that are unique to a replicate group are plotted.
#' @export
setMethod("plot", "featureGroups", function(x, colourBy = c("none", "rGroups", "fGroups"),
                                            onlyUnique = FALSE, retMin = FALSE,
                                            showLegend = TRUE, col = NULL,
                                            pch = NULL, ...)
{
    rGroups <- replicateGroups(x)
    if (is.null(which))
        which <- rGroups

    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ onlyUnique + retMin + showLegend, fixed = list(add = ac))
    colourBy <- checkmate::matchArg(colourBy, c("none", "rGroups", "fGroups"), add = ac)
    checkmate::reportAssertions(ac)

    if (length(x) == 0)
        plot(0, type = "n", ...)
    else
    {
        if (colourBy == "fGroups" || colourBy == "none")
        {
            if (is.null(col))
                col <- if (colourBy == "fGroups") colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(x)) else "black"
            if (is.null(pch))
                pch <- 16

            if (colourBy == "fGroups" && showLegend)
            {
                labels <- names(x)
                labCol <- rep(col, length.out = length(labels))
                labPch <- rep(pch, length.out = length(labels))
                names(labCol) <- labels; names(labPch) <- labels
            }
            else
                showLegend <- FALSE
        }
        else if (colourBy == "rGroups")
        {
            labels <- c(replicateGroups(x), "overlap")

            if (is.null(col))
                labCol <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(labels))
            else
                labCol <- rep(col, length.out = length(labels))
            names(labCol) <- labels

            if (is.null(pch))
            {
                # prefer closed symbols (15-25)
                ll <- length(labels)
                if (ll < (25-15))
                    labPch <- seq_len(ll) + 14
                else if (ll <= 25)
                    labPch <- seq_len(ll)
                else
                    labPch <- rep(16, ll) # just stick with one
            }
            else
                labPch <- rep(pch, length.out = length(labels))
            names(labPch) <- labels

            # get averaged intensities for each rGroup and omit initial name/rt/mz columns
            gTable <- as.data.table(x, average = TRUE)[, replicateGroups(x), with = FALSE]

            for (r in seq_len(nrow(gTable)))
            {
                present <- which(gTable[r] != 0)
                set(gTable, r, "label", if (length(present) > 1) "overlap" else labels[present])
            }

            if (onlyUnique)
            {
                keep <- gTable$label != "overlap"
                gTable <- gTable[keep]; x <- x[, keep]
            }

            # remove unassigned (e.g. rGroups without unique feature groups)
            labels <- labels[labels %in% gTable$label]

            col <- labCol[gTable$label]; pch <- labPch[gTable$label]
        }

        oldp <- par(no.readonly = TRUE)
        if (showLegend)
        {
            makeLegend <- function(x, y, ...)
            {
                return(legend(x, y, labels, col = labCol[labels], pch = labPch[labels],
                              text.col = labCol[labels],
                              xpd = NA, ncol = 1, cex = 0.75, bty = "n", ...))
            }

            plot.new()
            leg <- makeLegend(0, 0, plot = FALSE)
            lw <- (grconvertX(leg$rect$w, to = "ndc") - grconvertX(0, to = "ndc"))
            par(omd = c(0, 1 - lw, 0, 1), new = TRUE)
        }

        plot(if (retMin) x@groupInfo$rts / 60 else x@groupInfo$rts, x@groupInfo$mzs,
             xlab = if (retMin) "Retention time (min)" else "Retention time (sec.)",
             ylab = "m/z", col = col, pch = pch, ...)

        if (showLegend)
            makeLegend(par("usr")[2], par("usr")[4])

        par(oldp)
    }
})

#' @describeIn featureGroups Generates a line plot for the (averaged) intensity
#'   of feature groups within all analyses
#' @export
setMethod("plotInt", "featureGroups", function(obj, average = FALSE, pch = 20, type = "b", lty = 3, col = NULL, ...)
{
    checkmate::assertFlag(average)

    if (length(obj) == 0)
    {
        plot(0, type = "n")
        invisible(return(NULL))
    }

    anaInfo <- analysisInfo(obj)

    if (average)
    {
        gTable <- averageGroups(obj)
        snames <- unique(anaInfo$group)
    }
    else
    {
        gTable <- groupTable(obj)
        snames <- anaInfo$analysis
    }

    nsamp <- length(snames)

    plot(x = c(0, nsamp), y = c(0, max(gTable)), type = "n", xlab = "", ylab = "Intensity", xaxt = "n")
    axis(1, seq_len(nsamp), snames, las = 2)

    if (is.null(col))
        col <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(gTable))

    px <- seq_len(nsamp)
    for (i in seq_along(gTable))
        lines(x = px, y = gTable[[i]], type = type, pch = pch, lty = lty, col = col[i], ...)
})

setMethod("plotIntHash", "featureGroups", function(obj, average = FALSE, ...) makeHash(allArgs()))

#' @describeIn featureGroups Generates a chord diagram which can be used to
#'   visualize shared presence of feature groups between analyses or replicate
#'   groups. In addition, analyses/replicates sharing similar properties
#'   (\emph{e.g.} location, age, type) may be grouped to enhance visualization
#'   between these 'outer groups'.
#'
#' @param outerGroups Character vector of names to be used as outer groups. The
#'   values in the specified vector should be named by analysis names
#'   (\code{average} set to \code{FALSE}) or replicate group names
#'   (\code{average} set to \code{TRUE}), for instance: \code{c(analysis1 =
#'   "group1", analysis2 = "group1", analysis3 = "group2")}. Set to \code{NULL}
#'   to disable outer groups.
#' @param addIntraOuterGroupLinks If \code{TRUE} then links will be added within
#'   outer groups.
#'
#' @template plotChord-args
#'
#' @references \addCitations{circlize}{1}
#'
#' @export
setMethod("plotChord", "featureGroups", function(obj, addSelfLinks = FALSE, addRetMzPlots = TRUE, average = FALSE,
                                                 outerGroups = NULL, addIntraOuterGroupLinks = FALSE, ...)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ addSelfLinks + addRetMzPlots + average + addIntraOuterGroupLinks,
           fixed = list(add = ac))
    checkmate::assertCharacter(outerGroups, min.chars = 1, min.len = 2, names = "unique", null.ok = TRUE)
    checkmate::reportAssertions(ac)

    if (length(obj) == 0)
        stop("Can't plot empty feature groups object")

    hasOuter <- !is.null(outerGroups)

    obj <- removeEmptyAnalyses(obj)

    anaInfo <- analysisInfo(obj)
    gInfo <- groupInfo(obj)

    if (average)
        snames <- unique(anaInfo$group)
    else
        snames <- anaInfo$analysis

    if (length(snames) < 2)
        stop(sprintf("Nothing to compare: need multiple (non-empty) %s!", if (average) "replicate groups" else "analyses"))

    if (hasOuter && !all(snames %in% names(outerGroups)))
        stop(sprintf("The following %s have no outerGroups assigned: %s", if (average) "replicate groups" else "analyses",
             paste0(setdiff(snames, names(outerGroups)), collapse = ", ")))

    nsamp <- length(snames)

    chordTable <- rbindlist(lapply(seq_along(snames),
                                   function(sni) data.table(from = snames[sni], to = snames[seq(sni, length(snames))])))

    getGTable <- function(snlist = c())
    {
        if (average)
        {
            if (length(snlist) > 0)
                fgf <- replicateGroupFilter(obj, snlist, verbose = FALSE)
            else
                fgf <- obj
            if (length(fgf) == 0)
                return(data.table())
            return(averageGroups(fgf))
        }
        else
        {
            if (length(snlist) > 0)
                return(groupTable(obj[snlist]))
            return(groupTable(obj))
        }
    }

    getLinkScore <- function(sn1, sn2)
    {
        if (sn1 == sn2)
            return(0)

        gTable <- getGTable(c(sn1, sn2))
        if (nrow(gTable) == 0)
            return(0)

        # Count all feature groups that are present in both samples/groups
        return(sum(gTable[, sapply(.SD, function(rows) all(rows > 0))]))
    }

    chordTable[, value := as.integer(Vectorize(getLinkScore)(from, to))]

    if (addSelfLinks)
    {
        gt <- getGTable()
        uniqueLinkCount <- sapply(seq_along(snames),
                                  function(sni) sum(sapply(gt, function(ints) ints[sni] > 0 && all(ints[-sni] == 0))))
        chordTable[from == to, value := uniqueLinkCount[.GRP], by = from]
    }

    if (hasOuter)
    {
        chordTable[, groupFrom := outerGroups[from]]
        chordTable[, groupTo := outerGroups[to]]
        if (!addIntraOuterGroupLinks)
            chordTable[groupFrom == groupTo & from != to, value := 0] # clear links within same groups (except self links)
        setorder(chordTable, groupFrom)

        remainingSN <- unique(unlist(chordTable[value != 0, .(from, to)])) # assigned samples, others will be removed
        og <- outerGroups[remainingSN] # outer groups assigned to each remaining sample
        gaps <- rep(1, length(og)) # initialize gaps
        gaps[cumsum(sapply(unique(og), function(x) length(og[og == x])))] <- 8 # make gap bigger after each outer group
        circlize::circos.par(gap.after = gaps)
    }

    if (all(chordTable$value == 0))
        stop("Did not found any overlap! Nothing to plot.")

    tracks <- NULL
    if (hasOuter)
        tracks <- list(list(track.height = 0.1, track.margin = c(if (addRetMzPlots) 0.05 else 0.06, 0)))
    if (addRetMzPlots)
        tracks <- c(tracks, list(list(track.height = 0.1, track.margin = c(0.08, 0))))

    maxv <- max(if (hasOuter) chordTable[groupFrom != groupTo, value] else chordTable$value)
    colFunc <- circlize::colorRamp2(maxv * seq(0, 1, 0.25),
                                    c("blue4", "deepskyblue1", "green", "orange", "red"),
                                    transparency = 0.5)

    if (hasOuter && addIntraOuterGroupLinks)
    {
        colFuncWithin <- circlize::colorRamp2(range(chordTable[groupFrom == groupTo, value]),
                                    c("grey80", "grey60"), transparency = 0.7)
        linkColors <- chordTable[, ifelse(groupFrom == groupTo, colFuncWithin(value), colFunc(value))]
    }
    else
        linkColors <- chordTable[, colFunc(value)]

    cdf <- circlize::chordDiagram(chordTable[, 1:3], annotationTrack = c("grid", "axis"),
                                  preAllocateTracks = tracks,
                                  grid.col = getBrewerPal(nsamp, "Dark2"),
                                  col = linkColors,
                                  annotationTrackHeight = c(0.1, 0.09),
                                  ...)

    circlize::circos.track(track.index = length(tracks) + 1, panel.fun = function(x, y)
    {
        sector.index <- circlize::get.cell.meta.data("sector.index")
        xlim <- circlize::get.cell.meta.data("xlim")
        ylim <- circlize::get.cell.meta.data("ylim")
        circlize::circos.text(mean(xlim), mean(ylim), sector.index, col = "white", cex = 1, niceFacing = TRUE)
    }, bg.border = NA)

    if (addRetMzPlots)
    {
        retMz <- rbindlist(sapply(unique(cdf$rn), function(sn)
        {
            ftgrps <- colnames(getGTable(sn))
            return(gInfo[ftgrps, ])
        }, simplify = FALSE), idcol = "sname")
        retMz$rts <- retMz$rts / max(retMz$rts) # normalize

        circlize::circos.track(fa = retMz$sname, x = retMz$rts, y = retMz$mzs, ylim = c(0, max(retMz$mzs)), track.index = length(tracks),
                               panel.fun = function(x, y)
                               {
                                   x <- x / (max(x) / circlize::get.cell.meta.data("xrange"))
                                   circlize::circos.points(x, y, cex = 0.5, col = "blue", pch = 16)
                               })
    }

    if (hasOuter)
    {
        finalChordTable <- chordTable[from %in% cdf$rn]
        finalOuterGroups <- unique(finalChordTable$groupFrom)

        ogcol <- getBrewerPal(length(finalOuterGroups), "Paired")
        for (ogi in seq_along(finalOuterGroups))
            circlize::highlight.sector(unique(finalChordTable[groupFrom == finalOuterGroups[ogi], from]), track.index = 1, col = ogcol[ogi],
                                       text = finalOuterGroups[ogi], cex = 1, text.col = "white", niceFacing = TRUE)
    }

    circlize::circos.clear()
})

#' @describeIn featureGroups Plots extracted ion chromatograms (EICs) of feature
#'   groups.
#'
#' @param rtWindow Retention time (in seconds) that will be subtracted/added to
#'   respectively the minimum and maximum retention time of the plotted feature
#'   groups. Thus, setting this value to a positive value will 'zoom out' on the
#'   retention time axis.
#' @param mzWindow The \emph{m/z} value (in Da) which will be subtracted/added
#'   to a feature group \emph{m/z} value to determine the width of its EIC.
#' @param topMost Only plot EICs from features within this number of top most
#'   intense analyses. If \code{NULL} then all analyses are used for plotted.
#' @param topMostByRGroup If set to \code{TRUE} and \code{topMost} is set: only
#'   plot EICs for the top most features in each replicate group. For instance,
#'   when \code{topMost=1} and \code{topMostByRGroup=TRUE}, then EICs will be
#'   plotted for the most intense feature of each replicate group.
#' @param EICs Internal parameter for now and should be kept at \code{NULL}
#'   (default).
#' @param showPeakArea Set to \code{TRUE} to display integrated chromatographic
#'   peak ranges by filling (shading) their areas.
#' @param showFGroupRect Set to \code{TRUE} to mark the full retention/intensity
#'   range of all features within a feature group by drawing a rectangle around
#'   it.
#' @param title Character string used for title of the plot. If \code{NULL} a
#'   title will be automatically generated.
#' @param onlyPresent If \code{TRUE} then EICs will only be generated for
#'   analyses in which a particular feature group was detected. Disabling this
#'   option might be useful to see if any features were 'missed'.
#' @param annotate If set to \code{"ret"} and/or \code{"mz"} then retention
#'   and/or \emph{m/z} values will be drawn for each plotted feature group.
#' @param showProgress if set to \code{TRUE} then a text progressbar will be
#'   displayed when all EICs are being plot. Set to \code{"none"} to disable any
#'   annotation.
#'
#' @template plot-lim
#'
#' @export
setMethod("plotChroms", "featureGroups", function(obj, rtWindow = 30, mzExpWindow = 0.001, retMin = FALSE, topMost = NULL,
                                                  topMostByRGroup = FALSE, EICs = NULL, showPeakArea = FALSE,
                                                  showFGroupRect = TRUE, title = NULL,
                                                  colourBy = c("none", "rGroups", "fGroups"),
                                                  showLegend = TRUE, onlyPresent = TRUE,
                                                  annotate = c("none", "ret", "mz"), showProgress = FALSE,
                                                  xlim = NULL, ylim = NULL, ...)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertNumber, . ~ rtWindow + mzExpWindow, lower = 0, finite = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertFlag, . ~ retMin + topMostByRGroup + showPeakArea + showFGroupRect + showLegend +
               onlyPresent + showProgress,
           fixed = list(add = ac))
    checkmate::assertCount(topMost, positive = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertString(title, null.ok = TRUE, add = ac)

    colourBy <- checkmate::matchArg(colourBy, c("none", "rGroups", "fGroups"), add = ac)
    annotate <- checkmate::matchArg(annotate, c("none", "ret", "mz"), several.ok = TRUE, add = ac)

    assertXYLim(xlim, ylim, add = ac)

    checkmate::reportAssertions(ac)

    if (length(obj) == 0)
    {
        plot(0, type = "n")
        invisible(return(NULL))
    }

    if (showLegend && colourBy == "none")
        showLegend <- FALSE

    fTable <- featureTable(obj)
    gTable <- groupTable(obj)
    gInfo <- groupInfo(obj)
    gCount <- nrow(gInfo)
    gNames <- names(obj)
    anaInfo <- analysisInfo(obj)
    ftind <- groupFeatIndex(obj)

    rGroups <- unique(anaInfo$group)

    if (is.null(EICs))
        EICs <- getEICsForFGroups(obj, rtWindow, mzExpWindow, topMost, onlyPresent)
    EICFGroups <- unique(unlist(sapply(EICs, names)))

    if (colourBy == "rGroups")
    {
        EICColors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(rGroups))
        names(EICColors) <- rGroups
    }
    else if (colourBy == "fGroups")
    {
        EICColors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(gCount)
        names(EICColors) <- gNames
    }
    else
        EICColors <- "blue"

    fillColors <- adjustcolor(EICColors, alpha.f = 0.35)
    names(fillColors) <- names(EICColors)

    # get overall retention/intensity limits
    plotLimits <- list(rtRange = c(0, 0), maxInt = 0)

    for (grpi in seq_len(gCount))
    {
        rtrs <- unlist(sapply(seq_len(nrow(anaInfo)), function(anai)
        {
            ana <- anaInfo$analysis[anai]
            fti <- ftind[[grpi]][anai]
            if (fti == 0)
                return(NA)
            return(unlist(fTable[[ana]][fti, .(retmin, retmax)]))
        }))

        if (any(!is.na(rtrs)))
        {
            rtmin <- min(rtrs, na.rm = TRUE)
            rtmax <- max(rtrs, na.rm = TRUE)
            if (plotLimits$rtRange[1] == 0 || plotLimits$rtRange[1] > rtmin)
                plotLimits$rtRange[1] <- rtmin
            if (plotLimits$rtRange[2] < rtmax)
                plotLimits$rtRange[2] <- rtmax
        }

        plotLimits$maxInt <- max(plotLimits$maxInt, max(gTable[[grpi]]))
    }

    plotLimits$rtRange <- plotLimits$rtRange + c(-rtWindow, rtWindow)
    if (retMin)
        plotLimits$rtRange <- plotLimits$rtRange / 60

    if (is.null(title))
    {
        if (gCount == 1)
            title <- sprintf("Group '%s'\nrt: %.1f - m/z: %.4f", names(gTable)[1],
                             if (retMin) gInfo[1, "rts"] / 60 else gInfo[1, "rts"],
                             gInfo[1, "mzs"])
        else
            title <- sprintf("%d feature groups", gCount)
    }

    anaIndsToPlot <- sapply(gNames, function(grp)
    {
        anasWGroup <- names(EICs)[sapply(EICs, function(e) !is.null(e[[grp]]))]
        ret <- match(anasWGroup, anaInfo$analysis)
        if (onlyPresent)
            ret <- ret[gTable[[grp]][ret] != 0]
        return(ret)
    }, simplify = FALSE)

    rGroupsInPlot <- unique(anaInfo$group[unlist(anaIndsToPlot)])
    fGroupsInPlot <- gNames[gNames %in% EICFGroups]

    oldp <- par(no.readonly = TRUE)
    if (showLegend)
    {
        makeLegend <- function(x, y, ...)
        {
            texts <- if (colourBy == "rGroups") rGroupsInPlot else fGroupsInPlot
            return(legend(x, y, texts, col = EICColors[texts],
                          text.col = EICColors[texts], lty = 1,
                          xpd = NA, ncol = 1, cex = 0.75, bty = "n", ...))
        }

        plot.new()
        leg <- makeLegend(0, 0, plot = FALSE)
        lw <- (grconvertX(leg$rect$w, to = "ndc") - grconvertX(0, to = "ndc"))
        par(omd = c(0, 1 - lw, 0, 1), new = TRUE)
    }

    if (is.null(xlim))
        xlim <- plotLimits$rtRange
    if (is.null(ylim))
        ylim <- c(0, plotLimits$maxInt * 1.1)

    plot(0, type = "n", main = title, xlab = sprintf("Retention time (%s)", if (retMin) "min." else "sec."), ylab = "Intensity",
         xlim = xlim, ylim = ylim, ...)

    if (showProgress)
        prog <- openProgBar(0, gCount)

    for (grp in fGroupsInPlot)
    {
        grpi <- match(grp, fGroupsInPlot)

        fts <- rbindlist(sapply(anaInfo$analysis, function(ana)
        {
            fti <- ftind[[grpi]][match(ana, anaInfo$analysis)]
            if (fti == 0)
                return(data.table(ret = NA, retmin = NA, retmax = NA, mz = NA))
            return(fTable[[ana]][fti])
        }, simplify = F), fill = TRUE)

        rtRange <- c(min(fts[anaIndsToPlot[[grp]], retmin], na.rm = TRUE), max(fts[anaIndsToPlot[[grp]], retmax], na.rm = TRUE))
        # rtEICRange <- rtRange + c(-rtWindow, rtWindow)

        if (retMin)
            rtRange <- rtRange / 60

        for (ana in names(EICs))
        {
            anai <- match(ana, anaInfo$analysis)
            if (is.na(anai))
                next

            if (onlyPresent && gTable[[grp]][anai] == 0)
                next

            EIC <- EICs[[ana]][[grp]]
            if (is.null(EIC))
                next

            if (colourBy == "rGroups")
                colInd <- anaInfo$group[anai]
            else if (colourBy == "fGroups")
                colInd <- grp
            else
                colInd <- 1

            points(if (retMin) EIC$time / 60 else EIC$time, EIC$intensity, type = "l", col = EICColors[colInd])

            if (showPeakArea && !is.na(fts[["mz"]][anai]))
            {
                EICFill <- EIC[numGTE(EIC$time, fts[anai, retmin]) & numLTE(EIC$time, fts[anai, retmax]), ]
                if (retMin)
                    EICFill$time <- EICFill$time / 60
                polygon(c(EICFill$time, rev(EICFill$time)), c(EICFill$intensity, rep(0, length(EICFill$intensity))),
                        col = fillColors[colInd], border = NA)
            }
        }

        if (showFGroupRect || !"none" %in% annotate)
        {
            intRange <- c(0, max(fts[anaIndsToPlot[[grp]], intensity], na.rm = TRUE))

            if (showFGroupRect)
                rect(rtRange[1], intRange[1], rtRange[2], intRange[2], border = "red", lty = "dotted")

            if (!"none" %in% annotate)
            {
                antxt <- character()
                rt <- mean(fts[anaIndsToPlot[[grp]], ret], na.rm = TRUE)
                if (retMin)
                    rt <- rt / 60

                if ("ret" %in% annotate)
                    antxt <- sprintf("%.1f", rt)
                if ("mz" %in% annotate)
                    antxt <- paste(antxt, sprintf("%.4f", gInfo[grp, "mzs"]), sep = "\n")

                if (nzchar(antxt))
                    text(rt, intRange[2] + ylim[2] * 0.02, antxt)
            }
        }

        if (showProgress)
            setTxtProgressBar(prog, grpi)
    }

    if (showLegend)
        makeLegend(par("usr")[2], par("usr")[4])

    if (showProgress)
    {
        setTxtProgressBar(prog, gCount)
        close(prog)
    }

    par(oldp)
})

setMethod("plotChromsHash", "featureGroups", function(obj, rtWindow = 30, mzExpWindow = 0.001, retMin = FALSE,
                                                      topMost = NULL, topMostByRGroup = FALSE, EICs = NULL,
                                                      showPeakArea = FALSE, showFGroupRect = TRUE,
                                                      title = NULL, colourBy = c("none", "rGroups", "fGroups"),
                                                      showLegend = TRUE, onlyPresent = TRUE,
                                                      annotate = c("none", "ret", "mz"), showProgress = FALSE,
                                                      xlim = NULL, ylim = NULL, ...)
{
    colourBy <- checkmate::matchArg(colourBy, c("none", "rGroups", "fGroups"))
    annotate <- checkmate::matchArg(annotate, c("none", "ret", "mz"), several.ok = TRUE)
    if ("none" %in% annotate)
        annotate <- "none"
    makeHash(allArgs(FALSE))
})

#' @describeIn featureGroups plots a Venn diagram (using
#'   \pkg{\link{VennDiagram}}) outlining unique and shared feature groups
#'   between up to five replicate groups.
#' @template plotvenn-ret
#' @export
setMethod("plotVenn", "featureGroups", function(obj, which = NULL, ...)
{
    rGroups <- replicateGroups(obj)
    if (is.null(which))
        which <- rGroups
    
    checkmate::assert(checkmate::checkSubset(which, rGroups, empty.ok = FALSE),
                      checkmate::checkList(which, "character", any.missing = FALSE),
                      .var.name = "which")

    fGroupsList <- lapply(which, function(w) replicateGroupFilter(obj, w, verbose = FALSE))
    
    if (is.list(which))
    {
        if (!checkmate::testNamed(which))
            names(which) <- sapply(which, paste0, collapse = "+", USE.NAMES = FALSE)
    }
    else
        names(which) <- which
    
    makeVennPlot(lapply(fGroupsList, names), names(which), lengths(fGroupsList),
                 function(obj1, obj2) intersect(obj1, obj2),
                 length, ...)
})

#' @describeIn featureGroups plots an UpSet diagram (using the
#'   \code{\link[UpSetR]{upset}} function) outlining unique and shared feature
#'   groups between given replicate groups.
#'
#' @param nsets,nintersects See \code{\link[UpSetR]{upset}}.
#'
#' @references \insertRef{Conway2017}{patRoon} \cr\cr
#'   \insertRef{Lex2014}{patRoon}
#'
#' @export
setMethod("plotUpSet", "featureGroups", function(obj, which = NULL, nsets = length(which),
                                                 nintersects = NA, ...)
{
    rGroups <- replicateGroups(obj)
    if (is.null(which))
        which <- rGroups

    ac <- checkmate::makeAssertCollection()
    checkmate::assertSubset(which, rGroups, empty.ok = FALSE, add = ac)
    checkmate::assertCount(nsets, positive = TRUE)
    checkmate::assertCount(nintersects, positive = TRUE, na.ok = TRUE)
    checkmate::reportAssertions(ac)

    if (length(obj) == 0)
        stop("Can't plot empty feature groups object")

    obj <- replicateGroupFilter(obj, which, verbose = FALSE)

    gt <- as.data.table(obj, average = TRUE)
    gt <- gt[, which, with = FALSE] # isolate relevant columns
    gt[, (which) := lapply(.SD, function(x) as.integer(x > 0))]

    if (sum(sapply(gt[, which, with = FALSE], function(x) any(x>0))) < 2)
        stop("Need at least two replicate groups with non-zero intensities")

    UpSetR::upset(gt, nsets = nsets, nintersects = nintersects, ...)
})

#' @describeIn featureGroups Obtain a subset with unique feature groups
#'   present in one or more specified replicate group(s).
#' @param relativeTo A character vector with replicate groups that should be
#'   used for unique comparison. If \code{NULL} then all replicate groups are
#'   used for comparison. Replicate groups specified in \code{which} are
#'   ignored.
#' @param outer If \code{TRUE} then only feature groups are kept which do not
#'   overlap between the specified replicate groups for the \code{which}
#'   parameter.
#' @export
setMethod("unique", "featureGroups", function(x, which, relativeTo = NULL, outer = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCharacter(which, min.len = 1, min.chars = 1, any.missing = FALSE, add = ac)
    checkmate::assertSubset(which, replicateGroups(x), empty.ok = FALSE, add = ac)
    checkmate::assertCharacter(relativeTo, min.len = 1, min.chars = 1, any.missing = FALSE,
                               null.ok = TRUE, add = ac)
    checkmate::assertFlag(outer, add = ac)
    checkmate::reportAssertions(ac)

    if (is.null(relativeTo))
        relativeTo <- replicateGroups(x)
    else
    {
        relativeTo <- union(which, relativeTo)
        x <- replicateGroupFilter(x, relativeTo, verbose = FALSE)
    }

    anaInfo <- analysisInfo(x)
    rGroups <- unique(anaInfo$group)

    if (all(rGroups %in% which) && !outer)
        return(x) # nothing to do...

    # Split by selected and other replicate groups
    selFGroups <- replicateGroupFilter(x, which, verbose = FALSE)
    
    otherFGNames <- setdiff(relativeTo, which)
    # NOTE: check length here to avoid ending up with an empty 'otherFGroups'
    # below, which yields warnings with XCMS.
    if (length(otherFGNames) > 0)
    {
        # pick out all feature groups NOT present in others
        otherFGroups <- replicateGroupFilter(x, otherFGNames, verbose = FALSE)
        ret <- selFGroups[, setdiff(names(selFGroups), names(otherFGroups))]
    }
    else
        ret <- selFGroups
    
    # remove all that is in at least 2 replicate groups
    if (outer && length(which) > 1)
        ret <- minReplicatesFilter(ret, absThreshold = 2, negate = TRUE, verbose = FALSE)

    return(ret)
})

#' @describeIn featureGroups Obtain a subset with feature groups that overlap
#'   between a set of specified replicate group(s).
#' @param exclusive If \code{TRUE} then all feature groups are removed that are
#'   not unique to the given replicate groups.
#' @aliases overlap
#' @export
setMethod("overlap", "featureGroups", function(fGroups, which, exclusive)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCharacter(which, min.len = 2, min.chars = 1, any.missing = FALSE, add = ac)
    checkmate::assertSubset(which, replicateGroups(fGroups), empty.ok = FALSE, add = ac)
    checkmate::assertFlag(exclusive, add = ac)
    checkmate::reportAssertions(ac)

    anaInfo <- analysisInfo(fGroups)
    rGroups <- unique(anaInfo$group)

    if (length(which) < 2 || length(fGroups) == 0)
        return(fGroups) # nothing to do...

    if (exclusive)
        ret <- unique(fGroups, which = which)
    else
        ret <- replicateGroupFilter(fGroups, which, verbose = FALSE)

    ret <- minReplicatesFilter(ret, relThreshold = 1, verbose = FALSE)

    return(ret)
})

setMethod("calculatePeakQualities", "featureGroups", function(fGroups, flatnessFactor)
{
    EICs <- getEICsForFGroups(fGroups, 0, 0, NULL, FALSE, TRUE)
    ftind <- groupFeatIndex(fGroups)
    anas <- analyses(fGroups)
    gNames <- names(fGroups)

    checkPackage("MetaClean")
    
    # HACK HACK HACK: MetaClean::calculateGaussianSimilarity needs to have
    # xcms::SSgauss attached
    # based on https://stackoverflow.com/a/36611896
    withr::local_environment(list(SSgauss = xcms::SSgauss))
        
    featQualities <- list(
        ApexBoundraryRatio = list(func = MetaClean::calculateApexMaxBoundaryRatio, HQ = "LV", range = c(0, 1)),
        FWHM2Base = list(func = MetaClean::calculateFWHM, HQ = "HV", range = c(0, 1)),
        Jaggedness = list(func = MetaClean::calculateJaggedness, HQ = "LV", range = Inf),
        Modality = list(func = MetaClean::calculateModality, HQ = "LV", range = Inf),
        Symmetry = list(func = MetaClean::calculateSymmetry, HQ = "HV", range = c(-1, 1)),
        GaussianSimilarity = list(func = MetaClean::calculateGaussianSimilarity, HQ = "HV", range = c(0, 1)),
        Sharpness = list(func = MetaClean::calculateSharpness, HQ = "HV", range = Inf),
        TPASR = list(func = MetaClean::calculateTPASR, HQ = "LV", range = Inf),
        ZigZag = list(func = MetaClean::calculateZigZagIndex, HQ = "LV", range = Inf)
    )
    groupQualities <- list(
        ElutionShift = list(func = MetaClean::calculateElutionShift, HQ = "LV", range = Inf),
        RetentionTimeCorrelation = list(func = MetaClean::calculateRetentionTimeConsistency, HQ = "LV", range = Inf)
    )
    
    calcFeatQualities <- function(ret, retmin, retmax, intensity, EIC)
    {
        args <- list(c(rt = ret, rtmin = retmin, rtmax = retmax, maxo = intensity), as.matrix(EIC))
        return(sapply(names(featQualities), function(q)
        {
            a <- args
            if (q %in% c("Jaggedness", "Modality"))
                a <- c(a, flatnessFactor)
            return(do.call(featQualities[[q]]$func, a))
        }, simplify = FALSE))
    }
    
    for (ana in names(EICs))
    {
        feat <- copy(featureTable(fGroups)[[ana]])
        anai <- match(ana, anas)
        featInds <- unlist(ftind[anai])
        groups <- gNames[featInds != 0]
        featInds <- featInds[featInds != 0]
        feat[featInds, (names(featQualities)) := rbindlist(Map(calcFeatQualities, ret, retmin, retmax,
                                                               intensity, EICs[[ana]][groups]))]
        fGroups@features@features[[ana]] <- feat
    }
    
    grpScores <- rbindlist(lapply(gNames, function(grp)
    {
        featInds <- ftind[[grp]]
        doAna <- anas[featInds != 0]
        featInds <- featInds[featInds != 0]
        fList <- rbindlist(Map(doAna, featInds, f = function(ana, row) fGroups@features[[ana]][row]))
        avgFeatScores <- sapply(names(featQualities), function(q)
        {
            if (all(is.na(fList[[q]])))
                return(NA_real_)
            # UNDONE: allow to specify average function?
            return(mean(fList[[q]], na.rm = TRUE))
        }, simplify = FALSE)
        
        pdata <- lapply(seq_len(nrow(fList)), function(fti) list(rtmin = fList$retmin[fti],
                                                                 rtmax = fList$retmax[fti]))
        eic <- lapply(doAna, function(a) EICs[[a]][[grp]])
        
        grpQualityScores <- sapply(lapply(groupQualities, "[[", "func"), do.call, list(pdata, eic), simplify = FALSE)
        
        return(c(avgFeatScores, grpQualityScores))
    }))
    
    # normalize, invert if necessary to get low (worst) to high (best) order
    fixGroupQuality <- function(quality, values)
    {
        qi <- if (quality %in% names(featQualities)) featQualities[[quality]] else groupQualities[[quality]]
        if (all(is.finite(qi$range)))
        {
            if (!isTRUE(all.equal(qi$range, c(0, 1)))) # no need to normalize 0-1
                values <- normalize(values, minMax = qi$range[1] < 0, xrange = qi$range)
        }
        else
            values <- normalize(values, TRUE)
        
        if (qi$HQ == "LV")
            values <- 1 - values
        
        return(values)
    }
    grpScores[, (names(grpScores)) := Map(fixGroupQuality, names(.SD), .SD)]
    grpScores[, Quality := rowSums(.SD, na.rm = TRUE)]
    
    fGroups@groupInfo <- cbind(fGroups@groupInfo, grpScores)
    
    return(fGroups)
})

#' @templateVar func groupFeatures
#' @templateVar what group features
#' @templateVar ex1 groupFeaturesOpenMS
#' @templateVar ex2 groupFeaturesXCMS3
#' @templateVar algos openms,xcms,xcms3
#' @template generic-algo
#'
#' @rdname feature-grouping
#' @aliases groupFeatures
#' @export
setMethod("groupFeatures", "features", function(feat, algorithm, ..., verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertChoice(algorithm, c("openms", "xcms", "xcms3"), add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    f <- switch(algorithm,
                openms = groupFeaturesOpenMS,
                xcms = groupFeaturesXCMS,
                xcms3 = groupFeaturesXCMS3)

    f(feat, ..., verbose = verbose)
})

#' @details \code{importFeatureGroups} is a generic function to import feature
#'   groups produced by other software. The actual functionality is provided by
#'   specific functions such as \code{importFeatureGroupsBrukerPA} and
#'   \code{importFeatureGroupsEnviMass}.
#'
#' @param path The path that should be used for importing. For
#'   \code{importFeatureGroupsBrukerPA} an exported 'bucket table' \file{.txt}
#'   file from Bruker ProfileAnalysis, for \code{importFeatureGroupsBrukerTASQ}
#'   an exported global result table (converted to \file{.csv}) and for
#'   \code{importFeatureGroupsEnviMass} the path of the enviMass project.
#'
#' @rdname feature-grouping
#' @export
importFeatureGroups <- function(path, type, ...)
{
    f <- switch(type,
                brukerpa = importFeatureGroupsBrukerPA,
                brukertasq = importFeatureGroupsBrukerTASQ,
                envimass = importFeatureGroupsEnviMass,
                stop("Invalid algorithm! Should be: brukerpa, brukertasq or envimass"))

    f(path, ...)
}
