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
#'   groups \samp{"influent-A/B"} against those in \samp{"effluent-A/B"}. The
#'   names of the list are used for labelling in the plot, and will be made
#'   automatically if not specified.
#' @param colourBy Sets the automatic colour selection: \code{"none"} for a
#'   single colour or \code{"rGroups"}/\code{"fGroups"} for a distinct colour
#'   per replicate/feature group.
#' @param showLegend If \code{TRUE} a legend will be shown with either replicate
#'   groups (\code{colourBy == "rGroups"}) or feature groups (\code{colourBy ==
#'   "fGroups"}, only for \code{plotChroms}). If \code{colourBy} is
#'   \code{"none"} no legend will be shown.
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
#' @slot annotations A (\code{\link{data.table}}) with adduct/isotope
#'   annotations for each group (see \code{\link{mergeIons}}).
#'
#' @templateVar class featureGroups
#' @template class-hierarchy
#'
#' @export
featureGroups <- setClass("featureGroups",
                          slots = c(groups = "data.table", analysisInfo = "data.frame", groupInfo = "data.frame",
                                    features = "features", ftindex = "data.table", groupQualities = "data.table",
                                    groupScores = "data.table", annotations = "data.table"),
                          contains = c("VIRTUAL", "workflowStep"))

setMethod("initialize", "featureGroups", function(.Object, ...)
{
    args <- list(...)

    # data.table's don't seem to initialize well (gives error that slot is init as list)
    if (is.null(args[["groups"]]))
        args$groups <- data.table()
    if (is.null(args[["ftindex"]]))
        args$ftindex <- data.table()
    if (is.null(args[["groupQualities"]]))
        args$groupQualities <- data.table()
    if (is.null(args[["groupScores"]]))
        args$groupScores <- data.table()
    if (is.null(args[["annotations"]]))
        args$annotations <- data.table()
    
    .Object <- do.call(callNextMethod, c(list(.Object), args))
    
    if (nrow(.Object@ftindex) > 0)
    {
        ftitr <- transpose(.Object@ftindex)
        gNames <- names(.Object)
        .Object@features@features <- Map(.Object@features@features, ftitr, f = function(feat, inds)
        {
            wh <- which(inds != 0)
            feat <- copy(feat)
            feat[inds[wh], group := gNames[wh]][]
            return(feat)
        })
        # remove unassigned features (eg in case the grouping algorithm already did some cleanup)
        oldfn <- length(.Object@features)
        .Object@features <- delete(getFeatures(.Object), j = function(ft, ...) is.na(ft$group))
        if (oldfn != length(.Object@features))
            .Object <- reGenerateFTIndex(.Object)
    }
    
    return(.Object)
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

setReplaceMethod("featureTable", "featureGroups", function(obj, value)
{
    featureTable(obj@features) <- value
    return(obj)
})

#' @describeIn featureGroups Accessor for \code{features} slot.
#' @export
setMethod("getFeatures", "featureGroups", function(obj) obj@features)

#' @describeIn featureGroups Accessor for \code{ftindex} slot.
#' @aliases groupFeatIndex
#' @export
setMethod("groupFeatIndex", "featureGroups", function(fGroups) fGroups@ftindex)

#' @describeIn featureGroups Accessor for \code{groupQualities} slot.
#' @aliases groupQualities
#' @export
setMethod("groupQualities", "featureGroups", function(fGroups) fGroups@groupQualities)

#' @describeIn featureGroups Accessor for \code{groupScores} slot.
#' @aliases groupScores
#' @export
setMethod("groupScores", "featureGroups", function(fGroups) fGroups@groupScores)

#' @describeIn featureGroups Accessor for \code{annotations} slot.
#' @aliases annotations
#' @export
setMethod("annotations", "featureGroups", function(fGroups) fGroups@annotations)

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
        x <- delete(x, setdiff(analyses(x), i))
    }

    if (!missing(j))
    {
        j <- assertSubsetArgAndToChr(j, names(x))
        x <- delete(x, j = setdiff(names(x), j))
    }

    return(x)
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

#' @export
setMethod("delete", "featureGroups", function(obj, i = NULL, j = NULL, ...)
{
    anas <- analyses(obj)
    gNames <- names(obj)
    
    ac <- checkmate::makeAssertCollection()
    i <- assertDeleteArgAndToChr(i, anas, add = ac)
    if (is.data.table(j))
    {
        # checkmate::assertDataTable(j, types = c("logical", "numeric"), col.names = "unique", add = ac)
        checkmate::assertNames(names(j), identical.to = gNames, what = "colnames", add = ac)
    }
    else if (!is.function(j))
        j <- assertDeleteArgAndToChr(j, gNames, add = ac)
    checkmate::reportAssertions(ac)
    
    # cases:
    # i = vector; j = NULL: subset analyses
    # i = NULL; j = vector: subset groups
    # i = vector; j = vector: remove the same features from analyses i in groups j
    # i = NULL/vector; j = function: use supplied function to remove specific features from each group

    if (length(i) == 0 || length(j) == 0 || length(obj) == 0)
        return(obj) # nothing to remove...
    
    ftind <- groupFeatIndex(obj)
    gTable <- groupTable(obj)
    jByIndex <- !is.function(j) && !is.data.table(j)
    isAnaSubSet <- jByIndex && setequal(j, gNames)
    isGrpSubSet <- jByIndex && setequal(i, anas)
    
    # remove features first
    if (isAnaSubSet)
        obj@features <- delete(getFeatures(obj), i = i)
    else if (jByIndex)
        obj@features <- delete(getFeatures(obj), i = i, j = function(ft, ...) which(ft$group %chin% j))
    else
    {
        if (is.data.table(j))
        {
            featsToRemove <- setnames(transpose(j), i)
            set(featsToRemove, j = "group", value = names(j))
        }
        else
        {
            gt <- gTable[chmatch(i, anas)]
            
            featsToRemove <- Map(j, gt, gNames, MoreArgs = list(...))
            # equalize lengths
            ol <- length(i)
            featsToRemove <- lapply(featsToRemove, function(x)
            {
                # use as.vector(.., "list") as it's a bit faster than as.list
                if (is.logical(x))
                    return(as.vector(rep(x, length.out = ol), "list"))
                if (is.numeric(x))
                    return(as.vector(seq_along(anas) %in% x, "list"))
                return(as.vector(anas %chin% x, "list"))
            })
            featsToRemove <- setnames(rbindlist(featsToRemove), i)
            set(featsToRemove, j = "group", value = gNames)
        }
        
        obj@features <- delete(getFeatures(obj), i = i, j = function(ft, ana)
        {
            return(featsToRemove[chmatch(ft$group, group), ana, with = FALSE][[1]] == TRUE)
        })
    }
    
    if (length(getFeatures(obj)) == 0)
    {
        # all features were removed, just clear out slots
        obj@groups <- obj@ftindex <- obj@groupQualities <- obj@groupScores <- data.table()
        obj@groupInfo <- obj@groupInfo[FALSE, ]
        obj@analysisInfo <- obj@analysisInfo[FALSE, ]
    }
    else
    {
        # remove analyses
        removedAnas <- setdiff(anas, analyses(getFeatures(obj)))
        if (length(removedAnas) > 0)
        {
            ainds <- chmatch(removedAnas, anas)
            if (length(obj) > 0)
            {
                obj@groups <- obj@groups[-ainds]
                obj@ftindex <- obj@ftindex[-ainds]
            }
            obj@analysisInfo <- obj@analysisInfo[-ainds, , drop = FALSE]
        }
        
        # remove deleted and empty groups
        removedGroups <- character()
        if (isGrpSubSet)
            removedGroups <- j
        else
            removedGroups <- setdiff(gNames, unique(unlist(lapply(featureTable(obj), "[[", "group"))))
        if (length(removedGroups) > 0)
        {
            ginds <- chmatch(removedGroups, gNames)
            if (length(obj) > 0)
            {
                obj@groups <- obj@groups[, -ginds, with = FALSE]
                obj@ftindex <- obj@ftindex[, -ginds, with = FALSE]
            }
            obj@groupInfo <- obj@groupInfo[-ginds, ]
            if (hasFGroupScores(obj))
            {
                obj@groupQualities <- setkey(obj@groupQualities[names(obj@groups)], "group")
                obj@groupScores <- setkey(obj@groupScores[names(obj@groups)], "group")
            }
        }
        
        if (!isAnaSubSet)
        {
            # UNDONE: can we skip updating things based on i/j?
            
            # re-generate feat index table by matching group names
            obj@ftindex <- reGenerateFTIndex(obj)
            
            # update group intensities: zero missing features
            ftind <- groupFeatIndex(obj) # update var
            gNames <- names(obj) # update var
            # NOTE: if j is a function it's assumed that all groups are affected
            affectedGrps <- if (jByIndex) intersect(j, gNames) else gNames
            obj@groups <- copy(obj@groups)
            # NOTE: assignment with by seems to be the fastest, as it allows some DT optimizations apparently...
            obj@groups[, (affectedGrps) := lapply(.SD, function(x) x), by = rep(1, nrow(obj@groups)),
                       .SDcols = affectedGrps]
        }
    }
    
    return(obj)
})

# UNDONE: make this public?
setMethod("removeEmptyAnalyses", "featureGroups", function(fGroups)
{
    if (length(fGroups) > 0)
    {
        trGT <- transpose(groupTable(fGroups))

        empty <- trGT[, sapply(.SD, sum) == 0]
        if (any(empty))
            fGroups <- delete(fGroups, empty)
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

#' @describeIn featureGroups Exports feature groups to a \file{.csv} file that
#'   is readable to Bruker ProfileAnalysis (a 'bucket table'), Bruker TASQ (an
#'   analyte database) or that is suitable as input for the \verb{Targeted peak
#'   detection} functionality of \href{http://mzmine.github.io/}{MZmine}.
#' @param out The destination file for the exported data.
#' @export
setMethod("export", "featureGroups", function(obj, type, out)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertChoice(type, c("brukerpa", "brukertasq", "mzmine"), add = ac)
    checkmate::assertPathForOutput(out, overwrite = TRUE, add = ac) # NOTE: assert doesn't work on Windows...
    checkmate::reportAssertions(ac)

    if (length(obj) == 0)
        stop("Cannot export empty feature groups object")

    if (type == "brukerpa")
    {
        # UNDONE: do we need this?
        #files <- sapply(bucketInfo$fInfo$analysis, function(f) file.path(bucketInfo$dataPath, paste0(f, ".d")), USE.NAMES = F)
        files <- obj@analysisInfo$analysis

        # col.names: if NA an empty initial column is added
        write.table(obj@groups, out, na = "", sep = "\t", quote = FALSE, row.names = files, col.names = NA)
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

        df <- data.frame(matrix(ncol=length(hdr), nrow=ncol(obj@groups)))
        colnames(df) <- hdr

        df["name"] <- colnames(obj@groups)
        df["m/z"] <- obj@groupInfo$mzs
        df["rt"] <- obj@groupInfo$rts / 60

        write.csv(df, out, row.names = FALSE, na = "")
    }
    else if (type == "mzmine")
    {
        df <- obj@groupInfo
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
#' @param normFunc Function that should be used for normalization of data. The
#'   function is called for all intensities/areas of a feature group and these
#'   quantities are divided by the result of the function call. For example,
#'   when \code{\link{max}} is used normalized intensities will be between zero
#'   and one. If all quantities are zero then the function will not be called.
#'   Set to \code{NULL} to perform no normalization.
#' @export
setMethod("as.data.table", "featureGroups", function(x, average = FALSE, areas = FALSE, features = FALSE,
                                                     qualities = FALSE, regression = FALSE, normFunc = NULL)
{
    # NOTE: keep args in sync with as.data.table() method for featureGroupsSet
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(average, add = ac)
    checkmate::assertFlag(areas, add = ac)
    checkmate::assertFlag(features, add = ac)
    checkmate::assertFlag(regression, add = ac)
    checkmate::assertFunction(normFunc, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    checkmate::assert(checkmate::checkFALSE(qualities),
                      checkmate::checkChoice(qualities, c("quality", "score", "both")),
                      .var.name = "qualities")
    
    if (length(x) == 0)
        return(data.table(mz = numeric(), ret = numeric(), group = character()))

    if (features && average && regression)
        stop("Cannot add regression data for averaged features.")
    if (features && average && !is.null(normFunc))
        stop("Cannot normalize data for averaged features.")
    
    anaInfo <- analysisInfo(x)
    gNames <- names(x)
    gInfo <- groupInfo(x)
    doConc <- regression && !is.null(anaInfo[["conc"]]) && sum(!is.na(anaInfo[["conc"]]) > 1)
    addQualities <- !isFALSE(qualities) && qualities %in% c("both", "quality") && hasFGroupScores(x)
    addScores <- !isFALSE(qualities) && qualities %in% c("both", "score") && hasFGroupScores(x)

    if (regression && is.null(anaInfo[["conc"]]))
        warning("No concentration information specified in the analysis information (i.e. conc column, see ?`analysis-information`)")

    if (features)
    {
        ftindex <- groupFeatIndex(x)
        fTable <- featureTable(x)
        
        ret <- rbindlist(fTable, idcol = "analysis")
        setorder(ret, "group")
        if (doConc)
            ret[, conc := anaInfo$conc[match(analysis, anaInfo$analysis)]]

        if (!is.null(ret[["adduct"]]))
            ret[, adduct := NULL] # we already include group annotations below

        if (addQualities)
        {
            gq <- groupQualities(x)[ret$group, -"group"]
            ret[, (paste0("group_", names(gq))) := gq]
        }
        else if (hasFGroupScores(x))
            ret[, (intersect(featureQualityNames(), names(ret))) := NULL]
        if (addScores)
        {
            gs <- groupScores(x)[ret$group, -"group"]
            ret[, (paste0("group_", names(gs))) := gs]
        }
        else if (hasFGroupScores(x))
            ret[, (intersect(c(featureScoreNames(), "totalScore"), names(ret))) := NULL]
        
        if (average)
        {
            ret <- ret[, -c("isocount", "analysis", "ID")]
            numCols <- setdiff(names(ret), c("group"))
            ret[, (numCols) := lapply(.SD, mean), .SDcols = numCols, by = "group"]
            ret <- unique(ret, by = "group")
        }
        else
        {
            if (!is.null(normFunc))
            {
                ret[, c("area", "intensity") := .(if (all(area == 0)) area else area / normFunc(area),
                                                  if (all(intensity == 0)) intensity else intensity / normFunc(intensity)),
                    by = "group"]
            }
            
            doConc <- doConc && nrow(anaInfo) > 1
            if (doConc)
            {
                ret[, c("RSQ", "intercept", "slope") := {
                    notna <- !is.na(conc)
                    if (!any(notna))
                        NA_real_
                    else
                    {
                        suppressWarnings(reg <- summary(lm(intensity[notna] ~ conc[notna])))
                        slope <- if (nrow(reg[["coefficients"]]) > 1) reg[["coefficients"]][2, 1] else NA_real_
                        list(reg[["r.squared"]], reg[["coefficients"]][1, 1], slope)
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
        
        if (!is.null(normFunc))
        {
            normv <- gTable[, lapply(.SD, normFunc)]
            gTable <- copy(gTable)
            for (g in seq_along(gTable))
            {
                if (!all(gTable[[g]] == 0))
                    set(gTable, j = g, value = gTable[[g]] / normv[[g]])
            }
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
        
        if (addQualities)
            ret <- cbind(ret, groupQualities(x)[ret$group, -"group"])
        if (addScores)
            ret <- cbind(ret, groupScores(x)[ret$group, -"group"])
    }

    annTable <- annotations(x)
    if (nrow(ret) > 0 && nrow(annTable) > 0)
        ret <- merge(ret, annTable)
    
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
        EICs <- getEICsForFGroups(obj, rtWindow, mzExpWindow, topMost, topMostByRGroup, onlyPresent)
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
        # NOTE: plotChroms() for sets override default
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
    args <- allArgs(FALSE)
    makeHash(args[setdiff(names(args), "obj")], featureTable(obj), groupInfo(obj), analysisInfo(obj))
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

setMethod("calculatePeakQualities", "featureGroups", function(obj, weights, flatnessFactor, avgFunc = mean,
                                                              parallel = TRUE)
{
    allScores <- c(featureScoreNames(), featureGroupScoreNames())
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumeric(weights, finite = TRUE, any.missing = FALSE, min.len = 1, names = "unique",
                             null.ok = TRUE, add = ac)
    if (!is.null(weights))
        checkmate::assertNames(names(weights), subset.of = allScores, add = ac)
    checkmate::assertNumber(flatnessFactor, add = ac)
    checkmate::assertFunction(avgFunc, add = ac)
    checkmate::assertFlag(parallel, add = ac)
    checkmate::reportAssertions(ac)
    
    hash <- makeHash(obj, weights, flatnessFactor, avgFunc)
    cd <- loadCacheData("calculatePeakQualities", hash)
    if (!is.null(cd))
        return(cd)
    
    checkPackage("MetaClean")
    
    if (length(obj) == 0)
        return(obj)
    
    if (!is.null(weights))
    {
        weights[setdiff(allScores, names(weights))] <- 1
        weights <- weights[allScores]
    }
    
    w <- if (!is.null(weights) && any(names(weights) %in% featureScoreNames())) weights[featureScoreNames()] else NULL
    obj@features <- calculatePeakQualities(getFeatures(obj), weights = w, flatnessFactor = flatnessFactor,
                                           parallel = parallel)
    
    ftind <- groupFeatIndex(obj)
    anas <- analyses(obj)
    gNames <- names(obj)
    gCount <- length(obj)
    EICs <- getEICsForFGroups(obj, 0, 0.001, NULL, FALSE, TRUE)
    fgQualities <- featureGroupQualities()
    
    printf("Calculating group peak qualities and scores...\n")
    prog <- openProgBar(0, gCount)
    
    groupQualitiesScores <- rbindlist(sapply(gNames, function(grp)
    {
        featInds <- ftind[[grp]]
        doAna <- anas[featInds != 0]
        featInds <- featInds[featInds != 0]
        fList <- rbindlist(Map(doAna, featInds, f = function(ana, row) obj@features[[ana]][row]))
        featAvgs <- sapply(c(featureQualityNames(), featureScoreNames()), function(q)
        {
            if (all(is.na(fList[[q]])))
                return(NA_real_)
            return(avgFunc(fList[[q]], na.rm = TRUE))
        }, simplify = FALSE)
        
        pdata <- lapply(seq_len(nrow(fList)), function(fti) list(rtmin = fList$retmin[fti],
                                                                 rtmax = fList$retmax[fti]))
        # NOTE: MetaClean expects EIC matrices
        eic <- lapply(doAna, function(a) as.matrix(EICs[[a]][[grp]]))
        gq <- sapply(lapply(fgQualities, "[[", "func"), do.call, list(pdata, eic), simplify = FALSE)
        
        setTxtProgressBar(prog, match(grp, gNames))
        
        return(c(featAvgs, gq))
    }, simplify = FALSE), idcol = "group")
    
    groupQualitiesScores[, (featureGroupScoreNames()) := Map(scoreFeatQuality, fgQualities, .SD),
                         .SDcols = featureGroupQualityNames()]
    setkeyv(groupQualitiesScores, "group")
    
    obj@groupQualities <- groupQualitiesScores[, c("group", featureQualityNames(), featureGroupQualityNames()), with = FALSE]
    obj@groupScores <- groupQualitiesScores[, c("group", allScores), with = FALSE]
    
    if (is.null(weights))
        obj@groupScores[, totalScore := rowSums(.SD, na.rm = TRUE), .SDcols = allScores][]
    else
    {
        wsc <- obj@groupScores[, allScores, with = FALSE]
        wsc[, (names(wsc)) := Map("*", .SD, weights)]
        set(obj@groupScores, j = "totalScore", value = rowSums(wsc, na.rm = TRUE))
    }
    
    setTxtProgressBar(prog, gCount)
    
    saveCacheData("calculatePeakQualities", obj, hash)
    
    return(obj)
})

#' @export
setMethod("mergeIons", "featureGroups", function(fGroups, components, prefAdduct, selectIsoBy = "mono")
{
    # UNDONE is intensity_rel a proper measure? ie does it allow comparison if
    # isotopes/adducts are taken from different analyses?
    # UNDONE: add logging to see what happens
    # UNDONE: include isotope charge?
    # UNDONE: does it make sense to keep isotope other than mono? can't use it for annotation
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(components, "components", add = ac)
    checkmate::assertChoice(selectIsoBy, c("mono", "intensity"))
    checkmate::reportAssertions(ac)
    
    prefAdduct <- as.character(checkAndToAdduct(prefAdduct))
    
    if (is.null(componentInfo(components)[["neutral_mass"]]))
        stop("No adduct/isotope information available in given components!")
    
    cTab <- as.data.table(components)
    cTab <- cTab[group %in% names(fGroups) & (!is.na(isonr) | !is.na(adduct_ion))]
    cTab[, remove := FALSE]
    cTab[!is.na(isonr), remove := {
        if (.N == 1)
            FALSE
        else if (selectIsoBy == "mono")
            isonr != 0
        else
            !numEQ(intensity_rel, max(intensity_rel))
    }, by = c("name", "isogroup")]
    
    cTab[!is.na(adduct_ion) & remove == FALSE, remove := {
        if (.N == 1)
            FALSE
        # UNDONE: allow >1 pref adducts?
        else if (prefAdduct %in% adduct_ion)
            adduct_ion != prefAdduct
        else # fall back to most intense
            !numEQ(intensity_rel, max(intensity_rel))
    }, by = "name"]
    
    # remove unwanted isotopes/adducts
    fGroups <- fGroups[, setdiff(names(fGroups), cTab[remove == TRUE]$group)]
    
    # annotate remaining
    
    cTabIso <- cTab[!is.na(isonr) & !remove]
    cTabAdd <- cTab[!is.na(adduct_ion) & !remove]
    
    if (nrow(cTabIso) == 0 && nrow(cTabAdd) == 0)
    {
        fGroups@annotations <- data.table()
        cat("No adduct or isotope annotations found!\n")
    }
    else
    {
        fGroups@annotations <- data.table(group = names(fGroups))

        if (nrow(cTabIso) > 0)
        {
            fGroups@annotations <- merge(fGroups@annotations, cTabIso[, c("group", "isonr"), with = FALSE],
                                         by = "group", all.x = TRUE)
            fGroups@annotations[is.na(isonr), isonr := 0]
        }
        else
            fGroups@annotations[, isonr := 0]
        
        if (nrow(cTabAdd) > 0)
        {
            setnames(cTabAdd, "adduct_ion", "adduct")
            fGroups@annotations <- merge(fGroups@annotations, cTabAdd[, c("group", "adduct"), with = FALSE],
                                         by = "group", all.x = TRUE)
            fGroups@annotations[is.na(adduct), adduct := prefAdduct]
        }
        else
            fGroups@annotations[, adduct := prefAdduct]

        fGroups@annotations[, neutralMass := groupInfo(fGroups)[group, "mzs"] -
                                sapply(lapply(adduct, as.adduct), adductMZDelta)]
                
        # retain correct order
        fGroups@annotations <- fGroups@annotations[match(names(fGroups), group)]
        
        printf("Removed %d feature groups detected as unwanted adducts/isotopes\n", sum(cTab$remove))
        printf("Annotated %d feature groups with isotope information\n", nrow(cTabIso))
        printf("\tRemaining %d feature groups set as default isotope 0\n", length(fGroups) - nrow(cTabIso))
        printf("Annotated %d feature groups with adducts\n", nrow(cTabAdd))
        printf("\tRemaining %d feature groups set as default adduct %s\n", length(fGroups) - nrow(cTabAdd), prefAdduct)
    }
    
    return(fGroups)
})

#' @templateVar func groupFeatures
#' @templateVar what group features
#' @templateVar ex1 groupFeaturesOpenMS
#' @templateVar ex2 groupFeaturesXCMS3
#' @templateVar algos openms,xcms,xcms3,kpic2
#' @template generic-algo
#'
#' @rdname feature-grouping
#' @aliases groupFeatures
#' @export
setMethod("groupFeatures", "features", function(feat, algorithm, ..., verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertChoice(algorithm, c("openms", "xcms", "xcms3", "kpic2"), add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    f <- switch(algorithm,
                openms = groupFeaturesOpenMS,
                xcms = groupFeaturesXCMS,
                xcms3 = groupFeaturesXCMS3,
                kpic2 = groupFeaturesKPIC2)

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
