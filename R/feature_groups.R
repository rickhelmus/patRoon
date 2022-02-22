#' @include main.R
#' @include features.R
#' @include workflow-step.R
NULL

#' Base class for grouped features.
#'
#' This class holds all the information for grouped features.
#'
#' The \code{featureGroup} class is the workhorse of \pkg{patRoon}: almost all functionality operate on its instantiated
#' objects. The class holds all information from grouped features (obtained from \code{\link{features}}). This class
#' itself is \code{virtual}, hence, objects are not created directly from it. Instead, 'feature groupers' such as
#' \code{\link{groupFeaturesXCMS}} return a \code{featureGroups} derived object after performing the actual grouping of
#' features across analyses.
#'
#' @param fGroups,obj,x,object \code{featureGroups} object to be accessed.
#' @param rGroups For \code{[}: An optional \code{character} vector: if specified only keep results for the given
#'   replicate groups (equivalent to the \code{rGroups} argument to \code{\link[=filter,featureGroups-method]{filter}}).
#' @param \dots For the \code{"["} operator: ignored.
#'
#'   For \code{delete}: passed to the function specified as \code{j}.
#'
#'   \setsPassedArgs1{featureGroups}
#' @param average If \code{TRUE} then data within replicate groups are averaged.
#'
#'   For \code{as.data.table}: if \code{features=TRUE} other feature properties are also averaged.
#' @param averageFunc Function used for averaging. Only used when \code{average=TRUE} or \code{FCParams != NULL}.
#' @param areas If set to \code{TRUE} then areas are considered instead of peak intensities.
#'
#'   For \code{as.data.table}: ignored if \code{features=TRUE}, as areas of features are always reported.
#' @param which A character vector with replicate groups used for comparison.
#' @param FCParams A parameter list to calculate Fold change data. See \code{getFCParams} for more details. Set to
#'   \code{NULL} to not perform FC calculations.
#'
#' @templateVar seli analyses
#' @templateVar selOrderi analyses()
#' @templateVar selj feature groups
#' @templateVar selOrderj names()
#' @templateVar optionalji TRUE
#' @templateVar del TRUE
#' @templateVar deli analyses
#' @templateVar delj feature groups
#' @templateVar deljtype numeric index, logical or character
#' @templateVar delfwhat feature group
#' @templateVar delfa1 a vector of the group intensities
#' @templateVar delfa2 the group name
#' @templateVar delfr the analyses of the features in the group to be removed (same format as \code{i})
#' @templateVar dollarOpName feature group
#' @template sub_sel_del-args
#'
#' @slot groups Matrix (\code{\link{data.table}}) with intensities for each feature group (columns) per analysis (rows).
#'   Access with \code{groups} method.
#' @slot analysisInfo,features \link[=analysis-information]{Analysis info} and \code{\link{features}} class associated
#'   with this object. Access with \code{analysisInfo} and \code{featureTable} methods, respectively.
#' @slot groupInfo \code{data.frame} with retention time (\code{rts} column, in seconds) and \emph{m/z} (\code{mzs}
#'   column) for each feature group. Access with \code{groupInfo} method.
#' @slot ftindex Matrix (\code{\link{data.table}}) with feature indices for each feature group (columns) per analysis
#'   (rows). Each index corresponds to the row within the feature table of the analysis (see
#'   \code{\link{featureTable}}).
#' @slot groupQualities,groupScores A \code{\link{data.table}} with qualities/scores for each feature group (see the
#'   \code{calculatePeakQualities} method).
#' @slot annotations A \code{\link{data.table}} with adduct annotations for each group (see the \code{selectIons}
#'   method).
#'
#' @templateVar class featureGroups
#' @template class-hierarchy
#'
#' @seealso \code{\link{groupFeatures}} for generating feature groups, \link{feature-filtering} and
#'   \link{feature-plotting} for more advanced \code{featureGroups} methods.
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
    else
    {
        .Object@features@features <- lapply(.Object@features@features, function(feat)
        {
            feat <- copy(feat)
            feat[, group := character()]
        })
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
    fCount <- length(getFeatures(object)); gCount <- length(object)
    printf("Feature groups: %s (%d total)\n", getStrListWithMax(names(object), 6, ", "), gCount)
    printf("Features: %d (%.1f per group)\n", fCount, if (gCount > 0) fCount / gCount)
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
#' @export
setMethod("annotations", "featureGroups", function(obj) obj@annotations)

#' @describeIn featureGroups Returns a named \code{character} with adduct annotations assigned to each feature group (if
#'   available).
#' @export
setMethod("adducts", "featureGroups", function(obj)
{
    if (nrow(annotations(obj)) == 0)
        return(character())
    return(setNames(annotations(obj)$adduct, annotations(obj)$group))
})

#' @describeIn featureGroups Sets adduct annotations for feature groups.
#' @param value For \code{adducts<-}: A \code{character} with adduct annotations assigned to each feature group. The
#'   length should equal the number of feature groups. Can be named with feature group names to customize the assignment
#'   order.
#' @export
setReplaceMethod("adducts", "featureGroups", function(obj, value)
{
    checkmate::assertCharacter(value, min.chars = 1, any.missing = FALSE, len = length(obj))
    
    if (checkmate::testNamed(value))
    {
        checkmate::assertNames(names(value), permutation.of = names(obj), .var.name = "value")
        value <- value[names(obj)] # ensure correct order
    }
    else
        names(value) <- names(obj)
    
    obj@annotations <- updateAnnAdducts(annotations(obj), groupInfo(obj), value)
    
    return(obj)
})

#' @describeIn featureGroups Subset on analyses/feature groups.
#' @param results Optional argument. If specified only feature groups with results in the specified object are kept. The
#'   class of \code{results} should be \code{\link{featureAnnotations}} or \code{\link{components}}. Multiple objects
#'   can be specified in a \code{list}: in this case a feature group is kept if it has a result in \emph{any} of the
#'   objects (equivalent to the \code{results} argument to \code{\link[=filter,featureGroups-method]{filter}}).
#' @export
setMethod("[", c("featureGroups", "ANY", "ANY", "missing"), function(x, i, j, ..., rGroups, results, drop = TRUE)
{
    if (!missing(rGroups))
        x <- filter(x, rGroups = rGroups)
    if (!missing(results))
        x <- filter(x, results = results)
    
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

#' @templateVar where featureGroups
#' @templateVar what feature groups
#' @template delete
#' @export
setMethod("delete", "featureGroups", function(obj, i = NULL, j = NULL, ...)
{
    anas <- analyses(obj)
    gNames <- names(obj)
    iNULL <- is.null(i); jNULL <- is.null(j)
    
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
    # i = NULL/vector; j = function/data.table: remove specific features from each group (all analyses if i=NULL)

    if ((!iNULL && length(i) == 0) || (!jNULL && length(j) == 0))
        return(obj) # nothing to remove...
    
    ftind <- groupFeatIndex(obj)
    gTable <- groupTable(obj)
    jByIndex <- !is.function(j) && !is.data.table(j)
    isAnaSubSet <- jNULL
    isGrpSubSet <- jByIndex && iNULL
    
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
            obj@groupQualities <- obj@groupQualities[group %in% names(obj@groups)]
            obj@groupScores <- obj@groupScores[group %in% names(obj@groups)]
        }
        if (nrow(obj@annotations) > 0)
            obj@annotations <- obj@annotations[group %in% names(obj@groups)]
    }
    
    if (!isAnaSubSet)
    {
        # UNDONE: can we skip updating things based on i/j?
        
        # re-generate feat index table by matching group names
        obj <- reGenerateFTIndex(obj)
        
        # update group intensities: zero missing features
        ftind <- groupFeatIndex(obj) # update var
        gNames <- names(obj) # update var
        # NOTE: if j is a function it's assumed that all groups are affected
        affectedGrps <- if (jByIndex) intersect(j, gNames) else gNames
        obj@groups <- copy(obj@groups)
        # NOTE: assignment with by seems to be the fastest, as it allows some DT optimizations apparently...
        obj@groups[, (affectedGrps) := Map(.SD, affectedGrps, f = function(x, g) fifelse(ftind[[g]] == 0, 0, x)),
                   by = rep(1, nrow(obj@groups)), .SDcols = affectedGrps]
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

setMethod("averageGroups", "featureGroups", function(fGroups, areas, func)
{
    gTable <- copy(groupTable(fGroups, areas))
    if (nrow(gTable) == 0)
        return()

    gNames <- names(fGroups)
    anaInfo <- analysisInfo(fGroups)

    gTable[, sgroup := anaInfo$group]

    gTable[, (gNames) := lapply(.SD, function(v) { if (any(v > 0)) func(v[v>0]) else 0 }), by = sgroup, .SDcols = gNames]
    gTable <- unique(gTable, by = "sgroup")
    gTable[, sgroup := NULL]

    return(gTable[])
})

#' @describeIn featureGroups Exports feature groups to a \file{.csv} file that is readable to Bruker ProfileAnalysis (a
#'   'bucket table'), Bruker TASQ (an analyte database) or that is suitable as input for the \verb{Targeted peak
#'   detection} functionality of \href{http://mzmine.github.io/}{MZmine}.
#' @param type The export type: \code{"brukerpa"} (Bruker ProfileAnalysis), \code{"brukertasq"} (Bruker TASQ) or
#'   \code{"mzmine"} (MZmine).
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

#' @describeIn featureGroups Obtain a summary table (a \code{\link{data.table}}) with retention, \emph{m/z}, intensity
#'   and optionally other feature data.
#' @param features If \code{TRUE} then feature specific data will be added. If \code{average=TRUE} this data will be
#'   averaged for each feature group.
#' @param qualities Adds feature (group) qualities (\code{qualities="quality"}), scores (\code{qualities="score"}) or
#'   both (\code{qualities="both"}), if this data is available (\emph{i.e.} from \code{calculatePeakQualities}). If
#'   \code{qualities=FALSE} then nothing is reported.
#' @param regression Set to \code{TRUE} to add regression data for each feature group. For this a linear model is
#'   created (intensity/area [depending on \code{areas} argument] \emph{vs} concentration). The model concentrations
#'   (e.g. of a set of standards) is derived from the \code{conc} column of the \link[=analysis-information]{analysis
#'   information}. From this model the intercept, slope and R2 is added to the output. In addition, when
#'   \code{features=TRUE}, concentrations for each feature are added. Note that no regression information is added when
#'   no \code{conc} column is present in the analysis information or when less than two concentrations are specified
#'   (\emph{i.e.} the minimum amount).
#' @param normFunc Function that should be used for normalization of data. The function is called for all
#'   intensities/areas of a feature group and these quantities are divided by the result of the function call. For
#'   example, when \code{\link{max}} is used normalized intensities will be between zero and one. If all quantities are
#'   zero then the function will not be called. Set to \code{NULL} to perform no normalization.
#' @export
setMethod("as.data.table", "featureGroups", function(x, average = FALSE, areas = FALSE, features = FALSE,
                                                     qualities = FALSE, regression = FALSE, averageFunc = mean,
                                                     normFunc = NULL, FCParams = NULL)
{
    # NOTE: keep args in sync with as.data.table() method for featureGroupsSet
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(average, add = ac)
    checkmate::assertFlag(areas, add = ac)
    checkmate::assertFlag(features, add = ac)
    checkmate::assertFlag(regression, add = ac)
    checkmate::assertFunction(averageFunc, add = ac)
    checkmate::assertFunction(normFunc, null.ok = TRUE, add = ac)
    assertFCParams(FCParams, x, null.ok = TRUE, add = ac)
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
    if (features && !is.null(FCParams))
        stop("Cannot calculate fold-changes with features=TRUE")
    
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
            gq <- groupQualities(x)[match(ret$group, group), -"group"]
            ret[, (paste0("group_", names(gq))) := gq]
        }
        else if (hasFGroupScores(x))
            ret[, (intersect(featureQualityNames(group = FALSE), names(ret))) := NULL]
        if (addScores)
        {
            gs <- groupScores(x)[match(ret$group, group), -"group"]
            ret[, (paste0("group_", names(gs))) := gs]
        }
        else if (hasFGroupScores(x))
            ret[, (intersect(featureQualityNames(group = FALSE, scores = TRUE), names(ret))) := NULL]
        
        if (average)
        {
            ret <- ret[, -c("isocount", "analysis", "ID")]
            numCols <- setdiff(names(ret), c("group"))
            ret[, (numCols) := lapply(.SD, averageFunc), .SDcols = numCols, by = "group"]
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
        gTableAvg <- averageGroups(x, areas, func = averageFunc)
        gTableNonAvg <- groupTable(x, areas)

        if (!is.null(normFunc))
        {
            doNorm <- function(gt)
            {
                normv <- gt[, lapply(.SD, normFunc)]
                gt <- copy(gt)
                for (g in seq_along(gt))
                {
                    if (!all(gt[[g]] == 0))
                        set(gt, j = g, value = gt[[g]] / normv[[g]])
                }
                return(gt)
            }
            gTableAvg <- doNorm(gTableAvg); gTableNonAvg <- doNorm(gTableNonAvg)
        }
        
        if (average)
        {
            gTable <- gTableAvg
            snames <- unique(anaInfo$group)
            if (doConc)
                concs <- anaInfo[!duplicated(anaInfo$group), "conc"] # conc should be same for all replicates
        }
        else
        {
            gTable <- gTableNonAvg
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

        if (!is.null(FCParams))
        {
            calcFC <- function(x, y)
            {
                fixZeros <- function(x)
                {
                    zx <- which(x == 0)
                    if (FCParams$zeroMethod == "add")
                        x[zx] <- x[zx] + FCParams$zeroValue
                    else if (FCParams$zeroMethod == "fixed")
                        x[zx] <- FCParams$zeroValue
                    else # "omit"
                        x <- x[!zx]
                    return(x)                    
                }
                return(fixZeros(y) / fixZeros(x))
            }
            
            repInds <- match(FCParams$rGroups, replicateGroups(x))
            for (i in seq_along(gTableAvg))
                set(ret, i, "FC", do.call(calcFC, as.list(gTableAvg[[i]][repInds])))
            ret[, FC_log := log2(FC)]
            
            anaInds1 <- which(anaInfo$group %in% FCParams$rGroups[1])
            anaInds2 <- which(anaInfo$group %in% FCParams$rGroups[2])
            ret[, PV := mapply(gTableNonAvg[anaInds1, ], gTableNonAvg[anaInds2, ], FUN = FCParams$PVTestFunc)]
            ret[, PV := FCParams$PVAdjFunc(PV)]
            ret[, PV_log := -log10(PV)]
            
            isSignificant <- ret$PV < FCParams$thresholdPV
            ret[, classification := "insignificant"] # by default
            ret[isSignificant & numGTE(FC_log, FCParams$thresholdFC), classification := "increase"]
            ret[isSignificant & numLTE(FC_log, FCParams$thresholdFC), classification := "decrease"]
            ret[!isSignificant & numGTE(abs(FC_log), FCParams$thresholdFC), classification := "FC"]
            ret[isSignificant & numLTE(abs(FC_log), FCParams$thresholdFC), classification := "significant"]
        }
        
        ret[, c("group", "ret", "mz") := .(gNames, gInfo$rts, gInfo$mzs)]
        setcolorder(ret, c("group", "ret", "mz"))
        
        if (addQualities)
            ret <- cbind(ret, groupQualities(x)[match(ret$group, group), -"group"])
        if (addScores)
            ret <- cbind(ret, groupScores(x)[match(ret$group, group), -"group"])
    }

    annTable <- annotations(x)
    if (nrow(ret) > 0 && nrow(annTable) > 0)
        ret <- merge(ret, annTable, sort = FALSE)
    
    return(ret[])
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

#' @describeIn featureGroups Calculates peak and group qualities for all features and feature groups. The peak qualities
#'   (and scores) are calculated with the \link[=calculatePeakQualities,features-method]{features method of this
#'   function}, and subsequently averaged per feature group. Then, \pkg{MetaClean} is used to calculate the
#'   \verb{Elution Shift} and \verb{Retention Time Consistency} group quality metrics (see the \pkg{MetaClean}
#'   publication cited below for more details). Similarly to the \code{\link{features}} method, these metrics are scored
#'   by normalizing qualities among all groups and scaling them from \samp{0} (worst) to \samp{1} (best). The
#'   \verb{totalScore} for each group is then calculated as the weighted sum from all feature (group) scores. The
#'   \code{\link{getMCTrainData}} and \code{\link{predictCheckFeaturesSession}} functions can be used to train and apply
#'   Pass/Fail ML models from \pkg{MetaClean}.
#'
#' @inheritParams calculatePeakQualities,features-method
#' @param avgFunc The function used to average the peak qualities and scores for each feature group.
#'
#' @template parallel-arg
#' 
#' @references \insertRef{Chetnik2020}{patRoon}
#'
#' @return \code{calculatePeakQualities} returns a modified object amended with peak qualities and scores.
#'
#' @export
setMethod("calculatePeakQualities", "featureGroups", function(obj, weights, flatnessFactor, avgFunc = mean,
                                                              parallel = TRUE)
{
    allScores <- featureQualityNames(scores = TRUE, totScore = FALSE)
    
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
    
    fs <- featureQualityNames(group = FALSE, scores = TRUE, totScore = FALSE)
    w <- if (!is.null(weights) && any(names(weights) %in% fs)) weights[names(weights) %in% fs] else NULL
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
        featAvgs <- sapply(c(featureQualityNames(group = FALSE),
                             featureQualityNames(group = FALSE, scores = TRUE, totScore = FALSE)), function(q)
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
    
    groupQualitiesScores[, (featureQualityNames(feat = FALSE, scores = TRUE, totScore = FALSE)) :=
                             Map(scoreFeatQuality, fgQualities, .SD),
                         .SDcols = featureQualityNames(feat = FALSE)]
    
    obj@groupQualities <- groupQualitiesScores[, c("group", featureQualityNames()), with = FALSE]
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

#' @describeIn featureGroups uses \link[=generateComponents]{componentization} results to select feature groups with
#'   preferred adduct ion and/or isotope annotation. Typically, this means that only feature groups are kept if they are
#'   (de-)protonated adducts and are monoisotopic. The adduct annotation assignments for the selected feature groups are
#'   copied from the components to the \code{annotations} slot. If the adduct for a feature group is unknown, its
#'   annotation is defaulted to the 'preferred' adduct, and hence, the feature group will never be removed. Furthermore,
#'   if a component does not contain an annotation with the preferred adduct, the most intense feature group is selected
#'   instead. Similarly, if no isotope annotation is available, the feature group is assumed to be monoisotopic and thus
#'   not removed. An important advantage of \code{selectIons} is that it may considerably simplify your dataset.
#'   Furthermore, the adduct assignments allow formula/compound annotation steps later in the workflow to improve their
#'   annotation accuracy. On the other hand, it is important the componentization results are reliable. Hence, it is
#'   highly recommended that, prior to calling \code{selectIons}, the settings to \code{\link{generateComponents}} are
#'   optimized and its results are reviewed with \code{\link{checkComponents}}. Finally, the \code{adducts<-} method can
#'   be used to manually correct adduct assignments afterwards if necessary.
#'
#' @param components The \code{components} object that was generated for the given \code{featureGroups} object.
#'   Obviously, the components must be created with algorithms that support adduct/isotope annotations, such as those
#'   from \pkg{RAMClustR} and \pkg{cliqueMS}.
#' @param prefAdduct The 'preferred adduct' (see method description). This is often \code{"[M+H]+"} or \code{"[M-H]-"}.
#' @param onlyMonoIso Set to \code{TRUE} to only keep feature groups that were annotated as monoisotopic. Feature groups
#'   are never removed by this setting if no isotope annotations are available.
#' @param chargeMismatch Specifies how to deal with a mismatch in charge between adduct and isotope annotations. Valid
#'   values are: \code{"adduct"} (ignore isotope annotation), \code{"isotope"} (ignore adduct annotation), \code{"none"}
#'   (ignore both annotations) and \code{"ignore"} (don't check for charge mismatches). \emph{Important}: when
#'   \command{OpenMS} is used to find features, it already removes any detected non-monoisotopic features by default.
#'   Hence, in such case setting \code{chargeMismatch="adduct"} is more appropriate.
#'
#' @return \code{selectIons} returns a \code{featureGroups} object with only the selected feature groups and amended
#'   with adduct annotations.
#'
#' @aliases selectIons
#' @export
setMethod("selectIons", "featureGroups", function(fGroups, components, prefAdduct, onlyMonoIso = TRUE,
                                                  chargeMismatch = "adduct")
{
    # UNDONE is intensity a proper measure? ie does it allow comparison if
    # isotopes/adducts are taken from different analyses?
    # UNDONE: add logging to see what happens
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(components, "components", add = ac)
    checkmate::assertFlag(onlyMonoIso, add = ac)
    checkmate::assertChoice(chargeMismatch, c("isotope", "adduct", "none", "ignore"))
    checkmate::reportAssertions(ac)
    
    prefAdduct <- as.character(checkAndToAdduct(prefAdduct))
    
    if (is.null(componentInfo(components)[["neutral_mass"]]))
        stop("No adduct/isotope information available in given components!")
    
    cTab <- as.data.table(components)
    cTab <- cTab[group %in% names(fGroups)]
    hasIsos <- !is.null(cTab[["isonr"]]) & !all(is.na(cTab$isonr))
    if (hasIsos)
        cTab <- cTab[!is.na(isonr) | !is.na(adduct_ion)]
    else
        cTab <- cTab[!is.na(adduct_ion)]
    
    cTab[, remove := FALSE]
    
    if (hasIsos && chargeMismatch != "ignore")
    {
        cTab[!is.na(charge) & !is.na(adduct_ion),
             chMismatch := charge != sapply(lapply(adduct_ion, as.adduct), slot, "charge")]
        if (chargeMismatch == "isotope")
            cTab[chMismatch == TRUE, adduct_ion := NA_character_]
        else if (chargeMismatch == "adduct")
            cTab[chMismatch == TRUE, isonr := NA_integer_]
        else # "none"
            cTab[chMismatch == TRUE, remove := TRUE]
    }
    
    if (onlyMonoIso)
    {
        if (!hasIsos)
            cat("No isotope annotations available!\n")
        else
            cTab[remove == FALSE & !is.na(isonr), remove := isonr != 0]
    }
    
    cTab[!is.na(adduct_ion) & remove == FALSE, remove := {
        if (.N == 1)
            FALSE
        # UNDONE: allow >1 pref adducts?
        else if (prefAdduct %in% adduct_ion)
            adduct_ion != prefAdduct
        else # fall back to most intense
            !numEQ(intensity, max(intensity))
    }, by = "name"]
    
    # remove unwanted isotopes/adducts
    fGroups <- fGroups[, setdiff(names(fGroups), cTab[remove == TRUE]$group)]
    
    # annotate remaining
    
    cTabAdd <- cTab[!is.na(adduct_ion) & !remove]

    if (nrow(cTabAdd) == 0)
    {
        fGroups@annotations <- data.table()
        cat("No adduct annotations found!\n")
    }
    else
    {
        fGroups@annotations <- data.table(group = names(fGroups))
        
        if (nrow(cTabAdd) > 0)
        {
            setnames(cTabAdd, "adduct_ion", "adduct")
            fGroups@annotations <- merge(fGroups@annotations, cTabAdd[, c("group", "adduct"), with = FALSE],
                                         by = "group", all.x = TRUE)
            fGroups@annotations[is.na(adduct), adduct := prefAdduct]
        }
        else
            fGroups@annotations[, adduct := prefAdduct]

        fGroups@annotations[, neutralMass := calculateMasses(groupInfo(fGroups)[group, "mzs"],
                                                             lapply(adduct, as.adduct),
                                                             type = "neutral")]
                
        # retain correct order
        fGroups@annotations <- fGroups@annotations[match(names(fGroups), group)]
        
        printf("Removed %d feature groups detected as unwanted adducts/isotopes\n", sum(cTab$remove))
        printf("Annotated %d feature groups with adducts\n", nrow(cTabAdd))
        printf("\tRemaining %d feature groups set as default adduct %s\n", length(fGroups) - nrow(cTabAdd), prefAdduct)
    }
    
    return(fGroups)
})

#' Grouping of features
#'
#' Group equal features across analyses.
#'
#' After \link[=findFeatures]{features have been found}, the next step is to align and group them across analyses. This
#' process is necessary to allow comparison of features between multiple analyses, which otherwise would be difficult
#' due to small deviations in retention and mass data. Thus, algorithms of 'feature groupers' are used to collect
#' features with similar retention and mass data. In addition, advanced retention time alignment algorithms exist to
#' enhance grouping of features even with relative large retention time deviations (\emph{e.g.} possibly observed from
#' analyses collected over a long period). Like \link{findFeatures}, various algorithms are supported which may have
#' many parameters that can be fine-tuned. This fine-tuning is likely to be necessary, since optimal settings often
#' depend on applied methodology and instrumentation.
#'
#' @templateVar func groupFeatures
#' @templateVar what group features
#' @templateVar ex1 groupFeaturesOpenMS
#' @templateVar ex2 groupFeaturesXCMS3
#' @templateVar algos openms,xcms,xcms3,kpic2
#' @templateVar algosSuffix OpenMS,XCMS,XCMS3,KPIC2,SIRIUS
#' @templateVar noParam TRUE
#' @templateVar ret featureGroups
#' @template generic-algo
#' 
#' @param algorithm A \code{character} that specifies the algorithm to be used: either \code{"openms"}, \code{"xcms"},
#'   \code{"xcms3"} or \code{"kpic2"} (\code{features method}), or \code{"sirius"} (\code{data.frame} method).
#' @param obj Either a \code{\link{features}} object to be grouped, or a \code{data.frame} with
#'   \link[=analysis-information]{analysis info} to be passed to \code{groupFeaturesSIRIUS}
#' @param \dots Further parameters passed to the selected grouping algorithm.
#' @param verbose if \code{FALSE} then no text output will be shown.
#'  
#' @return An object of a class which is derived from \code{\link{featureGroups}}.
#' 
#' @templateVar what groupFeatures
#' @templateVar cl features
#' @template main-rd-method
#' @export
setMethod("groupFeatures", "features", function(obj, algorithm, ..., verbose = TRUE)
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

    f(obj, ..., verbose = verbose)
})

#' @details The \code{data.frame} method for \code{groupFeatures} is a special case that currently only supports the
#'   \code{"sirius"} algorithm.
#' @rdname groupFeatures
#' @export
setMethod("groupFeatures", "data.frame", function(obj, algorithm, ..., verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertChoice(algorithm, c("sirius"), add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    f <- switch(algorithm,
                sirius = groupFeaturesSIRIUS)
    
    f(obj, ..., verbose = verbose)
})

#' Import feature groups from files
#'
#' Generic function to import feature groups produced by other software from files.
#'
#' @templateVar func importFeatureGroups
#' @templateVar what import feature groups from files
#' @templateVar ex1 importFeatureGroupsBrukerTASQ
#' @templateVar ex2 importFeatureGroupsBrukerPA
#' @templateVar algosSuffix BrukerPA,BrukerTASQ,EnviMass
#' @templateVar ret featureGroups
#' @templateVar noParam TRUE
#' @template generic-algo
#'
#' @param path The path that should be used for importing. See the algorithm specific functions for more details.
#' @param type Which file type should be imported: \code{"brukerpa"} (Bruker ProfileAnalysis), \code{"brukertasq"}
#'   (Bruker TASQ) or \code{"envimass"} (\pkg{enviMass}).
#' @param \dots Further arguments passed to the selected import algorithm function.
#'
#' @inherit groupFeatures return
#'
#' @seealso \code{\link{groupFeatures}} to group features. Other import functions:
#'   \code{\link{importFeatureGroupsXCMS}}, \code{\link{importFeatureGroupsXCMS3}} and
#'   \code{\link{importFeatureGroupsKPIC2}}.
#'
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
