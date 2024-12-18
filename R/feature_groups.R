# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

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
#'   For \code{normInts}: passed to \code{\link{screenSuspects}} if \code{featNorm="istd"}.
#'
#'   \setsPassedArgs1{featureGroups}
#' @param average If \code{TRUE} then data within replicate groups are averaged.
#'
#'   For \code{as.data.table}: if \code{features=TRUE} other feature properties are also averaged.
#' @param averageFunc Function used for averaging. Only used when \code{average=TRUE} or \code{FCParams != NULL}.
#' @param areas If set to \code{TRUE} then areas are considered instead of peak intensities.
#'
#'   For \code{as.data.table}: ignored if \code{features=TRUE}, as areas of features are always reported.
#' @param normalized If \code{TRUE} then normalized intensity data is used (see the \verb{Feature intensity
#'   normalization} section.
#'
#'   For \code{as.data.table}: if no normalization data is available (\emph{e.g.} because \code{normInts} was not used)
#'   then an automatic group normalization is performed.
#' @param which A character vector with replicate groups used for comparison.
#'
#'   For \code{overlap}: can also be a \code{list} of \code{character} vectors with replicate groups to compare. For
#'   instance, \code{which=list(c("samp1", "samp2"), c("samp3", "samp4"))} returns the overlap between
#'   \code{"samp1"}+\code{"samp2"} and \code{"samp3"}+\code{"samp4"}.
#' @param FCParams A parameter list to calculate Fold change data. See \code{getFCParams} for more details. Set to
#'   \code{NULL} to not perform FC calculations.
#' @param MSLevel Integer vector with the ms levels (i.e., 1 for MS1 and 2 for MS2) to obtain TIC traces.
#' @param retentionRange Range of retention time (in seconds) to collect TIC traces. Should be a numeric vector with
#'   length of two containing the min/max values. Set to NULL to ignore.
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
#' @templateVar delfa a vector of the group intensities, the group name
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
#' @slot ISTDs A \code{data.table} with screening results for internal standards (filled in by the \code{normInts}
#'   method).
#' @slot ISTDAssignments A \code{list}, where each item is named by a feature group and consists of a vector with
#'   feature group names of the internal standards assigned to it (filled in by the \code{normInts} method).
#' @slot concentrations,toxicities A \code{data.table} with predicted concentrations/toxicities for each feature group.
#'   Assigned by the \code{\link{calculateConcs}}/\code{\link{calculateTox}} methods. Use the
#'   \code{concentratrions}/\code{toxicities} methods for access.
#'
#' @author Rick Helmus <\email{r.helmus@@uva.nl}> and Ricardo Cunha <\email{cunha@@iuta.de}> (\code{getTICs} and
#'   \code{getBPCs} functions)
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
                                    groupScores = "data.table", annotations = "data.table",
                                    ISTDs = "data.table", ISTDAssignments = "list", concentrations = "data.table",
                                    toxicities = "data.table"),
                          contains = c("VIRTUAL", "workflowStep"))

setMethod("initialize", "featureGroups", function(.Object, ...)
{
    args <- list(...)

    # data.table's don't seem to initialize well (gives error that slot is init as list)
    for (s in c("groups", "ftindex", "groupQualities", "groupScores", "annotations", "ISTDs", "concentrations",
                "toxicities"))
    {
        if (is.null(args[[s]]))
            args[[s]] <- data.table()
    }
    
    .Object@ISTDAssignments <- makeEmptyListNamed(.Object@ISTDAssignments)

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
    if (length(object) > 0)
    {
        printf("Has normalized intensities: %s\n", as.character(!is.null(featureTable(object)[[1]][["intensity_rel"]])))
        printf("Internal standards used for normalization: ")
        if (nrow(internalStandards(object)) == 0)
            printf("no\n")
        else
            printf("%s (%d assigned total)\n", getStrListWithMax(unique(internalStandards(object)$name), 6, ", "),
                   nrow(internalStandards(object)))
        printf("Predicted concentrations: ")
        concs <- concentrations(object)
        if (nrow(concs) == 0)
            printf("none\n")
        else
            printf("%d feature groups (%.2f%%)\n", uniqueN(concs$group), uniqueN(concs$group) / gCount * 100)
        printf("Predicted toxicities: ")
        tox <- toxicities(object)
        if (nrow(tox) == 0)
            printf("none\n")
        else
            printf("%d feature groups (%.2f%%)\n", uniqueN(tox$group), uniqueN(tox$group) / gCount * 100)
            
    }
    showAnaInfo(analysisInfo(object))
})

#' @describeIn featureGroups Accessor for \code{groups} slot.
#' @aliases groupTable
#' @export
setMethod("groupTable", "featureGroups", function(object, areas = FALSE, normalized = FALSE)
{
    checkmate::assertFlag(areas)
    checkmate::assertFlag(normalized)

    if (length(object) == 0 || (!areas && !normalized))
        return(object@groups)

    anaInfo <- analysisInfo(object)
    ret <- copy(object@groups)
    ftindex <- object@ftindex
    fTable <- featureTable(object)
    
    if (normalized && is.null(fTable[[1]][["intensity_rel"]]))
        stop("There is no normalized data, did you run normInts()?")
    
    colName <- if (areas && normalized)
        "area_rel"
    else if (areas)
        "area"
    else # if ("normalized")
        "intensity_rel"
    
    for (cl in seq_along(ret))
    {
        ftinds <- ftindex[[cl]]
        anainds <- seq_len(nrow(ret))[ftinds != 0]
        ftinds <- ftinds[ftinds != 0]
        as <- mapply(anainds, ftinds, SIMPLIFY = TRUE, FUN = function(a, i)
        {
            fTable[[anaInfo$analysis[a]]][[colName]][i]
        })
        set(ret, anainds, cl, as)
    }
    
    return(ret)
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

#' @describeIn featureGroups Accessor for \code{ISTDs} slot.
#' @aliases internalStandards
#' @export
setMethod("internalStandards", "featureGroups", function(fGroups) fGroups@ISTDs)

#' @describeIn featureGroups Accessor for \code{ISTDAssignments} slot.
#' @aliases internalStandardAssignments
#' @export
setMethod("internalStandardAssignments", "featureGroups", function(fGroups) fGroups@ISTDAssignments)

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

#' @describeIn featureGroups Accessor for \code{concentrations} slot.
#' @aliases concentrations
#' @export
setMethod("concentrations", "featureGroups", function(fGroups) fGroups@concentrations)

#' @describeIn featureGroups Accessor for \code{toxicities} slot.
#' @aliases toxicities
#' @export
setMethod("toxicities", "featureGroups", function(fGroups) fGroups@toxicities)

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
        if (length(obj@concentrations) > 0)
            obj@concentrations <- obj@concentrations[, setdiff(names(obj@concentrations), removedAnas), with = FALSE]
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
        if (nrow(obj@ISTDs) > 0)
        {
            obj@ISTDs <- obj@ISTDs[group %in% names(obj@groups)]
            obj@ISTDAssignments <- internalStandardAssignments(obj)[names(obj@ISTDAssignments) %chin% names(obj@groups)]
        }
        if (nrow(obj@concentrations) > 0)
            obj@concentrations <- obj@concentrations[group %in% names(obj@groups)]
        if (nrow(obj@toxicities) > 0)
            obj@toxicities <- obj@toxicities[group %in% names(obj@groups)]
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
        if (nrow(concentrations(obj)) > 0)
        {
            obj@concentrations <- copy(obj@concentrations)
            anas <- analyses(obj) # update var
            obj@concentrations[group %chin% affectedGrps, (anas) := {
                vals <- mget(anas)
                vals[ftind[[group]] == 0] <- NA_real_
                vals
            }, by = seq_len(sum(obj@concentrations$group %chin% affectedGrps))]
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

setMethod("averageGroups", "featureGroups", function(fGroups, areas, normalized, func)
{
    gTable <- copy(groupTable(fGroups, areas, normalized))
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
#' @param concAggrParams,toxAggrParams Parameters to aggregate calculated concentrations/toxicities (obtained with
#'   \code{\link{calculateConcs}}/\code{\link{calculateTox}}). See \link[=pred-aggr-params]{prediction aggregation
#'   parameters} for more information. Set to \code{NULL} to omit this data.
#' @param normConcToTox Set to \code{TRUE} to normalize concentrations to toxicities. Only relevant if this data is
#'   present (see \code{\link{calculateConcs}}/\code{\link{calculateTox}}).
#' @export
setMethod("as.data.table", "featureGroups", function(x, average = FALSE, areas = FALSE, features = FALSE,
                                                     qualities = FALSE, regression = FALSE, averageFunc = mean,
                                                     normalized = FALSE, FCParams = NULL,
                                                     concAggrParams = getDefPredAggrParams(),
                                                     toxAggrParams = getDefPredAggrParams(), normConcToTox = FALSE)
{
    return(doFGAsDataTable(x, average, areas, features, qualities, regression, averageFunc, normalized, FCParams,
                           concAggrParams, toxAggrParams, normConcToTox))
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
    rGroups <- replicateGroups(fGroups)
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assert(checkmate::checkSubset(which, rGroups, empty.ok = FALSE),
                      checkmate::checkList(which, "character", any.missing = FALSE),
                      .var.name = "which", add = ac)
    checkmate::assertFlag(exclusive, add = ac)
    checkmate::reportAssertions(ac)

    if (length(which) < 2 || length(fGroups) == 0)
        return(fGroups) # nothing to do...

    fGroupsList <- lapply(which, replicateGroupFilter, fGroups = fGroups, verbose = FALSE)
    ov <- Reduce(intersect, lapply(fGroupsList, names))
    ret <- fGroups[, ov]

    if (exclusive)
        ret <- unique(ret, which = unlist(which))
    
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
    checkPackage("MetaClean")
    
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
    EICs <- getEICsForFGroups(obj, EICParams = getDefEICParams(rtWindow = 0))
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
    {
        cat("No adduct/isotope information available in given components!\n")
        return(fGroups)
    }
    if (length(components) == 0)
    {
        cat("Components are empty, skipping...\n")
        return(fGroups)
    }
    
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

#' @describeIn featureGroups Provides various methods to normalizes feature intensities for each sample analysis or of
#'   all features within a feature group. See the \verb{Feature intensity normalization} section below.
#'
#' @param featNorm The method applied for feature normalization: \code{"istd"}, \code{"tic"}, \code{"conc"} or
#'   \code{"none"}. See the \verb{Feature intensity normalization} section for details.
#' @param groupNorm If \code{TRUE} then group normalization is performed. See the \verb{Feature intensity normalization}
#'   section for details.
#' @param normFunc A \code{function} to combine data for normalization. See the \verb{Feature intensity normalization}
#'   section for details.
#' @param standards A \code{data.table} (or \code{data.frame}) with all internal standards. Should follow the format of
#'   a \link[=suspect-screening]{suspect list}. Only used if \code{featNorm="istd"}. See the \verb{Feature intensity
#'   normalization} section for details.
#'
#'   \setsWF Can also be a \code{list} with internal standard lists.
#'
#'   See the \code{suspects} argument to \code{\link{screenSuspects}} for more details.
#' @param ISTDRTWindow,ISTDMZWindow The retention time and \emph{m/z} windows for IS selection. Only used if
#'   \code{featNorm="istd"}. See the \verb{Feature intensity normalization} section for details.
#' @param minISTDs The minimum number of IS that should be assigned to each feature (if possible). Only used if
#'   \code{featNorm="istd"}. See the \verb{Feature intensity normalization} section for details.
#'
#' @section Feature intensity normalization: The \code{normInts} method performs normalization of feature intensities
#'   (and areas). These values are amended in the \code{features} slot, while the original intensities/areas are kept.
#'   To use the normalized intensities set \code{normalized=TRUE} to methods such as \code{\link{plotInt}},
#'   \code{\link{generateComponentsIntClust}} and \code{as.data.table}. Please see the \code{normalized} argument
#'   documentation for these methods for more details.
#'
#'   The \code{normInts} method supports several methods to normalize intensities/areas of features within the same
#'   analysis. Most methods are influenced by the \emph{normalization concentration} (\code{norm_conc} in the
#'   \link[=analysis-information]{analysis information}) set for each sample analysis. For \code{NA} or zero values the
#'   output will be zero. If the \code{norm_conc} is completely absent from the analysis information, the normalization
#'   concentration is defaulted to one.
#'
#'   The different normalization methods are:
#'
#'   \enumerate{
#'
#'   \item \code{featNorm="istd"} Uses \emph{internal standards} (IS) for normalization. The IS are screened internally
#'   by the \code{\link{screenSuspects}} function. Hence, the IS specified by the \code{standards} argument should
#'   follow the format of a \link[=suspect-screening]{suspect list}. Note that labelled elements in IS formulae should
#'   be specified with the \CRANpkg{rcdk} format, \emph{e.g.} \code{"[13]C"} for 13C, \code{"[2]H"} for a deuterium etc.
#'   Example IS lists are provided with the \pkg{patRoonData} package.
#'
#'   The assignment of IS to features is automatically performed, using the following criteria: \enumerate{
#'
#'   \item Only analyses are considered with a defined normalization concentration.
#'
#'   \item The IS must be detected in all of the analyses in which the feature was detected.
#'
#'   \item The retention time and \emph{m/z} are reasonably close (\code{ISTDRTWindow}/\code{ISTDMZWindow} arguments).
#'   However, additional IS candidates outside these windows will be chosen if the number of candidates is less than the
#'   \code{minISTDs} argument. In this case the next close(st) candidate(s) will be chosen.
#'
#'   }
#'
#'   Normalization of features within the same feature group always occur with the same IS. If multiple IS are assigned
#'   to a feature then normalization occurs with the combined intensity (area), which is calculated with the function
#'   defined by the \code{normFunc} argument. The (combined) IS intensity is then normalized by the normalization
#'   concentration, and finally used for feature normalization.
#'
#'   \item \code{featNorm="tic"} Uses the Total Ion Current (TIC) to normalize intensities. The TIC is calculated by
#'   combining all intensities with the function defined by the \code{normFunc} argument. For this reason, you may need
#'   to take care to perform normalization before \emph{e.g.} suspect screening or other prioritization techniques. The
#'   TIC normalized intensities are finally divided by the normalization concentration.
#'
#'   \item \code{featNorm="conc"} Simply divides all intensities (areas) with the normalization concentration defined
#'   for the sample.
#'
#'   \item \code{featNorm="none"} Performs no normalization. The raw intensity values are simply copied. This is mainly
#'   useful if you only want to do group normalization (described below).
#'
#'   }
#'
#'   The meaning of the normalization concentration differs for each method: for \code{"istd"} it resembles the IS
#'   concentration of a sample analysis, whereas for \code{"tic"} and \code{"conc"} it is used to normalize different
#'   sample amounts (\emph{e.g.} injection volume).
#'
#'   If \code{groupNorm=TRUE} then feature intensities (areas) will be normalized by the combined values for its feature
#'   group (again, combination occurs with \code{normFunc}). This \emph{group normalization} always occurs \emph{after}
#'   aforementioned normalization methods. Group normalization was the only method with \pkg{patRoon} \samp{<2.1}, and
#'   still occurs automatically if \code{normInts} was not called when a method is executed that requests normalized
#'   data.
#'
#' @return \code{normInts} returns a \code{featureGroups} object, amended with data in the \code{ISTDs} and
#'   \code{ISTDAssignments} slots if \code{featNorm="istd"}.
#'
#' @aliases normInts
#' @export
setMethod("normInts", "featureGroups", function(fGroups, featNorm, groupNorm, normFunc, standards, ISTDRTWindow,
                                                ISTDMZWindow, minISTDs, ...)
{
    # NOTE: keep args in sync with sets method
    
    # UNDONE: default for minISTDs OK? (or no default for minISTDs and ISTDRTWindow?)
    # UNDONE: ISTD: doc that sorting occurs on both RT and m/z deviation
    # UNDONE: add adduct argument here, to make it more clear that it needs to be specified if featNorm=="istd"?
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertSubset(featNorm, c("tic", "istd", "conc", "none"))
    checkmate::assertFlag(groupNorm, add = ac)
    checkmate::assertFunction(normFunc, add = ac)
    checkmate::assertDataFrame(standards, null.ok = featNorm != "istd", add = ac) # more asserts in screenSuspects()
    aapply(checkmate::assertNumber, . ~ ISTDRTWindow + ISTDMZWindow, lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::assertCount(minISTDs, positive = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(fGroups) == 0)
        return(fGroups)
    
    anaInfo <- analysisInfo(fGroups)
    
    normConcs <- anaInfo[["norm_conc"]]
    if (is.null(normConcs))
    {
        printf("NOTE: No normalization concentrations defined (norm_conc column is absent in analysis information). Defaulting to 1.\n")
        normConcs <- rep(1, nrow(anaInfo))
    }
    
    # reset
    fGroups@ISTDs <- data.table()
    fGroups@ISTDAssignments <- list()
    
    if (featNorm == "istd")
    {
        # HACK: what we should do here for screening is exactly the same as screenSuspects(). So simply call that and use
        # its output...
        fGroupsScr <- screenSuspects(fGroups, suspects = standards, ...)
        fGroups@ISTDs <- screenInfo(fGroupsScr)
        origN <- uniqueN(fGroups@ISTDs$name)
        
        # only keep hits that are present in the analyses with non-NA conc
        fGroupsWithISTD <- fGroups[, fGroups@ISTDs$group]
        fGroupsWithISTD <- fGroupsWithISTD[!is.na(normConcs) & normConcs != 0]
        fGroupsWithISTD <- minAnalysesFilter(fGroupsWithISTD, relThreshold = 1, verbose = FALSE)
        
        fGroups@ISTDs <- fGroups@ISTDs[group %in% names(fGroupsWithISTD)]
        
        printf("Removed %d non-ubiquitous internal standards\n", origN - uniqueN(fGroups@ISTDs$name))
        
        gInfo <- groupInfo(fGroups)
        gInfoISTDs <- gInfo[unique(fGroups@ISTDs$group), ]
        fGroups@ISTDAssignments <- setNames(Map(gInfo$rts, gInfo$mzs, f = function(rt, mz)
        {
            # UNDONE: with configurable N, handle duplicate IS assignments differently? (ie select on suspect RT instead of group RT)
            
            # sort by closest eluting ISTDs with closest m/z
            gi <- gInfoISTDs[order(abs(gInfoISTDs$rts - rt), abs(gInfoISTDs$mzs - mz)), ]
            giInRange <- gi[numLTE(abs(gi$rts - rt), ISTDRTWindow) & numLTE(abs(gi$mzs - mz), ISTDMZWindow), ]
            if (nrow(giInRange) >= minISTDs)
                return(rownames(giInRange)) # only take those in range
            if (nrow(gi) > minISTDs)
                gi <- gi[seq_len(minISTDs), ] # just take minimum of ISTDs, even if some are out of range
            return(rownames(gi))
        }), names(fGroups))
        fGroups@ISTDAssignments <- fGroups@ISTDAssignments[lengths(fGroups@ISTDAssignments) > 0]
        fGroups@ISTDAssignments <- fGroups@ISTDAssignments[!names(fGroups@ISTDAssignments) %in% fGroups@ISTDs$group]
        
        fGroups@features@features <- Map(featureTable(fGroups), normConcs, f = function(ft, nconc)
        {
            ft <- copy(ft)
            
            if (is.na(nconc) || nconc == 0)
                ft[, c("intensity_rel", "area_rel") := 0]
            else
            {
                ft[group %in% names(fGroups@ISTDAssignments), intensity_rel := mapply(intensity, group, FUN = function(int, grp)
                {
                    iint <- normFunc(ft[group %in% fGroups@ISTDAssignments[[grp]]]$intensity)
                    return(int / (iint / nconc))
                })]
                ft[group %in% names(fGroups@ISTDAssignments), area_rel := mapply(area, group, FUN = function(ar, grp)
                {
                    iar <- normFunc(ft[group %in% fGroups@ISTDAssignments[[grp]]]$area)
                    return(ar / (iar / nconc))
                })]
                # HACK: don't want NA values
                ft[is.na(intensity_rel), c("intensity_rel", "area_rel") := 0]
            }
            return(ft)
        })
    }
    else if (featNorm == "tic")
    {
        fGroups@features@features <- Map(featureTable(fGroups), normConcs, f = function(ft, nconc)
        {
            ft <- copy(ft)
            nint <- normFunc(ft$intensity) * nconc
            narea <- normFunc(ft$area) * nconc
            if (is.na(nconc) || nconc == 0 || nint == 0)
                ft[, c("intensity_rel", "area_rel") := 0]
            else
                ft[, c("intensity_rel", "area_rel") := .(intensity / nint, area / narea)]
            return(ft)
        })
    }
    else if (featNorm == "conc")
    {
        fGroups@features@features <- Map(featureTable(fGroups), normConcs, f = function(ft, nconc)
        {
            ft <- copy(ft)
            if (is.na(nconc) || nconc == 0)
                ft[, c("intensity_rel", "area_rel") := 0]
            else
                ft[, c("intensity_rel", "area_rel") := .(intensity / nconc, area / nconc)]
            return(ft)
        })
    }
    else # "none"
    {
        fGroups@features@features <- lapply(featureTable(fGroups), function(ft)
        {
            ft <- copy(ft)
            ft[, c("intensity_rel", "area_rel") := .(intensity, area)]
            return(ft)
        })
    }
    
    if (groupNorm)
    {
        gNames <- names(fGroups)
        
        featsPerGroup <- split(rbindlist(featureTable(fGroups)), by = "group")
        nInts <- sapply(lapply(featsPerGroup, "[[", "intensity_rel"), normFunc)
        nAreas <- sapply(lapply(featsPerGroup, "[[", "area_rel"), normFunc)
        
        fGroups@features@features <- lapply(featureTable(fGroups), function(ft)
        {
            ft <- copy(ft)
            ft[, c("intensity_rel", "area_rel") := .(intensity_rel / nInts[group],
                                                     area_rel / nAreas[group])]
            return(ft)
        })
    }
    
    return(fGroups)
})

#' @aliases calculateConcs
#' @rdname pred-quant
#' @export
setMethod("calculateConcs", "featureGroups", function(fGroups, featureAnn, areas = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(featureAnn, "featureAnnotations", add = ac)
    checkmate::assertFlag(areas, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(fGroups) == 0)
    {
        cat("No feature groups, nothing to do...\n")
        return(fGroups)
    }
    
    featureAnn <- featureAnn[names(fGroups)]
    
    if (length(featureAnn) == 0)
    {
        cat("No (relevant) feature annotations, nothing to do...\n")
        fGroups@concentrations <- data.table()
        return(fGroups)
    }
    
    annTab <- as.data.table(featureAnn)
    if (is.null(annTab[["RF_SMILES"]]) && is.null(annTab[["RF_SIRFP"]]))
        stop("Feature annotations lack predicted response factors. Please call predictRespFactors() first!",
             call. = FALSE)
    
    concs <- data.table()
    
    if (!is.null(annTab[["RF_SMILES"]]) && any(!is.na(annTab$RF_SMILES)))
    {
        resp <- annTab[!is.na(RF_SMILES), c("group", "SMILES", "RF_SMILES"), with = FALSE]
        resp[, type := "compound"]
        setnames(resp, c("SMILES", "RF_SMILES"), c("candidate", "RF"))
        if (!is.null(annTab[["compoundName"]]))
            resp[, candidate_name := annTab$compoundName]
        concs <- calcFeatureConcs(fGroups, resp, areas)
    }
    
    if (!is.null(annTab[["RF_SIRFP"]]) && any(!is.na(annTab$RF_SIRFP)))
    {
        resp <- annTab[!is.na(RF_SIRFP), c("group", "neutral_formula", "RF_SIRFP"), with = FALSE]
        resp[, type := "SIRIUS_FP"]
        setnames(resp, c("neutral_formula", "RF_SIRFP"), c("candidate", "RF"))
        # UNDONE: collapse compoundName? Results mostly in very long columns.
        # if (!is.null(annTab[["compoundName"]]))
        #     resp[, candidate_name := annTab$compoundName]
        resp <- unique(resp, by = c("group", "candidate"))
        concs <- rbind(concs, calcFeatureConcs(fGroups, resp, areas), fill = TRUE)
    }

    if (nrow(concs) > 0)
        concs <- finalizeFeaturePredTab(concs)
    
    fGroups@concentrations <- concs
    
    return(fGroups)
})

#' @aliases calculateTox
#' @rdname pred-tox
#' @export
setMethod("calculateTox", "featureGroups", function(fGroups, featureAnn)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(featureAnn, "featureAnnotations", add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(fGroups) == 0)
    {
        cat("No feature groups, nothing to do...\n")
        return(fGroups)
    }
    
    featureAnn <- featureAnn[names(fGroups)]
    
    if (length(featureAnn) == 0)
    {
        cat("No (relevant) feature annotations, nothing to do...\n")
        fGroups@toxicities <- data.table()
        return(fGroups)
    }
    
    annTab <- as.data.table(featureAnn)
    if (is.null(annTab[["LC50_SMILES"]]) && is.null(annTab[["LC50_SIRFP"]]))
        stop("Feature annotations lack predicted LC50 values. Please call predictTox() first!", call. = FALSE)
    
    toxicities <- data.table()
    
    if (!is.null(annTab[["LC50_SMILES"]]) && any(!is.na(annTab$LC50_SMILES)))
    {
        LC50Tab <- annTab[!is.na(LC50_SMILES), c("group", "SMILES", "LC50_SMILES"), with = FALSE]
        LC50Tab[, type := "compound"]
        setnames(LC50Tab, c("SMILES", "LC50_SMILES"), c("candidate", "LC50"))
        if (!is.null(annTab[["compoundName"]]))
            LC50Tab[, candidate_name := annTab$compoundName]
        toxicities <- LC50Tab
    }
    
    if (!is.null(annTab[["LC50_SIRFP"]]) && any(!is.na(annTab$LC50_SIRFP)))
    {
        LC50Tab <- annTab[!is.na(LC50_SIRFP), c("group", "neutral_formula", "LC50_SIRFP"), with = FALSE]
        LC50Tab[, type := "SIRIUS_FP"]
        setnames(LC50Tab, c("neutral_formula", "LC50_SIRFP"), c("candidate", "LC50"))
        # UNDONE: collapse compoundName? Results mostly in very long columns.
        # if (!is.null(annTab[["compoundName"]]))
        #     LC50Tab[, candidate_name := annTab$compoundName]
        LC50Tab <- unique(LC50Tab, by = c("group", "candidate"))
        toxicities <- rbind(toxicities, LC50Tab, fill = TRUE)
    }
    
    if (nrow(toxicities) > 0)
        toxicities <- finalizeFeaturePredTab(toxicities)
    
    fGroups@toxicities <- toxicities
    
    return(fGroups)
})

#' @describeIn featureGroups Obtain the total ion chromatogram/s (TICs) of the analyses.
#' @export
setMethod("getTICs", "featureGroups", function(obj, retentionRange = NULL, MSLevel = c(1, 2))
{
    getTICs(obj@features, retentionRange, MSLevel)
})

#' @describeIn featureGroups Obtain the base peak chromatogram/s (BPCs) of the analyses.
#' @export
setMethod("getBPCs", "featureGroups", function(obj, retentionRange = NULL, MSLevel = c(1, 2))
{
    getBPCs(obj@features, retentionRange, MSLevel)
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
