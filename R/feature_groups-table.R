#' @include main.R
NULL

getADTIntCols <- function(wh) paste0(wh, "_intensity")
stripADTIntSuffix <- function(cols) sub("_intensity$", "", cols)
getADTRegCols <- function() c("RSQ", "intercept", "slope", "p")

mergeScreenInfoWithDT <- function(tab, scrInfo, collapseSuspects, onlyHits)
{
    scrInfo <- copy(scrInfo)
    setnames(scrInfo, names(scrInfo), paste0("susp_", names(scrInfo))) # add susp_ column prefixes
    
    if (!is.null(collapseSuspects))
    {
        scrInfo[, susp_name := paste0(susp_name, collapse = collapseSuspects), by = "susp_group"]
        # only keep unique and remove suspect specific columns
        scrInfo <- unique(scrInfo[, c("susp_group", "susp_name"), with = FALSE], by = "susp_group")
    }
    
    ret <- merge(tab, scrInfo, by.x = "group", by.y = "susp_group", all.x = !onlyHits, sort = FALSE,
                 allow.cartesian = is.null(collapseSuspects))
    return(ret)
}

doFGAADTGroups <- function(fGroups, intColNames, average, averageBy, areas, addQualities, addScores, regression,
                           regressionBy, averageFunc, normalized, FCParams, concAggrParams, toxAggrParams,
                           normConcToTox, collapseSuspects, onlyHits)
{
    anaInfo <- analysisInfo(fGroups)
    gNames <- names(fGroups)
    gInfo <- groupInfo(fGroups)
    
    if (normalized)
        fGroups <- maybeAutoNormalizeFGroups(fGroups)
    
    gTable <- if (isFALSE(average))
        groupTable(fGroups, areas, normalized)
    else
        averageGroups(fGroups, areas, normalized, by = averageBy, func = averageFunc)
    
    ret <- transpose(gTable)
    setnames(ret, intColNames)
    if (!isFALSE(regression) && length(intColNames) > 1)
    {
        if (is.null(regressionBy))
        {
            xvec <- anaInfo[, mean(get(regression)), by = averageBy][[2]]
            regr <- sapply(gTable, calcFeatureRegression, xvec = xvec, simplify = FALSE)
            regr <- lapply(regr, "[", getADTRegCols())
            ret <- cbind(ret, rbindlist(regr))
        }
        else
        {
            for (rb in unique(anaInfo[[regressionBy]]))
            {
                rbAnaInfo <- anaInfo[get(regressionBy) == rb]
                rbxvec <- rbAnaInfo[, mean(get(regression)), by = averageBy][[2]] # average x values if needed
                rbIntCols <- getADTIntCols(unique(rbAnaInfo[[averageBy]]))
                rbRows <- match(rbIntCols, intColNames)
                
                regr <- sapply(gTable[rbRows], function(ints) calcFeatureRegression(rbxvec, ints)[getADTRegCols()],
                               simplify = FALSE)
                regr <- rbindlist(regr)
                setnames(regr, paste0(names(regr), "_", rb))
                ret <- cbind(ret, regr)
            }
        }
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
        
        gTableNonAvg <- groupTable(fGroups, areas, normalized)
        gTableAvg <- averageGroups(fGroups, areas, normalized, by = "group", func = averageFunc) # UNDONE: support averageBy
        
        repInds <- match(FCParams$rGroups, replicateGroups(fGroups))
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
        ret <- cbind(ret, groupQualities(fGroups)[match(ret$group, group), -"group"])
    if (addScores)
        ret <- cbind(ret, groupScores(fGroups)[match(ret$group, group), -"group"])
    
    annTable <- annotations(fGroups)
    if (nrow(ret) > 0 && nrow(annTable) > 0)
    {
        if (isFGSet(fGroups))
        {
            # split adduct/ion_mz into set specific columns
            annTable <- dcast(annTable, group + neutralMass ~ set, value.var = c("adduct", "ion_mz"), sep = "-",
                              drop = c(TRUE, FALSE))
            
            # in case a set doesn't have any features then it will have no columns --> add them manually
            missingSets <- setdiff(sets(fGroups), annotations(fGroups)$set)
            for (ms in missingSets)
                set(annTable, j = c(paste0("adduct-", ms), paste0("ion_mz-", ms)), value = list(NA_character_, NA_real_))
            
            # fixup order
            setcolorder(annTable, c("group", paste0("ion_mz-", sets(fGroups)), paste0("adduct-", sets(fGroups))))
        }
        ret <- merge(ret, annTable, by = "group", sort = FALSE)
    }
    
    ISTDAssign <- internalStandardAssignments(fGroups)
    if (length(ISTDAssign) > 0 && nrow(ret) > 0)
    {
        colISTDs <- function(ia) paste0(ia, collapse = ",")
        
        if (isFGSet(fGroups))
        {
            if (!is.null(ret[["set"]]))
                ret[, ISTD_assigned := sapply(ISTDAssign[[set[1]]][group], colISTDs), by = "set"]
            else
            {
                for (s in sets(fGroups))
                    set(ret, j = paste0("ISTD_assigned-", s), value = sapply(ISTDAssign[[s]][ret$group], colISTDs))
            }
        }
        else
            ret[, ISTD_assigned := sapply(ISTDAssign[group], function(ia) paste0(ia, collapse = ","))]
    }
    
    # NOTE: do this _before_ adding concs/tox
    if (isScreening(fGroups))
        ret <- mergeScreenInfoWithDT(ret, screenInfo(fGroups), collapseSuspects, onlyHits)
    
    splitSusps <- isScreening(fGroups) && is.null(collapseSuspects)
    mergePred <- function(tab, pred, mby, typesCol)
    {
        if (splitSusps && !is.null(pred[["candidate_name"]]))
        {
            # UNDONE: until incomparables arg becomes available for data.table::merge(), we use a dummy column so merges
            # between suspects/suspects and non-suspects/non-suspects can be done in one go.
            tab <- copy(tab); pred <- copy(pred)
            pred[, susp_merge := fifelse(get(typesCol) == "suspect", candidate_name, "nosusp")]
            tab[, susp_merge := fifelse(!is.na(susp_name) & susp_name %in% pred[get(typesCol) == "suspect"]$candidate_name,
                                        susp_name, "nosusp")]
            tab <- merge(tab, pred, by = c(mby, "susp_merge"), all.x = TRUE, sort = FALSE)
            tab[, susp_merge := NULL]
        }
        else
            tab <- merge(tab, pred, by = mby, all.x = TRUE, sort = FALSE)
        return(tab)
    }
    
    if (!is.null(concAggrParams) && nrow(concentrations(fGroups)) > 0)
    {
        # merge concentrations:
        #   - if suspects are collapsed or concs are not from suspect data, then merging should be done by fGroup
        #   - else concs from suspect data should be merged by group/suspect in ret
        
        concs <- aggregateConcs(concentrations(fGroups), anaInfo, concAggrParams, splitSusps)
        setnames(concs, "type", "conc_types")
        
        if (!isFALSE(average))
        {
            for (icol in stripADTIntSuffix(intColNames))
            {
                # NOTE: intColNames will just be "intensity" if averaging by fGroups
                anas <- if (averageBy == "fGroups") anaInfo$analysis else anaInfo[get(averageBy) == icol]$analysis
                concs[, (paste0(icol, "_conc")) := aggrVec(unlist(.SD), averageFunc), .SDcols = anas, by = seq_len(nrow(concs))]
            }
            concs[, (anaInfo$analysis) := NULL]
        }
        else
            setnames(concs, anaInfo$analysis, paste0(anaInfo$analysis, "_conc"))
        
        setcolorder(concs, setdiff(names(concs), "conc_types")) # move to end
        
        ret <- mergePred(ret, concs, "group", "conc_types")
        ret <- removeDTColumnsIfPresent(ret, "candidate_name")
    }
    
    if (!is.null(toxAggrParams) && nrow(toxicities(fGroups)) > 0)
    {
        # merge like concs, but a bit simpler since toxicities are assigned per fGroup instead of feature
        
        tox <- aggregateTox(toxicities(fGroups), toxAggrParams, splitSusps)
        setnames(tox, "type", "LC50_types")
        setcolorder(tox, setdiff(names(tox), "LC50_types")) # move to end
        
        ret <- mergePred(ret, tox, "group", "LC50_types")
        ret <- removeDTColumnsIfPresent(ret, "candidate_name")
        
        if (normConcToTox && !is.null(ret[["conc_types"]])) # conc_types should be present if concentration data is
        {
            cols <- grep("_conc$", names(ret), value = TRUE)
            ret[, (cols) := lapply(.SD, "/", LC50), .SDcols = cols]
        }
    }

    return(ret[])
}

doFGAADTFeatures <- function(fGroups, fgTab, intColNames, average, averageBy, addQualities, addScores, regression,
                             regressionBy, averageFunc, anaInfoCols)
{
    anaInfo <- analysisInfo(fGroups)
    
    # prepare feature properties to merge
        
    featTab <- as.data.table(getFeatures(fGroups))
    
    # remove adduct column from sets data: the fGroup adduct is already present, and adducts are assigned per
    # group/set anyway
    featTab <- removeDTColumnsIfPresent(featTab, "adduct")
    
    # if feature qualities/scores are present, then they are already available in featTab. Hence
    # 1. remove them if they should _not_ be reported
    # 2. otherwise replace the feature score specific properties from the group table, as these are feature averaged
    
    if (!addQualities)
        featTab <- removeDTColumnsIfPresent(featTab, featureQualityNames(group = FALSE))
    else
        fgTab <- removeDTColumnsIfPresent(fgTab, featureQualityNames(group = FALSE))
    if (!addScores)
        featTab <- removeDTColumnsIfPresent(featTab, featureQualityNames(group = FALSE, scores = TRUE))
    else
        fgTab <- removeDTColumnsIfPresent(fgTab, featureQualityNames(group = FALSE, scores = TRUE))
    
    by <- "group"
    if (averageBy != "fGroups")
        by <- c(by, "average_group")
    
    if (isFALSE(average))
    {
        featTab[, replicate_group := anaInfo$group[match(analysis, anaInfo$analysis)]]
        featTab[, average_group := analysis]
    }
    else
    {
        if (averageBy != "fGroups")
            featTab[, average_group := anaInfo[[averageBy]][match(analysis, anaInfo$analysis)]]
        
        featTab <- removeDTColumnsIfPresent(featTab, c("isocount", "analysis", "ID", "set"))
        numCols <- setdiff(names(which(sapply(featTab, is.numeric))), "average_group")
        featTab[, (numCols) := lapply(.SD, averageFunc), .SDcols = numCols, by = by]
        featTab <- unique(featTab, by = by)
    }
    
    # prepare main table for merge
    
    if (averageBy != "fGroups")
    {
        # Melt by intensity column to get the proper format. Afterwards, we remove the dummy intensity column, as we
        # want the raw feature intensity data.
        mCols <- list(intensity = intColNames)
        
        concCols <- paste0(stripADTIntSuffix(intColNames), "_conc")
        if (all(concCols %in% names(fgTab)))
            mCols <- c(mCols, list(conc = concCols))
        
        fgTab <- melt(fgTab, measure.vars = mCols, variable.name = "average_group", variable.factor = FALSE,
                      value.name = "intensity")
        fgTab <- fgTab[intensity != 0]
        
        if (length(mCols) > 1)
        {
            # DT changes the analyses names to (character) indices with >1 measure vars: https://github.com/Rdatatable/data.table/issues/4047
            fgTab[, average_group := intColNames[as.integer(average_group)]]
        }
        
        fgTab[, average_group := stripADTIntSuffix(average_group)]
    }
    fgTab <- removeDTColumnsIfPresent(fgTab, "intensity")
    setnames(fgTab, c("ret", "mz"), c("group_ret", "group_mz"))
    annCols <- grep("^(neutralMass|ion_mz|adduct)", names(fgTab), value = TRUE)
    if (length(annCols) > 0)
        setnames(fgTab, annCols, paste0("group_", annCols))
    
    # merge
    fgTab <- merge(fgTab, featTab, by = by, sort = FALSE)
    
    # fixup merged table
    
    if (!isFALSE(regression) && length(intColNames) > 0)
    {
        if (!is.null(regressionBy))
        {
            # combine split regression columns
            fgTab[, (c(getADTRegCols(), "regression_group")) := {
                # get corresponding regressionBy value from average group
                rb <- anaInfo[get(averageBy) == average_group][[regressionBy]][1]
                c(mget(paste0(getADTRegCols(), "_", rb)), rb)
            }, by = "average_group"]
            # remove specific regression columns
            regByCols <- grep(sprintf("^(%s)_(%s)$", paste0(getADTRegCols(), collapse = "|"),
                                      paste0(unique(anaInfo[[regressionBy]]), collapse = "|")),
                              names(fgTab), value = TRUE)
            fgTab[, (regByCols) := NULL]
        }
        fgTab[, x_reg := (intensity - intercept) / slope] # y = ax+b
    }
    
    # set nice column order
    qualCols <- c(featureQualityNames(), featureQualityNames(scores = TRUE))
    colord <- c("group", "set", "analysis", "average_group", "replicate_group", "group_ret", "group_mz",
                "ID", "ret", "mz", "ion_mz", "intensity", "area", "intensity_rel", "area_rel")
    colord <- c(colord, setdiff(names(featTab), c(colord, qualCols)))
    colord <- c(colord, grep("group_(ion_mz|adduct)", names(fgTab), value = TRUE), "neutralMass", "x_reg",
                getADTRegCols(), "regression_group", featureQualityNames(),
                featureQualityNames(scores = TRUE))
    setcolorder(fgTab, intersect(colord, names(fgTab)))
    
    # restore order
    fgTab[, gorder := match(group, names(fGroups))]
    if (averageBy != "fGroups")
    {
        fgTab[, aorder := match(average_group, unique(anaInfo[[averageBy]]))]
        setorderv(fgTab, c("gorder", "aorder"))
    }
    else
        setorderv(fgTab, "gorder")
    fgTab <- removeDTColumnsIfPresent(fgTab, c("gorder", "aorder"))
    
    if (length(anaInfoCols) > 0)
    {
        ai <- anaInfo[, c(averageBy, anaInfoCols), with = FALSE]
        if (!isFALSE(average))
            ai[, (anaInfoCols) := lapply(.SD, averageFunc), .SDcols = anaInfoCols, by = averageBy]
        ai <- unique(ai, by = averageBy)
        ai <- ai[match(fgTab$average_group, get(averageBy))][, (averageBy) := NULL]
        anaInfoCols <- paste0("anaInfo_", anaInfoCols)
        fgTab[, (anaInfoCols) := ai]
    }
    
    if (averageBy == "analysis") # no averaging
        fgTab[, average_group := NULL] # no need for this

    return(fgTab[])
}


# this combines all functionality from all fGroup as.data.table methods, a not so pretty but pragmatic solution...
doFGAsDataTable <- function(fGroups, average = FALSE, areas = FALSE, features = FALSE, qualities = FALSE,
                            regression = FALSE, regressionBy = NULL, averageFunc = mean, normalized = FALSE,
                            FCParams = NULL, concAggrParams = getDefPredAggrParams(),
                            toxAggrParams = getDefPredAggrParams(), normConcToTox = FALSE, anaInfoCols = NULL,
                            collapseSuspects = ",", onlyHits = FALSE)
{
    anaInfo <- analysisInfo(fGroups)
    if (isTRUE(regression))
        regression <- "conc" # legacy
    
    assertFGAsDataTableArgs(fGroups, areas, features, qualities, regression, regressionBy, averageFunc, normalized,
                            FCParams, concAggrParams, toxAggrParams, normConcToTox, anaInfoCols, collapseSuspects,
                            onlyHits)
    averageBy <- assertAndPrepareAnaInfoBy(average, anaInfo, TRUE)
    
    if (length(fGroups) == 0)
        return(data.table(mz = numeric(), ret = numeric(), group = character()))

    intColNames <- if (averageBy == "fGroups") "intensity" else getADTIntCols(unique(anaInfo[[averageBy]]))    
    if (!isFALSE(regression) && (sum(!is.na(anaInfo[[regression]]) < 2) || averageBy == "fGroups"))
        regression <- FALSE
    addQualities <- !isFALSE(qualities) && qualities %in% c("both", "quality") && hasFGroupScores(fGroups)
    addScores <- !isFALSE(qualities) && qualities %in% c("both", "score") && hasFGroupScores(fGroups)
    
    if (!isFALSE(regression))
    {
        if (averageBy == "fGroups")
            stop("Cannot perform regression if averageBy=\"fGroups\"", call. = FALSE)
        if (!is.null(regressionBy) && !isFALSE(average))
            checkAnaInfoAggrGrouping(anaInfo, "averaged", averageBy, regressionBy)
    }
    
    if (length(anaInfoCols) > 0)
    {
        if (!features)
            warning("The anaInfoCols argument is only supported if features = TRUE", call. = FALSE)
        else if (!isFALSE(average))
        {
            notNum <- which(!sapply(anaInfo[, anaInfoCols, with = FALSE], is.numeric))
            if (length(notNum) > 0)
            {
                stop("Cannot average because the following columns from anaInfoCols are not numeric: ",
                     paste0(names(notNum), collapse = ", "), call. = FALSE)
            }
        }
    }

    ret <- doFGAADTGroups(fGroups, intColNames, average, averageBy, areas, addQualities, addScores, regression,
                          regressionBy, averageFunc, normalized, FCParams, concAggrParams, toxAggrParams,
                          normConcToTox, collapseSuspects, onlyHits)    

    if (features)
    {
        ret <- doFGAADTFeatures(fGroups, ret, intColNames, average, averageBy, addQualities, addScores, regression,
                                regressionBy, averageFunc, anaInfoCols)
    }
    
    return(ret[])
}

# both non-sets and sets
doFGScrAsDataTable <- function(x, ..., collapseSuspects = ",", onlyHits = FALSE)
{
    return(doFGAsDataTable(x, ..., collapseSuspects = collapseSuspects, onlyHits = onlyHits))
}


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
                                                     qualities = FALSE, regression = FALSE, regressionBy = NULL,
                                                     averageFunc = mean, normalized = FALSE, FCParams = NULL,
                                                     concAggrParams = getDefPredAggrParams(),
                                                     toxAggrParams = getDefPredAggrParams(), normConcToTox = FALSE,
                                                     anaInfoCols = NULL)
{
    return(doFGAsDataTable(x, average, areas, features, qualities, regression, regressionBy, averageFunc, normalized,
                           FCParams, concAggrParams, toxAggrParams, normConcToTox, anaInfoCols))
})

#' @describeIn featureGroupsScreening Obtain a summary table (a \code{\link{data.table}}) with retention, \emph{m/z},
#'   intensity and optionally other feature data. Furthermore, the output table will be merged with information from
#'   \code{screenInfo}, such as suspect names and other properties and annotation data.
#'
#' @param collapseSuspects If a \code{character} then any suspects that were matched to the same feature group are
#'   collapsed to a single row and suspect names are separated by the value of \code{collapseSuspects}. If \code{NULL}
#'   then no collapsing occurs, and each suspect match is reported on a single row. See the \verb{Suspect collapsing}
#'   section below for additional details.
#'
#' @section {Suspect collapsing}: The \code{as.data.table} method fir \code{featureGroupsScreening} supports an
#'   additional format where each suspect hit is reported on a separate row (enabled by setting
#'   \code{collapseSuspects=NULL}). In this format the suspect
#'   properties from the \code{screenInfo} method are merged with each suspect row. Alternatively, if \emph{suspect
#'   collapsing} is enabled (the default) then the regular \code{as.data.table} format is used, and amended with the
#'   names of all suspects matched to a feature group (separated by the value of the \code{collapseSuspects} argument).
#'
#'   Suspect collapsing also influences how calculated feature concentrations/toxicities are reported (\emph{i.e.}
#'   obtained with \code{\link{calculateConcs}}/\code{\link{calculateTox}}). If these values were directly predicted for
#'   suspects, \emph{i.e.} by using \code{\link{predictRespFactors}}/\code{\link{predictTox}} on the feature groups
#'   object, \emph{and} suspects are \emph{not} collapsed, then the calculated concentration/toxicity reported for each
#'   suspect row is not aggregated and specific for that suspect (unless not available). Hence, this allows you to
#'   obtain specific concentration/toxicity values for each suspect/feature group pair.
#'
#' @export
setMethod("as.data.table", "featureGroupsScreening", doFGScrAsDataTable)

#' @rdname featureGroupsScreening-class
#' @export
setMethod("as.data.table", "featureGroupsScreeningSet", doFGScrAsDataTable)
