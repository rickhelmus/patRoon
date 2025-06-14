#' @include main.R
NULL

#' Feature group to table conversion
#'
#' Convert feature group data to a \code{\link{data.table}} (or \code{data.frame}).
#'
#' The \code{as.data.table} generic function converts most feature group data to a highly customizable
#' \code{\link{data.table}}. If a classical \code{data.frame} is preferred, the \code{as.data.frame} generic function
#' can be used instead and accepts the exact same arguments. The methods defined for \link[=suspect-screening]{suspect
#' screening worfklows} will merge the information from \code{screenInfo}, such as
#' suspect names and other properties and annotation data.
#'
#' @param x The \code{featureGroups} like object to be exported as table.
#' @param average Controls the averaging of feature intensities. Averaging also influences the calculation of regression
#'   parameters. Set to \code{FALSE} to disable averaging, \code{TRUE} to average per replicate or to the name of a
#'   column in the \link[=analysis-information]{analysis information} to compare by a custom grouping of analyses.
#'
#'   If \code{features=TRUE} all numerical properties are averaged and non-numericals are collapsed. Each row represents
#'   the feature aggregated data and the groups used for aggregation (\emph{e.g.} replicates) are stored in the
#'   \code{average_group} column. It is also possible to average these data per feature group, by setting
#'   \code{average="fGroups"}.
#' @param averageFunc Function used for averaging. Only used when data is averaged or \code{FCParams != NULL}.
#' @param features If \code{TRUE} then feature specific data will be added. Also see the \code{average} argument.
#' @param areas If set to \code{TRUE} then areas are output instead of peak intensities. Ignored if
#'   \code{features=TRUE}, as areas of features are always reported.
#' @param normalized If \code{TRUE} then normalized intensities are used (see the \verb{Feature intensity normalization}
#'   section). If no normalization data is available (\emph{e.g.} because \code{normInts} was not used) then an
#'   automatic group normalization is performed.
#' @param FCParams A parameter list to calculate Fold change data. See \code{getFCParams} for more details. Set to
#'   \code{NULL} to not perform FC calculations.
#' @param qualities Adds feature (group) qualities (\code{qualities="quality"}), scores (\code{qualities="score"}) or
#'   both (\code{qualities="both"}), if this data is available (\emph{i.e.} from \code{calculatePeakQualities}). If
#'   \code{qualities=FALSE} then nothing is reported.
#' @param regression,regressionBy Used for regression calculations. See the \verb{Regression calculation} section below.
#'   Set to \code{NULL} to ignore.
#' @param concAggrParams,toxAggrParams Parameters to aggregate calculated concentrations/toxicities (obtained with
#'   \code{\link{calculateConcs}}/\code{\link{calculateTox}}). See \link[=pred-aggr-params]{prediction aggregation
#'   parameters} for more information. Set to \code{NULL} to omit this data.
#' @param normConcToTox Set to \code{TRUE} to normalize concentrations to toxicities. Only relevant if this data is
#'   present (see \code{\link{calculateConcs}}/\code{\link{calculateTox}}).
#' @param anaInfoCols A \code{character} with any additional columns from the \link[=analysis-information]{analysis
#'   information} table. Only supported if \code{features=TRUE}. If averaging is performed then the data in the
#'   specified columns should be numeric. Set to \code{NULL} to ignore.
#' @param collapseSuspects If a \code{character} then any suspects that were matched to the same feature group are
#'   collapsed to a single row and suspect names are separated by the value of \code{collapseSuspects}. If \code{NULL}
#'   then no collapsing occurs, and each suspect match is reported on a single row. See the \verb{Suspect collapsing}
#'   section below for additional details.
#' @param onlyHits If \code{TRUE} then only feature groups with suspect hits are reported.
#' @param \dots Passed to the parent \code{as.data.table} method.
#'
#' @templateVar consider to be returned
#' @template IMS-arg
#'
#' @section Regression calculation: The \code{regression} argument controls the calculation of regression parameters
#'   from a regression model calculated with feature intensities (or areas if \code{areas=TRUE}). Here, simple linear
#'   regression is used, \emph{i.e.} \samp{y=ax+b} with \samp{a} the slope and \samp{b} the intercept. The value for
#'   \code{regression} should be the name of a column in the \link[=analysis-information]{analysis information} table
#'   with numerical data to be used for x-values. Alternatively, if \code{regression=TRUE} then the \code{"conc"} column
#'   is used. Any \code{NA} x-values are ignored, and no regression will be calculated if less than two (non-NA)
#'   x-values are available. The output table will contain properties such as the slope and correlation coefficient
#'   (R-squared). Furthermore, if \code{features=TRUE} then x-values will be calculated from the model and stored in the
#'   \code{x_reg} column.
#'
#'   The \code{regressionBy} argument can be used to construct separate regression models for different groups of
#'   analysis. It should be set to the name of a column in the \link[=analysis-information]{analysis information} table
#'   which defines the grouping between samples. If \code{features=TRUE} then the grouping is stored in the
#'   \code{regression_group} column of the output table.
#'
#'   Please see the handbook for examples on how to use the regression functionality.
#'
#' @section {Suspect collapsing}: The \code{as.data.table} method for \code{featureGroupsScreening} supports an
#'   additional format where each suspect hit is reported on a separate row (enabled by setting
#'   \code{collapseSuspects=NULL}). In this format the suspect
#'   properties from the \code{screenInfo} method are merged with each suspect row. Alternatively, if \emph{suspect
#'   collapsing} is enabled (the default) then the regular \code{as.data.table} format is used, and amended with the
#'   names and estimated ID levels (if available) of the suspects matched to a feature group (each separated by the
#'   value of the \code{collapseSuspects} argument).
#'
#'   Suspect collapsing also influences the reporting of predicted \link[=pred-quant]{feature concentrations} and
#'   \link[=pred-tox]{toxicities}. In the case that (1) suspects are \emph{not} collapsed in the output table and (2)
#'   predictions are available for a specific suspect hit (\emph{i.e.} if \code{\link{predictRespFactors}} or
#'   \code{\link{predictTox}} was called on the feature groups object), then only the suspect specific data is reported
#'   and no aggregation is performed. Hence, this allows you to obtain specific concentration/toxicity values for each
#'   suspect/feature group pair.
#'
#' @section IMS workflows: If the \code{IMS} argument is set to \code{"both"} or \code{"maybe"} then
#'   \code{"mobility_collapsed"} and \code{"CCS_collapsed"} columns will be added that summarize all
#'   mobility/\acronym{CCS} values of the mobility features (or feature groups) assigned to this IMS parent. These
#'   numbers are currently rounded to \samp{3} decimals.
#'
#'
#' @section Sets workflows: In a \link[=sets-workflow]{sets workflow} normalization of feature intensities occur per
#'   set.
#'
#'   In sets workflows the \link[=analysis-information]{analysis information} contains an additional \code{"set"}
#'   column, which can be used for arguments that involve grouping of analyses. For instance, if
#'   \code{regressionBy="set"} then regression models will be calculated for each set.
#'
#' @name feature-table
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
        if (!is.null(scrInfo[["susp_estIDLevel"]]))
            scrInfo[, susp_bestEstIDLevel := getBestIDLevel(scrInfo$susp_estIDLevel), by = "susp_group"]
        # only keep unique and remove suspect specific columns
        scrInfo <- unique(subsetDTColumnsIfPresent(scrInfo, c("susp_group", "susp_name", "susp_bestEstIDLevel")),
                          by = "susp_group")
    }
    
    ret <- merge(tab, scrInfo, by.x = "group", by.y = "susp_group", all.x = !onlyHits, sort = FALSE,
                 allow.cartesian = is.null(collapseSuspects))
    return(ret)
}

doFGAADTGroups <- function(fGroups, intColNames, average, averageBy, areas, addQualities, addScores, regression,
                           regressionBy, averageFunc, normalized, FCParams, concAggrParams, toxAggrParams,
                           normConcToTox, IMS, collapseSuspects, onlyHits)
{
    fGroupsOrig <- fGroups
    if (IMS != "both")
        fGroups <- selectIMSFilter(fGroups, IMS = IMS, verbose = FALSE)
    
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
        gTableAvg <- averageGroups(fGroups, areas, normalized, by = "replicate", func = averageFunc) # UNDONE: support averageBy
        
        repInds <- match(FCParams$replicates, replicates(fGroups))
        for (i in seq_along(gTableAvg))
            set(ret, i, "FC", do.call(calcFC, as.list(gTableAvg[[i]][repInds])))
        ret[, FC_log := log2(FC)]
        
        anaInds1 <- which(anaInfo$replicate %in% FCParams$replicates[1])
        anaInds2 <- which(anaInfo$replicate %in% FCParams$replicates[2])
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
    
    ret[, c("group", "ret", "mz") := .(gNames, gInfo$ret, gInfo$mz)]
    if (hasMobilities(fGroups))
    {
        cols <- intersect(c("mobility", "CCS", "ims_parent_group"), names(gInfo))
        ret[, (cols) := gInfo[, cols, with = FALSE]]
        if (!isTRUE(IMS))
        {
            gInfoOrig <- groupInfo(fGroupsOrig)
            doCollapse <- function(coln, digits)
            {
                if (!is.null(ret[[coln]]))
                {
                    ret[!is.na(get(coln)), (paste0(coln, "_collapsed")) := as.character(round(get(coln), digits))]
                    ret[is.na(get(coln)), (paste0(coln, "_collapsed")) := sapply(group, function(g)
                    {
                        paste0(round(gInfoOrig[ims_parent_group == g][[coln]], digits), collapse = ",")
                    })]
                }
            }
            # UNDONE: is the rounding with a good number?
            # NOTE: update docs if this is changed
            doCollapse("mobility", 3)
            doCollapse("CCS", 3)
        }
    }
    setcolorder(ret, intersect(c("group", "ims_parent_group", "ret", "mz", "mobility", "mobility_collapsed", "CCS",
                                 "CCS_collapsed"), names(ret)))
    
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
                             regressionBy, averageFunc, anaInfoCols, IMS)
{
    fGroupsOrig <- NULL
    if (IMS != "both")
    {
        fGroupsOrig <- fGroups
        fGroups <- selectIMSFilter(fGroups, IMS = IMS, verbose = FALSE)
    }
    
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
        featTab[, replicate := anaInfo$replicate[match(analysis, anaInfo$analysis)]]
        featTab[, average_group := analysis]
    }
    else
    {
        if (averageBy != "fGroups")
            featTab[, average_group := anaInfo[[averageBy]][match(analysis, anaInfo$analysis)]]
        
        # average all numeric columns
        numCols <- setdiff(names(which(sapply(featTab, is.numeric))), "average_group")
        # avoid DT warning when averaging int values that lead to decimal values
        # UNDONE: or collapse int values instead: for piek they are scan numbers, OpenMS isocounts, so could make sense
        featTab[, (numCols) := lapply(.SD, as.numeric), .SDcols = numCols]
        featTab[, (numCols) := lapply(.SD, averageFunc), .SDcols = numCols, by = by]
        
        # collapse all character columns
        chCols <- names(which(sapply(featTab, is.character)))
        featTab[, (chCols) := lapply(.SD, function(x) paste0(unique(x), collapse = ",")), .SDcols = chCols, by = by]
        featTab <- unique(featTab, by = by)
    }
    
    if (IMS != "both" && hasMobilities(fGroups))
    {
        fTabOrig <- featureTable(fGroupsOrig)
        
        # UNDONE: is this a good number?
        # NOTE: update docs if this is changed
        digits <- 3
        
        featTab[!is.na(mobility), mobility_collapsed := as.character(round(mobility, digits))]
        featTab[is.na(mobility), mobility_collapsed := mapply(ID, analysis, FUN = function(i, a)
        {
            paste0(round(fTabOrig[[a]][ims_parent_ID == i]$mobility, digits), collapse = ",")
        })]
    }
    
    # prepare main table for merge
    
    if (averageBy != "fGroups")
    {
        # Melt by intensity column to get the proper format. Afterwards, we remove the dummy intensity column, as we
        # want the raw feature intensity data.
        
        # If there is also concentration data, then we melt both the intensity and concentration columns.
        concCols <- paste0(stripADTIntSuffix(intColNames), "_conc")
        if (all(concCols %in% names(fgTab)))
        {
            fgTab <- melt(fgTab, measure.vars = list(intensity = intColNames, conc = concCols),
                          variable.name = "average_group", variable.factor = FALSE)
            # DT changes the analyses names to (character) indices with >1 measure vars: https://github.com/Rdatatable/data.table/issues/4047
            fgTab[, average_group := intColNames[as.integer(average_group)]]
        }
        else
        {
            # NOTE: we recent data.table versions, we cannot use measure.vars with a list of length 1:
            # https://github.com/Rdatatable/data.table/issues/5209
            fgTab <- melt(fgTab, measure.vars = intColNames, variable.name = "average_group", variable.factor = FALSE,
                          value.name = "intensity")
        }
        
        fgTab <- fgTab[intensity != 0]
        fgTab[, average_group := stripADTIntSuffix(average_group)]
    }
    fgTab <- removeDTColumnsIfPresent(fgTab, "intensity")
    setnames(fgTab, c("ret", "mz", "mobility", "mobility_collapsed", "CCS", "CCS_collapsed"),
             c("group_ret", "group_mz", "group_mobility", "group_mobility_collapsed", "group_CCS", "group_CCS_collapsed"),
             skip_absent = TRUE)
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
    colord <- c("group", "ims_parent_group", "set", "analysis", "average_group", "replicate", "group_ret",
                "group_mz", "group_mobility", "group_mobility_collapsed", "ID", "ims_parent_ID", "ret", "mz", "ion_mz",
                "mobility", "mobility_collapsed", "intensity", "area", "intensity_rel", "area_rel")
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
                            IMS = "both", collapseSuspects = ",", onlyHits = FALSE)
{
    anaInfo <- analysisInfo(fGroups)
    if (isTRUE(regression))
        regression <- "conc" # legacy
    
    assertFGAsDataTableArgs(fGroups, areas, features, qualities, regression, regressionBy, averageFunc, normalized,
                            FCParams, concAggrParams, toxAggrParams, normConcToTox, anaInfoCols, IMS, collapseSuspects,
                            onlyHits)
    averageBy <- assertAndPrepareAnaInfoBy(average, anaInfo, TRUE)
    
    if (length(fGroups) == 0)
        return(data.table(mz = numeric(), ret = numeric(), group = character()))

    intColNames <- if (averageBy == "fGroups") "intensity" else getADTIntCols(unique(anaInfo[[averageBy]]))    
    if (!isFALSE(regression) && (sum(!is.na(anaInfo[[regression]])) < 2 || averageBy == "fGroups"))
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
                          normConcToTox, IMS, collapseSuspects, onlyHits)    

    if (features)
    {
        ret <- doFGAADTFeatures(fGroups, ret, intColNames, average, averageBy, addQualities, addScores, regression,
                                regressionBy, averageFunc, anaInfoCols, IMS)
    }
    
    return(ret[])
}

# both non-sets and sets
doFGScrAsDataTable <- function(x, ..., collapseSuspects = ",", onlyHits = FALSE)
{
    return(doFGAsDataTable(x, ..., collapseSuspects = collapseSuspects, onlyHits = onlyHits))
}


#' @rdname feature-table
#' @export
setMethod("as.data.table", "featureGroups", function(x, average = FALSE, areas = FALSE, features = FALSE,
                                                     qualities = FALSE, regression = FALSE, regressionBy = NULL,
                                                     averageFunc = mean, normalized = FALSE, FCParams = NULL,
                                                     concAggrParams = getDefPredAggrParams(),
                                                     toxAggrParams = getDefPredAggrParams(), normConcToTox = FALSE,
                                                     anaInfoCols = NULL, IMS ="both")
{
    return(doFGAsDataTable(x, average, areas, features, qualities, regression, regressionBy, averageFunc, normalized,
                           FCParams, concAggrParams, toxAggrParams, normConcToTox, anaInfoCols))
})

#' @rdname feature-table
#' @export
setMethod("as.data.table", "featureGroupsScreening", doFGScrAsDataTable)

#' @rdname feature-table
#' @export
setMethod("as.data.table", "featureGroupsScreeningSet", doFGScrAsDataTable)
