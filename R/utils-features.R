# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

makeFGroupName <- function(id, ret, mz) sprintf("M%.0f_R%.0f_%d", mz, ret, id)

showAnaInfo <- function(anaInfo)
{
    rGroups <- unique(anaInfo$group)
    blGroups <- unique(anaInfo$blank)
    printf("Analyses: %s (%d total)\n", getStrListWithMax(anaInfo$analysis, 6, ", "), nrow(anaInfo))
    printf("Replicate groups: %s (%d total)\n", getStrListWithMax(rGroups, 8, ", "), length(rGroups))
    printf("Replicate groups used as blank: %s (%d total)\n", getStrListWithMax(blGroups, 8, ", "), length(blGroups))
}

reGenerateFTIndex <- function(fGroups)
{
    gNames <- names(fGroups)
    fGroups@ftindex <- setnames(rbindlist(lapply(featureTable(fGroups),
                                                 function(ft) as.list(chmatch(gNames, ft$group, 0)))), gNames)
    return(fGroups)
}
isFGSet <- function(fGroups) inherits(fGroups, "featureGroupsSet")

featureQualities <- function()
{
    checkPackage("MetaClean")
    list(ApexBoundaryRatio = list(func = MetaClean::calculateApexMaxBoundaryRatio, HQ = "LV", range = c(0, 1)),
         FWHM2Base = list(func = MetaClean::calculateFWHM, HQ = "HV", range = c(0, 1)),
         Jaggedness = list(func = MetaClean::calculateJaggedness, HQ = "LV", range = Inf),
         Modality = list(func = MetaClean::calculateModality, HQ = "LV", range = Inf),
         Symmetry = list(func = MetaClean::calculateSymmetry, HQ = "HV", range = c(-1, 1)),
         GaussianSimilarity = list(func = MetaClean::calculateGaussianSimilarity, HQ = "HV", range = c(0, 1)),
         Sharpness = list(func = MetaClean::calculateSharpness, HQ = "HV", range = Inf),
         TPASR = list(func = MetaClean::calculateTPASR, HQ = "LV", range = Inf),
         ZigZag = list(func = MetaClean::calculateZigZagIndex, HQ = "LV", range = Inf))
}

featureGroupQualities <- function()
{
    checkPackage("MetaClean")
    list(
        ElutionShift = list(func = MetaClean::calculateElutionShift, HQ = "LV", range = Inf),
        RetentionTimeCorrelation = list(func = MetaClean::calculateRetentionTimeConsistency, HQ = "LV", range = Inf)
    )
}

# normalize, invert if necessary to get low (worst) to high (best) order
scoreFeatQuality <- function(quality, values)
{
    if (all(is.finite(quality$range)))
    {
        if (!isTRUE(all.equal(quality$range, c(0, 1)))) # no need to normalize 0-1
            values <- normalize(values, minMax = quality$range[1] < 0, xrange = quality$range)
    }
    else
        values <- normalize(values, TRUE)
    
    if (quality$HQ == "LV")
        values <- 1 - values
    
    return(values)
}

hasFGroupScores <- function(fGroups) nrow(groupScores(fGroups)) > 0

getFilteredFGroups <- function(fGroups, retFilter) 
{
    if (!is.null(retFilter$subset))
    {
        if (is.null(retFilter$subset$i) && !is.null(retFilter$subset$j))
            fGroups <- fGroups[, retFilter$subset$j]
        else if (!is.null(retFilter$subset$i) && is.null(retFilter$subset$j))
            fGroups <- fGroups[retFilter$subset$i]
        else
            fGroups <- fGroups[retFilter$subset$i, retFilter$subset$j]  
    }
    if (!is.null(retFilter$delete))
        fGroups <- delete(fGroups, i = retFilter$delete$i, j = retFilter$delete$j)
    return(fGroups)
}

doFGroupsFilter <- function(fGroups, what, hashParam, func, cacheCateg = what, verbose = TRUE)
{
    if (verbose)
    {
        printf("Applying %s filter... ", what)
        oldFCount <- length(getFeatures(fGroups)); oldGCount <- length(fGroups)
    }
    
    cacheName <- sprintf("filterFGroups_%s", cacheCateg)
    hash <- makeHash(fGroups, what, hashParam)
    ret <- loadCacheData(cacheName, hash)
    if (is.null(ret))
    {
        ret <- if (length(fGroups) > 0) func(fGroups) else NULL
        if (!is.null(ret))
            saveCacheData(cacheName, ret, hash)
    }

    fGroups <- getFilteredFGroups(fGroups, ret)
    
    if (verbose)
    {
        newFCount <- length(getFeatures(fGroups)); newGCount <- length(fGroups)
        newn <- length(ret)
        printf("Done! Filtered %d (%.2f%%) features and %d (%.2f%%) feature groups. Remaining: %d features in %d groups.\n",
               oldFCount - newFCount, if (oldFCount > 0) (1 - (newFCount / oldFCount)) * 100,
               oldGCount - newGCount, if (oldGCount > 0) (1 - (newGCount / oldGCount)) * 100,
               newFCount, newGCount)
    }
    
    return(fGroups)
}

# this combines all functionality from all fGroup as.data.table methods, a not so pretty but pragmatic solution...
doFGAsDataTable <- function(fGroups, average = FALSE, areas = FALSE, features = FALSE, qualities = FALSE,
                            regression = FALSE, averageFunc = mean, normalized = FALSE, FCParams = NULL,
                            concAggrParams = getDefPredAggrParams(), toxAggrParams = getDefPredAggrParams(),
                            normConcToTox = FALSE, collapseSuspects = ",", onlyHits = FALSE)
{
    assertFGAsDataTableArgs(fGroups, average, areas, features, qualities, regression, averageFunc, normalized, FCParams,
                            concAggrParams, toxAggrParams, normConcToTox, collapseSuspects, onlyHits)
    
    if (features && average && regression)
        stop("Cannot add regression data for averaged features.")
    if (features && !is.null(FCParams))
        stop("Cannot calculate fold-changes with features=TRUE")

    if (length(fGroups) == 0)
        return(data.table(mz = numeric(), ret = numeric(), group = character()))
    
    anaInfo <- analysisInfo(fGroups)
    gNames <- names(fGroups)
    gInfo <- groupInfo(fGroups)
    doConc <- regression && !is.null(anaInfo[["conc"]]) && sum(!is.na(anaInfo[["conc"]]) > 1)
    addQualities <- !isFALSE(qualities) && qualities %in% c("both", "quality") && hasFGroupScores(fGroups)
    addScores <- !isFALSE(qualities) && qualities %in% c("both", "score") && hasFGroupScores(fGroups)
    
    if (regression && is.null(anaInfo[["conc"]]))
        warning("No concentration information specified in the analysis information (i.e. conc column, see ?`analysis-information`)")
    
    if (features)
    {
        ftindex <- groupFeatIndex(fGroups)
        fTable <- featureTable(fGroups)
        
        ret <- rbindlist(fTable, idcol = "analysis")
        setorder(ret, "group")
        
        if (doConc)
            ret[, conc := anaInfo$conc[match(analysis, anaInfo$analysis)]]
        
        if (!is.null(ret[["adduct"]]))
            ret[, adduct := NULL] # we already include group annotations below
        
        if (addQualities)
        {
            gq <- groupQualities(fGroups)[match(ret$group, group), -"group"]
            ret[, (paste0("group_", names(gq))) := gq]
        }
        else if (hasFGroupScores(fGroups))
            ret[, (intersect(featureQualityNames(group = FALSE), names(ret))) := NULL]
        if (addScores)
        {
            gs <- groupScores(fGroups)[match(ret$group, group), -"group"]
            ret[, (paste0("group_", names(gs))) := gs]
        }
        else if (hasFGroupScores(fGroups))
            ret[, (intersect(featureQualityNames(group = FALSE, scores = TRUE), names(ret))) := NULL]
        
        if (average)
        {
            ret <- removeDTColumnsIfPresent(ret, c("isocount", "analysis", "ID"))
            numCols <- names(which(sapply(ret, is.numeric)))
            ret[, (numCols) := lapply(.SD, averageFunc), .SDcols = numCols, by = "group"]
            ret <- unique(ret, by = "group")
        }
        else
        {
            ret[, replicate_group := anaInfo$group[match(analysis, anaInfo$analysis)]]
            
            doConc <- doConc && nrow(anaInfo) > 1
            if (doConc)
            {
                ret[, c("RSQ", "intercept", "slope", "p") := {
                    notna <- !is.na(conc)
                    if (!any(notna))
                        NA_real_
                    else
                    {
                        ints <- intensity[notna]
                        ints[ints == 0] <- NA
                        if (all(is.na(ints)))
                            NA_real_
                        else
                        {
                            suppressWarnings(reg <- summary(lm(ints ~ conc[notna])))
                            slope <- pv <- NA_real_
                            if (nrow(reg[["coefficients"]]) > 1)
                            {
                                slope <- reg[["coefficients"]][2, 1]
                                pv <- reg[["coefficients"]][2, 4]
                            }
                            list(reg[["r.squared"]], reg[["coefficients"]][1, 1], slope, pv)
                        }
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
        if (normalized)
            fGroups <- maybeAutoNormalizeFGroups(fGroups)
        
        gTableAvg <- averageGroups(fGroups, areas, normalized, func = averageFunc)
        gTableNonAvg <- groupTable(fGroups, areas, normalized)
        
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
            regr <- sapply(gTable, function(grp)
            {
                grp[grp == 0] <- NA
                if (all(is.na(grp[notna])))
                    return(NULL)
                return(summary(lm(grp[notna] ~ notnaconcs)))
            }, simplify = FALSE)
            notNULL <- !sapply(regr, is.null)
            ret[notNULL, c("RSQ", "intercept", "slope", "p") :=
                    .(sapply(regr[notNULL], "[[", "r.squared"),
                      sapply(regr[notNULL], function(r) r$coefficients[1, 1]),
                      sapply(regr[notNULL], function(r) if (nrow(r$coefficients) > 1) r$coefficients[2, 1] else NA_real_),
                      sapply(regr[notNULL], function(r) if (nrow(r$coefficients) > 1) r$coefficients[2, 4] else NA_real_))]
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
    }

    if (isFGSet(fGroups) && features && !average) # add set column if feature data is present
    {
        ret[, set := anaInfo[match(analysis, anaInfo$analysis), "set"]]
        setcolorder(ret, c("group", "group_ret", "group_mz", "set", "analysis"))
    }
    
    annTable <- annotations(fGroups)
    if (nrow(ret) > 0 && nrow(annTable) > 0)
    {
        if (isFGSet(fGroups))
        {
            if (features && !average)
                ret <- merge(ret, annTable, by = c("group", "set"), sort = FALSE)
            else
            {
                # collapse annotation info for each group
                annTable <- copy(annTable)
                annTable[, adduct := paste0(adduct, collapse = ","), by = "group"]
                annTable <- unique(annTable, by = "group")[, -"set"]
                ret <- merge(ret, annTable, by = "group", sort = FALSE)
            }
        }
        else
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
        #   - if features==T then take care to melt the concs to the same format as ret, and then merge by group/analysis or group/suspect/analysis
        
        concs <- aggregateConcs(concentrations(fGroups), anaInfo, concAggrParams, splitSusps)
        setnames(concs, "type", "conc_types")

        if (features)
        {
            concs <- melt(concs, measure.vars = anaInfo$analysis, variable.name = "analysis", value.name = "conc")
            if (average)
            {
                concs[, conc := if (all(is.na(conc))) NA_real_ else averageFunc(conc[!is.na(conc)]), by = "group"]
                concs[, analysis := NULL]
                concs <- unique(concs, by = "group")
            }
        }
        else if (average)
        {
            for (rg in replicateGroups(fGroups))
            {
                anas <- anaInfo[anaInfo$group == rg, "analysis"]
                concs[, (paste0(rg, "_conc")) := aggrVec(unlist(.SD), averageFunc), .SDcols = anas, by = seq_len(nrow(concs))]
            }
            concs[, (anaInfo$analysis) := NULL]
        }
        else
            setnames(concs, anaInfo$analysis, paste0(anaInfo$analysis, "_conc"))
        
        setcolorder(concs, setdiff(names(concs), "conc_types")) # move to end

        mby <- "group"
        if (features && !average)
            mby <- c(mby, "analysis")
        ret <- mergePred(ret, concs, mby, "conc_types")
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

# used by adducts()<- methods
updateAnnAdducts <- function(annTable, gInfo, adducts)
{
    if (nrow(annTable) > 0)
    {
        # only consider changed
        adducts <- adducts[adducts != annTable$adduct]
        
        if (length(adducts) == 0)
            return(annTable) # nothing changed
    }
    
    adducts <- sapply(adducts, checkAndToAdduct, .var.name = "value", simplify = FALSE)
    adductsChr <- sapply(adducts, as.character) # re-make characters: standardize format
    nm <- calculateMasses(gInfo[names(adducts), "mzs"], adducts, type = "neutral")
    
    if (nrow(annTable) > 0)
    {
        # update table
        annTable <- copy(annTable)
        annTable[match(names(adducts), group), c("adduct", "neutralMass") := .(adductsChr, nm)][]
    }
    else # initialize table
        annTable <- data.table(group = names(adducts), adduct = adductsChr, neutralMass = nm)
    
    return(annTable)
}

unSetAnaInfo <- function(anaInfo) anaInfo[, setdiff(names(anaInfo), "set")]

maybeAutoNormalizeFGroups <- function(fGroups)
{
    if (length(fGroups) == 0 || !is.null(featureTable(fGroups)[[1]][["intensity_rel"]]))
        return(fGroups) # no features or already normalized
    printf("Automatically normalizing feature groups, see ?normInts for more options.\n")
    return(normInts(fGroups, featNorm = "none", groupNorm = TRUE))
}

getAnnotationsFromSetFeatures <- function(fGroups)
{
    if (length(fGroups) > 0)
    {
        anaInfo <- analysisInfo(fGroups)
        ftind <- groupFeatIndex(fGroups)
        fTable <- featureTable(fGroups)
        
        ret <- rbindlist(sapply(sets(fGroups), function(s)
        {
            anaInds <- which(anaInfo$set == s)
            anas <- anaInfo[anaInds, "analysis"]
            grps <- names(fGroups)[sapply(ftind[anaInds], function(x) any(x != 0))]
            firstFeats <- rbindlist(lapply(ftind[anaInds, grps, with = FALSE], function(x)
            {
                firstAna <- which(x != 0)[1]
                return(fTable[[anas[firstAna]]][x[firstAna]])
            }))
            
            return(data.table(group = grps, adduct = firstFeats$adduct))
        }, simplify = FALSE), idcol = "set", fill = TRUE) # set fill for empty objects
        ret[, neutralMass := groupInfo(fGroups)[ret$group, "mzs"]]
    }
    else
        ret <- data.table()
}

finishFGroupsForSets <- function(fGroups, ..., verbose)
{
    otherArgs <- list(...)
    feat <- getFeatures(fGroups) # may have been changed (eg in initialize())
    ret <- featureGroupsSet(groupAlgo = algorithm(fGroups), groupArgs = otherArgs, groupVerbose = verbose,
                            groups = groupTable(fGroups), groupInfo = groupInfo(fGroups),
                            analysisInfo = analysisInfo(fGroups), features = feat, ftindex = groupFeatIndex(fGroups),
                            algorithm = makeSetAlgorithm(list(fGroups)))
    ret@annotations <- getAnnotationsFromSetFeatures(ret)
    return(ret)
}

groupFeaturesSets <- function(feat, grouper, ..., verbose)
{
    fGroups <- selectMethod(grouper, "features")(feat = feat, ...)
    return(finishFGroupsForSets(fGroups, ..., verbose = verbose))
}

filterEICs <- function(EICs, fGroups, analysis = NULL, groupName = NULL, topMost = NULL, topMostByRGroup = FALSE,
                       onlyPresent = FALSE)
{
    if (!is.null(analysis))
        EICs <- EICs[names(EICs) %chin% analysis]
    else
        analysis <- analyses(fGroups)
    if (!is.null(groupName))
        EICs <- lapply(EICs, function(e) e[names(e) %chin% groupName])

    if (onlyPresent)
    {
        gTable <- groupTable(fGroups)
        EICs <- Map(names(EICs), EICs, f = function(ana, aeic)
        {
            anaInd <- match(ana, analyses(fGroups))
            absentFGs <- names(gTable)[gTable[anaInd] == 0]
            return(aeic[!names(aeic) %chin% absentFGs])
        })
    }
    
    if (!is.null(topMost))
    {
        gTable <- copy(groupTable(fGroups))
        gTable[, c("analysis", "rGroup") := analysisInfo(fGroups)[, c("analysis", "group")]]
        for (fg in names(fGroups))
        {
            anasWithFG <- Map(names(EICs), EICs, f = function(ana, aeic) if (fg %chin% names(aeic)) ana else character())
            anasWithFG <- pruneList(anasWithFG, checkEmptyElements = TRUE)
            anasWithFG <- unlist(anasWithFG)
            tab <- gTable[analysis %chin% anasWithFG, c(fg, "analysis", "rGroup"), with = FALSE]
            if (nrow(tab) > 0)
            {
                rmAnas <- if (topMostByRGroup)
                {
                    tab[, rank := frank(-get(fg), ties.method = "first"), by = "rGroup"]
                    tab[rank > topMost]$analysis
                }
                else if (nrow(tab) > topMost)
                {
                    setorderv(tab, fg, order = -1L)
                    tab[seq(topMost + 1, nrow(tab))]$analysis
                }
                else
                    character()
                EICs[rmAnas] <- lapply(EICs[rmAnas], "[[<-", fg, NULL)
            }
        }
    }

    return(pruneList(EICs, checkEmptyElements = TRUE))
}

convertConc <- function(conc, unitFrom, unitTo, MW)
{
    parseUnit <- function(unit)
    {
        logBase <- if (startsWith(unit, "log10 "))
            10
        else if (startsWith(unit, "log2 "))
            2
        else if (startsWith(unit, "log "))
            exp(1)
        else
            0
        if (logBase != 0)
            unit <- sub("^log(2|10)? ", "", unit)
        base <- switch(substring(unit, 1, 1),
                       n = 1E-9,
                       u = 1E-6,
                       m = 1E-3,
                       1)
        molar <- endsWith(unit, "M")
        return(list(logBase = logBase, base = base, molar = molar))
    }
    
    unitFromP <- parseUnit(unitFrom); unitToP <- parseUnit(unitTo)
    
    if (unitFromP$logBase != 0 && unitToP$logBase == 0)
        conc <- unitFromP$logBase^conc
    else if (unitFromP$logBase == 0 && unitToP$logBase)
        conc <- log(conc, unitToP$logBase)
    else if (unitFromP$logBase != unitToP$logBase)
        conc <- log(unitFromP$logBase^conc, unitToP$logBase)
    
    conc <- conc * (unitFromP$base / unitToP$base)
    
    if (unitFromP$molar && !unitToP$molar)
        conc <- conc * MW
    else if (!unitFromP$molar && unitToP$molar)
        conc <- conc / MW
    
    return(conc)
}

calcFeatureConcs <- function(fGroups, resp, areas)
{
    # get concentration data from response factors
    
    gt <- transpose(groupTable(fGroups, areas = areas)[, resp$group, with = FALSE])
    setnames(gt, analyses(fGroups))
    
    concs <- copy(resp)
    concs[, (analyses(fGroups)) := lapply(gt, function(ints) ints / RF)]
    concs[, RF := NULL]

    return(concs[])
}

finalizeFeaturePredTab <- function(pred)
{
    co <- intersect(c("group", "type", "candidate", "candidate_name", "sets"), names(pred))
    setcolorder(pred, co)
    setorderv(pred, c("group", "type", "candidate"))
    return(pred)
}

doCalcConcSets <- function(fGroups, featureAnn, areas)
{
    if (is.null(featureAnn))
    {
        # just from screening info
        usFGroups <- sapply(sets(fGroups), unset, obj = fGroups, simplify = FALSE)
        usFGroups <- sapply(usFGroups, calculateConcs, areas = areas, simplify = FALSE)
    }
    else
    {
        usFeatAnns <- checkAndUnSetOther(sets(fGroups), featureAnn, "featureAnn")
        usFeatAnns <- usFeatAnns[lengths(usFeatAnns) > 0]
        
        if (length(usFeatAnns) == 0)
        {
            cat("No feature annotations, nothing to do...\n")
            return(fGroups)
        }
        
        usFGroups <- sapply(names(usFeatAnns), unset, obj = fGroups, simplify = FALSE)
        usFGroups <- Map(usFGroups, usFeatAnns, f = calculateConcs, MoreArgs = list(areas = areas))
    }
    
    usFGroups <- usFGroups[sapply(usFGroups, function(ufg) nrow(concentrations(ufg)) > 0)]
    if (length(usFGroups) == 0)
        return(fGroups)
    
    if (length(usFGroups) == 1)
        fGroups@concentrations <- copy(concentrations(usFGroups[[1]]))
    else
    {
        fGroups@concentrations <- Reduce(lapply(usFGroups, concentrations), f = function(left, right)
        {
            cols <- intersect(c("group", "type", "candidate", "candidate_name"), c(names(left), names(right)))
            return(merge(left, right, by = cols, all = TRUE))
        })
        fGroups@concentrations <- finalizeFeaturePredTab(fGroups@concentrations)
    }
    
    # add missing analyses
    missingAnaCols <- setdiff(analyses(fGroups), names(fGroups@concentrations))
    if (length(missingAnaCols) > 0)
        fGroups@concentrations[, (missingAnaCols) := NA_real_][]
    
    return(fGroups)
}

doCalcToxSets <- function(fGroups, featureAnn)
{
    if (is.null(featureAnn))
    {
        # just from screening info
        usFGroups <- sapply(sets(fGroups), unset, obj = fGroups, simplify = FALSE)
        usFGroups <- sapply(usFGroups, calculateTox, simplify = FALSE)
    }
    else
    {
        usFeatAnns <- checkAndUnSetOther(sets(fGroups), featureAnn, "featureAnn")
        usFeatAnns <- usFeatAnns[lengths(usFeatAnns) > 0]
        
        if (length(usFeatAnns) == 0)
        {
            cat("No feature annotations, nothing to do...\n")
            return(fGroups)
        }
        
        usFGroups <- sapply(names(usFeatAnns), unset, obj = fGroups, simplify = FALSE)
        usFGroups <- Map(usFGroups, usFeatAnns, f = calculateTox)
    }
    
    usFGroups <- usFGroups[sapply(usFGroups, function(ufg) nrow(toxicities(ufg)) > 0)]
    if (length(usFGroups) == 0)
        return(fGroups)
    
    if (length(usFGroups) == 1)
    {
        fGroups@toxicities <- copy(toxicities(usFGroups[[1]]))
        fGroups@toxicities[, sets := names(usFGroups)]
    }
    else
    {
        allTox <- sapply(usFGroups, toxicities, simplify = FALSE)
        
        allToxTab <- rbindlist(allTox, idcol = "sets")
        
        # NOTE: only for SIRIUS_FP it makes sense to split sets data, the rest we can merge
        allToxTab[type != "SIRIUS_FP", sets := paste0(unique(sets), collapse = ","), by = c("group", "candidate", "type")]
        allToxTab <- unique(allToxTab, by = c("group", "type", "candidate", "sets"))
        
        fGroups@toxicities <- finalizeFeaturePredTab(allToxTab)
    }
    
    return(fGroups)
}

aggregateConcs <- function(concs, anaInfo, aggrParams, splitSuspects = FALSE)
{
    concs <- copy(concs)
    concs <- subsetDTColumnsIfPresent(concs, c("group", "type", "candidate", "candidate_name", anaInfo$analysis))
    
    if (splitSuspects && !any(concs$type == "suspect"))
        splitSuspects <- FALSE # no suspects, nothing to split
    
    if (aggrParams$preferType != "none")
    {
        concs[, keep := !aggrParams$preferType %in% type | type == aggrParams$preferType, by = "group"]
        concs <- concs[keep == TRUE][, keep := NULL]
    }
    
    ignoreFGs <- character()
    if (splitSuspects)
    {
        # fgs with suspect results will be reported for each suspect per row and results from other types will be
        # ignored
        fgsWithSuspect <- unique(concs[type == "suspect"]$group)
        concs <- concs[!group %chin% fgsWithSuspect | type == "suspect"] # remove all others
        ignoreFGs <- fgsWithSuspect # no further aggregation needed
    }
    else if (!is.null(concs[["candidate_name"]]))
        concs[, candidate_name := NULL]

    doAggr <- function(func, by, combineTypes = FALSE)
    {
        concs[!group %chin% ignoreFGs, (anaInfo$analysis) := lapply(.SD, aggrVec, func),
              .SDcols = anaInfo$analysis, by = by]
        dups <- duplicated(concs, by = by)
        if (combineTypes)
            concs[, type := paste0(unique(type), collapse = ","), by = by]            
        return(concs[group %chin% ignoreFGs | dups == FALSE])
    }
    concs <- doAggr(aggrParams$candidateFunc, c("group", "type", "candidate"))
    concs <- doAggr(aggrParams$typeFunc, c("group", "type"))
    concs <- doAggr(aggrParams$groupFunc, "group", combineTypes = TRUE)
    
    mby <- if (splitSuspects) c("group", "candidate_name") else "group"
    concs <- unique(concs, by = mby)
    
    concs[, candidate := NULL]
    
    return(concs[])
}

aggregateTox <- function(tox, aggrParams, splitSuspects = FALSE)
{
    tox <- copy(tox)
    tox <- subsetDTColumnsIfPresent(tox, c("group", "type", "candidate", "candidate_name", "LC50"))
    
    if (splitSuspects && !any(tox$type == "suspect"))
        splitSuspects <- FALSE # no suspects, nothing to split

    if (aggrParams$preferType != "none")
    {
        tox[, keep := !aggrParams$preferType %in% type | type == aggrParams$preferType, by = "group"]
        tox <- tox[keep == TRUE][, keep := NULL]
    }
    
    ignoreFGs <- character()
    if (splitSuspects)
    {
        # fgs with suspect results will be reported for each suspect per row and results from other types will be
        # ignored
        fgsWithSuspect <- unique(tox[type == "suspect"]$group)
        tox <- tox[!group %chin% fgsWithSuspect | type == "suspect"] # remove all others
        ignoreFGs <- fgsWithSuspect # no further aggregation needed
    }
    else if (!is.null(tox[["candidate_name"]]))
        tox[, candidate_name := NULL]

    doAggr <- function(func, by, combineTypes = FALSE)
    {
        tox[!group %chin% ignoreFGs, LC50 := aggrVec(LC50, func), by = by]
        dups <- duplicated(tox, by = by)
        if (combineTypes)
            tox[, type := paste0(unique(type), collapse = ","), by = by]            
        return(tox[group %chin% ignoreFGs | dups == FALSE])
    }
    tox <- doAggr(aggrParams$candidateFunc, c("group", "type", "candidate"))
    tox <- doAggr(aggrParams$typeFunc, c("group", "type"))
    tox <- doAggr(aggrParams$groupFunc, "group", combineTypes = TRUE)

    mby <- if (splitSuspects) c("group", "candidate_name") else "group"
    tox <- unique(tox, by = mby)
    
    tox[, candidate := NULL]
    
    return(tox[])
}
