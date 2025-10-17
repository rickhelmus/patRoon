# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

isFGSet <- function(fGroups) inherits(fGroups, "featureGroupsSet")

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
            anas <- anaInfo$analysis[anaInds]
            grps <- names(fGroups)[sapply(ftind[anaInds], function(x) any(x != 0))]
            firstFeats <- rbindlist(lapply(ftind[anaInds, grps, with = FALSE], function(x)
            {
                firstAna <- which(x != 0)[1]
                return(fTable[[anas[firstAna]]][x[firstAna]])
            }))
            
            return(data.table(group = grps, adduct = firstFeats$adduct))
        }, simplify = FALSE), idcol = "set", fill = TRUE) # set fill for empty objects
        ret[, neutralMass := groupInfo(fGroups)$mz[match(ret$group, groupInfo(fGroups)$group)]]
        adducts <- sapply(unique(ret$adduct), as.adduct)
        ret[, ion_mz := calculateMasses(neutralMass, adducts[adduct], type = "mz")]
    }
    else
        ret <- data.table()
    return(ret)
}

setMethod("doGroupFeatures", "featuresSet", function(feat, grouper, groupAlgo, ..., IMSWindow, verbose)
{
    fg <- callNextMethod()
    
    otherArgs <- list(...)
    feat <- getFeatures(fg) # may have been changed (eg in initialize())
    ret <- featureGroupsSet(groupAlgo = algorithm(fg), groupArgs = otherArgs, groupVerbose = verbose,
                            groups = groupTable(fg), groupInfo = groupInfo(fg), features = feat,
                            ftindex = groupFeatIndex(fg), algorithm = makeSetAlgorithm(list(fg)))
    ret@annotations <- getAnnotationsFromSetFeatures(ret)
    
    return(ret)
})

setMethod("getFeatureEIXInputTab", "featureGroupsSet", function(obj, type, analysis, groupName, EIXParams)
{
    ret <- callNextMethod()
    
    anaInfo <- analysisInfo(obj)
    featTab <- as.data.table(getFeatures(obj))
    
    # HACK: since feature tables store the character form, it's easier to keep it all the same
    EIXParams$setsAdductPos <- as.character(EIXParams$setsAdductPos)
    EIXParams$setsAdductNeg <- as.character(EIXParams$setsAdductNeg)
    
    # 'ionize' m/zs
    return(Map(names(ret), ret, f = function(grp, ranges)
    {
        featTabGrp <- featTab[group == grp]
        ranges[, adduct := featTabGrp[match(ranges$analysis, analysis)]$adduct]
        
        if (!EIXParams$onlyPresent && any(is.na(ranges$adduct))) # adduct will be NA for 'missing' features
        {
            # First try to get adduct from other features in the same set: assume that adduct per set for a single
            # feature group is always the same
            adductSets <- unique(featTabGrp[, c("adduct", "set"), with = FALSE])
            ranges[is.na(adduct), adduct := {
                s <- anaInfo$set[match(analysis, anaInfo$analysis)]
                adductSets[set == s]$adduct
            }, by = "analysis"]
            
            # Then fallback to default adducts for a polarity. For this we get the polarity from another feature in the
            # same analysis (if present)
            ranges[is.na(adduct), adduct := sapply(analysis, function(ana)
            {
                t <- featTab[analysis == ana]
                if (nrow(t) == 0)
                    NA # all features were removed
                else if (as.adduct(t$adduct[1])@charge > 0)
                    EIXParams$setsAdductPos
                else
                    EIXParams$setsAdductNeg
            })]
            
            if (any(is.na(ranges$adduct)))
            {
                # If all failed then simply omit
                warning(sprintf("Cannot get adduct information for group '%s' for features in analyses %s", grp,
                                paste0(ranges[is.na(adduct)]$analysis, collapse = ", ")), call. = FALSE)
                ranges <- ranges[!is.na(adduct)]
            }
        }
        
        # NOTE: this is mostly copied from unset features method        
        allAdducts <- sapply(unique(ranges$adduct), as.adduct)
        mzmins <- calculateMasses(ranges$mzmin, allAdducts[ranges$adduct], type = "mz")
        nmd <- mzmins - ranges$mzmin
        set(ranges, j = c("mzmin", "mzmax"), value = list(mzmins, ranges$mzmax + nmd))
        
        return(ranges)
    }))
})

setMethod("getFeatureEIXs", "featuresSet", function(obj, type, ...)
{
    unsetFeatList <- sapply(sets(obj), unset, obj = obj, simplify = FALSE)
    EIXList <- sapply(unsetFeatList, getFeatureEIXs, type = type, ..., simplify = FALSE)
    EIXs <- unlist(EIXList, recursive = FALSE, use.names = FALSE) # use.names gives combined set/ana name, we just want ana
    names(EIXs) <- unlist(lapply(EIXList, names))
    EIXs <- EIXs[intersect(analyses(obj), names(EIXs))] # sync order
    return(EIXs)
})

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
