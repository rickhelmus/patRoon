# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include components.R
NULL

#' @rdname components-class
#' @export
componentsFeatures <- setClass("componentsFeatures", slots = c(featureComponents = "list"),
                               contains = "components")

setMethod("initialize", "componentsFeatures", function(.Object, fGroups, minSize, absMzDev, relMinAdductAbundance,
                                                       adductConflictsUsePref, NMConflicts, prefAdducts,
                                                       featureComponents = list(), ...)
{
    if (length(fGroups) == 0)
        return(callNextMethod(.Object, featureComponents = data.table(), components = list(),
                              componentInfo = data.table(), ...))
    
    ftindex <- groupFeatIndex(fGroups)
    gNames <- names(fGroups)
    gInfo <- groupInfo(fGroups)
    anas <- analyses(fGroups)
    gTable <- groupTable(fGroups)
    
    featureComponents <- Map(featureComponents, split(ftindex, seq_len(nrow(ftindex))), f = function(fCmpL, fti)
    {
        # assign group names and prune features without groups
        fti <- unlist(fti)
        fCmpL <- lapply(fCmpL, function(cmp)
        {
            set(cmp, j = "group", value = gNames[match(cmp$ID, fti)])
            return(cmp[!is.na(group)])
        })
        fCmpL <- pruneList(fCmpL, checkZeroRows = TRUE)
        return(fCmpL)
    })
    
    # fCMP: unique feature component ID within an analysis
    cmpTab <- rbindlist(lapply(featureComponents, rbindlist, idcol = "fCMP"), idcol = "analysis")
    
    if (nrow(cmpTab) == 0)
        return(callNextMethod(.Object, featureComponents = featureComponents, components = list(),
                              componentInfo = data.table(), ...))
    
    # fCMPID: unique identifier throughout all analyses
    cmpTab[, fCMPID := paste0(match(analysis, analyses(fGroups)), "-", fCMP)]

    # NOTE: abundance only takes assigned features into account, as unassigned won't be present
    cmpTab[!is.na(adduct_ion), adduct_abundance := sapply(adduct_ion, function(a) sum(a == adduct_ion)) / .N, by = "group"]
    
    # Filter adducts not abundantly assigned to same feature group
    cmpTab <- cmpTab[is.na(adduct_abundance) | numGTE(adduct_abundance, relMinAdductAbundance)]
    
    if (adductConflictsUsePref && length(prefAdducts) > 0)
    {
        # for fGroups with features that have a preferential adduct: remove all others or ones that are lower ranked
        cmpTab[!is.na(adduct_ion), prefInd := match(adduct_ion, prefAdducts, nomatch = length(prefAdducts) + 1),
               by = "group"]
        # NOTE: below leaves features untouched if none of the adducts are preferential, since prefInd will be the same
        # for all and thus all are equal to min(prefInd)
        cmpTab[, keep := is.na(adduct_ion) | prefInd == min(prefInd), by = "group"]
        cmpTab <- cmpTab[is.na(adduct_ion) | keep == TRUE][, keep := NULL]
    }
    
    # Only the most abundantly assigned adduct for each feature group. NOTE: if preferential adducts were selected above
    # then these are now always the most abundant.
    # UNDONE: handle ties?
    cmpTab[!is.na(adduct_ion), keep :=
               uniqueN(adduct_ion) == 1 | adduct_ion == adduct_ion[which.max(adduct_abundance)], by = "group"]
    cmpTab <- cmpTab[is.na(adduct_ion) | keep == TRUE][, keep := NULL]

    # Start making group components; for each feature group:
    # - find all feature components that this "parent group" is in
    # - assume that this group and all groups in the feature components are related
    # - mark all these groups so they won't be used as parent groups next iterations
    # NOTE: we can do the latter since groups present in multiple components will be removed afterwards, hence nothing
    # will be 'missed'.
    usedGroups <- setNames(rep(FALSE, length(gNames)), gNames)
    linkedFGs <- rbindlist(sapply(gNames, function(gn)
    {
        if (usedGroups[gn])
            return(NULL)
        IDs <- unique(cmpTab[group == gn]$fCMPID)
        ct <- cmpTab[fCMPID %chin% IDs]
        
        # special case: if parent has NA neutralMass, ensure all others also have and vice versa. This ensures a split
        # between those with NAs and those without. If the parentGroup has mixed NAs, the ones with NA will be ignored.
        ct <- if (any(!is.na(ct[group == gn]$neutralMass))) ct[!is.na(neutralMass)] else ct[is.na(neutralMass)]
        
        usedGroups[ct$group] <<- TRUE
        
        if (!allSame(ct$neutralMass, function(x1, x2) numLTE(abs(x1 - x2), absMzDev)))
        {
            # conflict in neutral masses --> only retain those with most abundant neutral mass
            # UNDONE: handle ties?

            nrowPrior <- nrow(ct)
            colsPrior <- copy(names(ct)) # NOTE: need a copy: https://stackoverflow.com/a/15913648
            
            # since we must work with mass tolerances, first cluster presumably masses together
            hc <- fastcluster::hclust(dist(ct$neutralMass))
            ct[, clust := cutree(hc, h = absMzDev)]
            ct[, clust_size := .N, by = "clust"]
            
            for (cnf in NMConflicts)
            {
                bestCL <- integer()
                if (cnf == "preferential" && length(prefAdducts) > 0)
                {
                    ct[, prefAdductMatch := match(adduct_ion, prefAdducts, nomatch = length(prefAdducts) + 1),
                       by = "clust"]
                    ct[, topRankedMatch := min(prefAdductMatch), by = "clust"]
                    
                    if (uniqueN(ct$topRankedMatch) == max(ct$clust)) # all clusters ranked differently?
                        bestCL <- ct[which.min(topRankedMatch)]$clust
                }
                else if (cnf == "mostAbundant")
                    bestCL <- ct[which.max(clust_size)]$clust
                else # mostIntense
                {
                    ct[, intensity := mapply(analysis, group, FUN = function(a, g) gTable[[g]][match(a, anas)])]
                    ct[, maxClustInt := max(intensity) / clust_size, by = "clust"]
                    bestCL <- ct[which.max(intensity)]$clust
                }
                
                if (length(bestCL) == 1)
                {
                    ct <- ct[clust == bestCL]
                    break
                }
            }
            if (nrow(ct) == nrowPrior)
            {
                # Could not resolve neutral mass conflict --> just default to first cluster...
                ct <- ct[clust == 1]
            }
            
            ct <- ct[, colsPrior, with = FALSE] # remove temporary work columns
        }
        
        return(ct)
    }, simplify = FALSE), idcol = "parentGroup")

    # collapse features: only retain one row per feature group
    linkedFGs <- unique(linkedFGs, by = c("parentGroup", "group"))
    
    # remove feature groups that occur in multiple to be components
    linkedFGs <- linkedFGs[!group %chin% getDuplicatedStrings(group)]
    
    # prepare for components
    cols <- intersect(c("parentGroup", "group", "neutralMass", "isonr", "charge", "adduct_ion", "adduct_abundance"),
                      names(linkedFGs))
    linkedFGs <- linkedFGs[, cols, with = FALSE]
    linkedFGs[, c("ret", "mz") := gInfo[group, c("rts", "mzs")]]
    setcolorder(linkedFGs, c("group", "ret", "mz"))
    
    comps <- split(linkedFGs, by = "parentGroup", keep.by = FALSE)
    
    # Remove any fGroups from components with equal adducts (unless assigned to different isotope)
    if (!is.null(linkedFGs[["isonr"]]))
        comps <- lapply(comps, function(ct) ct[is.na(adduct_ion) | !paste0(adduct_ion, isonr) %chin% getDuplicatedStrings(paste0(adduct_ion, isonr))])
    else
        comps <- lapply(comps, function(ct) ct[is.na(adduct_ion) | !adduct_ion %chin% getDuplicatedStrings(adduct_ion)])
    
    # NOTE: minSize should be >= 1 to filter out empty components
    comps <- comps[sapply(comps, nrow) >= minSize]
    
    if (length(comps) > 0)
    {
        comps <- calculateComponentIntensities(comps, fGroups)
        names(comps) <- paste0("CMP", seq_along(comps))
    }
    
    cInfo <- data.table(name = names(comps), cmp_ret = sapply(comps, function(cmp) mean(cmp$ret)),
                        cmp_retsd = sapply(comps, function(cmp) sd(cmp$ret)),
                        neutral_mass = sapply(comps, function(cmp) mean(cmp$neutralMass)),
                        size = sapply(comps, nrow))
    
    
    return(callNextMethod(.Object, featureComponents = featureComponents, components = comps,
                          componentInfo = cInfo, ...))
})

#' @rdname components-class
#' @export
setMethod("show", "componentsFeatures", function(object)
{
    callNextMethod()
    
    printf("Feature components: %d total\n",
           if (length(object@featureComponents) == 0) 0 else sum(lengths(object@featureComponents)))
})
