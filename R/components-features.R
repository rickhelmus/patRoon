#' @include main.R
#' @include components.R
NULL

#' @export
componentsFeatures <- setClass("componentsFeatures", slots = c(featureComponents = "list"),
                               contains = "components")

setMethod("initialize", "componentsFeatures", function(.Object, fGroups, minSize, mzWindow, relMinAdductAbundance,
                                                       featureComponents = list(), ...)
{
    ftindex <- groupFeatIndex(fGroups)
    gNames <- names(fGroups)
    gInfo <- groupInfo(fGroups)
    
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
    cmpTab[!is.na(adduct_ion), abundance := sapply(adduct_ion, function(a) sum(a == adduct_ion)) / .N, by = "group"]
    
    # Filter adducts not abundantly assigned to same feature group
    cmpTab <- cmpTab[is.na(abundance) | numGTE(abundance, relMinAdductAbundance)]
    
    # Only keep the most abundantly assigned adduct for each feature group
    # UNDONE: handle ties?
    cmpTab[!is.na(adduct_ion), keep := adduct_ion == adduct_ion[which.max(abundance)], by = "group"]
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
        
        if (!allSame(ct$neutralMass, function(x1, x2) numLTE(abs(x1 - x2), mzWindow)))
        {
            # conflict in neutral masses --> only retain those with most abundant neutral mass
            # UNDONE: handle ties?
            
            # since we must work with mass tolerances, first cluster presumably masses together
            hc <- fastcluster::hclust(dist(ct$neutralMass))
            tr <- cutree(hc, h = mzWindow)
            mostAbundantCL <- as.integer(names(which.max(table(tr))))
            ct <- ct[tr == mostAbundantCL]
        }
        
        return(ct)
    }, simplify = FALSE), idcol = "parentGroup")

    # collapse features: only retain one row per feature group
    linkedFGs <- unique(linkedFGs, by = c("parentGroup", "group"))
    
    dups <- function(v) names(which(table(v) > 1))
    
    # remove feature groups that occur in multiple to be components
    linkedFGs <- linkedFGs[!group %chin% dups(group)]
    
    # prepare for components
    cols <- intersect(c("parentGroup", "group", "neutralMass", "isonr", "charge", "adduct_ion"),
                      names(linkedFGs))
    linkedFGs <- linkedFGs[, cols, with = FALSE]
    linkedFGs[, c("ret", "mz") := gInfo[group, c("rts", "mzs")]]
    setcolorder(linkedFGs, c("group", "ret", "mz"))
    
    comps <- split(linkedFGs, by = "parentGroup", keep.by = FALSE)
    
    # Remove any fGroups from components with equal adducts (unless assigned to different isotope)
    if (!is.null(linkedFGs[["isonr"]]))
        comps <- lapply(comps, function(ct) ct[is.na(adduct_ion) | !paste0(adduct_ion, isonr) %chin% dups(paste0(adduct_ion, isonr))])
    else
        comps <- lapply(comps, function(ct) ct[is.na(adduct_ion) | !adduct_ion %chin% dups(adduct_ion)])
    
    # NOTE: minSize should be >= 1 to filter out empty components
    comps <- comps[sapply(comps, nrow) >= minSize]
    
    comps <- calculateComponentIntensities(comps, fGroups)
    
    names(comps) <- paste0("CMP", seq_along(comps))
    
    cInfo <- data.table(name = names(comps), cmp_ret = sapply(comps, function(cmp) mean(cmp$ret)),
                        cmp_retsd = sapply(comps, function(cmp) sd(cmp$ret)),
                        neutral_mass = sapply(comps, function(cmp) mean(cmp$neutralMass)),
                        size = sapply(comps, nrow))
    
    
    return(callNextMethod(.Object, featureComponents = featureComponents, components = comps,
                          componentInfo = cInfo, ...))
})

#' @export
setMethod("show", "componentsFeatures", function(object)
{
    callNextMethod()
    
    printf("Feature components: %d total\n",
           if (length(object@featureComponents) == 0) 0 else sum(lengths(object@featureComponents)))
})
