#' @include main.R
#' @include components.R
NULL

#' @export
componentsFeatures <- setClass("componentsFeatures", slots = c(featureComponents = "list"),
                               contains = "components")

setMethod("initialize", "componentsFeatures", function(.Object, fGroups, featureComponents = list(), ...)
{
    ftindex <- groupFeatIndex(fGroups)
    gNames <- names(fGroups)
    gInfo <- groupInfo(fGroups)
    
    featureComponents <- Map(featureComponents, split(ftindex, seq_len(nrow(ftindex))), f = function(fCmpL, fti)
    {
        # prune unassigned features
        fCmpL <- Filter(function(cmp) nrow(cmp) > 1 || nzchar(cmp$adduct), fCmpL)
        # convert to data.tables, assign group names prune features without groups
        fti <- unlist(fti)
        fCmpL <- lapply(fCmpL, function(cmp)
        {
            set(setDT(cmp), j = "group", value = gNames[match(cmp$ID, fti)])
            return(cmp[!is.na(group)])
        })
        fCmpL <- pruneList(fCmpL, checkZeroRows = TRUE)
        return(fCmpL)
    })
    
    cmpTab <- rbindlist(lapply(featureComponents, rbindlist, idcol = "fCMP"), idcol = "analysis")
    cmpTab[, fCMPID := paste0(analysis, "-", fCMP)]
    
    # NOTE: this only takes assigned features into account, as unassigned won't be present
    cmpTab[, abundance := sapply(adduct, function(a) sum(a == adduct)) / .N, by = "group"]
    cmpTab <- cmpTab[numGTE(abundance, 0.7)] # UNDONE
    cmpTab[, keep := which.max(abundance), by = "group"]
    cmpTab <- cmpTab[keep == TRUE][, keep := NULL]
    
    usedGroups <- setNames(rep(FALSE, length(fGroups)), gNames)
    linkedFGs <- rbindlist(sapply(gNames, function(gn)
    {
        if (usedGroups[gn])
            return(NULL)
        IDs <- unique(cmpTab[group == gn]$fCMPID)
        ct <- cmpTab[fCMPID %chin% IDs]
        usedGroups[ct$group] <<- TRUE
        return(ct)
    }, simplify = FALSE), idcol = "parentGroup")
    
    # collapse features: only retain one row per feature group
    linkedFGs <- unique(linkedFGs, by = c("parentGroup", "group"))
    
    dups <- function(v) names(which(table(v) > 1))
    
    # remove feature groups that occur in multiple to be components
    linkedFGs <- linkedFGs[!group %chin% dups(group)]
    
    # prepare for components
    linkedFGs <- linkedFGs[, c("parentGroup", "group", "charge", "adduct"), with = FALSE]
    linkedFGs[, c("ret", "mz") := gInfo[group, c("rts", "mzs")]]
    
    comps <- split(linkedFGs, by = "parentGroup", keep.by = FALSE)
    
    # Remove any fGroups from components with equal adducts
    comps <- lapply(comps, function(ct) ct[!adduct %chin% dups(adduct)])
    
    comps <- pruneList(comps, checkZeroRows = TRUE)
    names(comps) <- paste0("CMP", seq_along(comps))
    
    # UNDONE: neutral masses
    cInfo <- data.table(name = names(comps), cmp_ret = sapply(comps, function(cmp) mean(cmp$ret)),
                        cmp_ret = sapply(comps, function(cmp) sd(cmp$ret)),
                        size = sapply(comps, nrow))
    
    return(callNextMethod(.Object, components = comps, componentInfo = cInfo, ...))
})

