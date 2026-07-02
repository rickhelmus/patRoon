# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
NULL

getNetCompHClust <- function(rmat)
{
    distm <- as.dist(1 - abs(rmat))
    hc <- fastcluster::hclust(distm, method = "complete") # UNDONE: make method configurable?
    ct <- cutree(hc, h = 0.95) # make configurable/same as applied to rmat
    clIDs <- sort(unique(ct))
    return(lapply(clIDs, \(id) names(ct[ct == id])))
}

getNetCompHCS <- function(graph)
{
    # UNDONE: make configurable
    return(RBGL::highlyConnSG(igraph::as_graphnel(graph), sat = 3, ldv = c(3, 2, 1))$clusters)
}

getNetCompCommunity <- function(graph, func = igraph::cluster_walktrap, ...)
{
    # UNDONE: default OK?
    return(unname(igraph::communities(func(graph, ...))))
}

getNetCompCliques <- function(graph)
{
    # UNDONE: make configurable
    cliques <- igraph::max_cliques(graph)
    
    # ensure cliques don't overlap
    cliques <- cliques[order(sapply(cliques, length), decreasing = TRUE)]
    assigned <- character()
    cliquesF <- list()
    for (clq in cliques)
    {
        take <- setdiff(names(clq), assigned)
        if (length(take) > 0)
        {
            cliquesF[[length(cliquesF) + 1]] <- take
            assigned <- c(assigned, take)
        }
    }
    
    # add singletons as their own cliques (UNDONE?)
    cliquesF <- c(cliquesF, lapply(setdiff(names(igraph::V(graph)), assigned), list))
    
    return(cliquesF)
}

makeCompNetFeatures <- function(fTable, EICs)
{
    # UNDONE: also allow cosine correlation
    # UNDONE: handle new deps
    
    eicm <- do.call(cbind, lapply(EICs, \(eic) eic[, "intensity"]))
    
    corr <- Hmisc::rcorr(eicm, type = "pearson")
    # rmat <- corr$r
    # rmat[rmat < 0.95 | corr$P >= 0.05] <- 0 # UNDONE: make configurable
    # diag(rmat) <- 0

    rmat <- proxy::simil(eicm, method = "cosine", by_rows = FALSE) |> as.matrix()
    rmat[rmat < 0.95] <- 0 # UNDONE: make configurable
    diag(rmat) <- 0
    
    graph <- igraph::graph_from_adjacency_matrix(rmat, mode = "undirected", weighted = TRUE, diag = FALSE)
    compList <- getNetCompCommunity(graph) #getNetCompCliques(graph)  #getNetCompHCS(graph) #getNetCompHClust(rmat)
    
    compTabs <- lapply(compList, function(grps)
    {
        if (length(grps) == 1)
            return(data.table(group = grps, degree = 0, corMin = NA_real_, corMax = NA_real_, corMean = NA_real_))
        
        subg <- igraph::induced_subgraph(graph, grps)
        
        tab <- data.table(group = grps, degree = igraph::degree(subg, normalized = TRUE))
        weights <- sapply(grps, \(g) igraph::E(subg)[igraph::incident(subg, g)]$weight, simplify = FALSE)
        weights <- lapply(weights, function(w) w[w > 0])
        tab[, c("corMin", "corMax", "corMean") := .(sapply(weights, min), sapply(weights, max), sapply(weights, mean))]
        return(tab)
    })
    
    return(list(graph = graph, components = compTabs))
}

annotateCompNetFM <- function(componList, ionization, ...)
{
    # UNDONE: more configuration and defaults
    componList <- lapply(componList, function(comp)
    {
        comp <- copy(comp)
        fm <- InterpretMSSpectrum::findMAIN(comp[, c("mz", "intensity"), with = FALSE], ionmode = ionization, ...)
        fmtab <- as.data.table(fm[[1]]) # HACK: [[1]] is how the print() method gets the table
        comp[, c("isogroup", "isonr", "charge", "adduct", "ppm") := .(fmtab$isogr, fmtab$iso, fmtab$charge, fmtab$adduct, fmtab$ppm)]
        comp[!is.na(adduct), neutalMass := mapply(mz, adduct, FUN = \(m, a) calculateMasses(m, as.adduct(a), type = "neutral"))]
        return(comp)
    })
    
    return(componList)
}

annotateCompNetNontarget <- function(componList, iso, add, ...)
{
    # UNDONE: check deps
    # UNDONE: increase default mass tolerances
    # UNDONE: ignore ret
    
    epEnv <- new.env()
    if (is.null(iso)) # UNDONE: doc that this is the default
    {
        data(isotopes, package = "enviPat", envir = epEnv)
        iso <- nontarget::make.isos(epEnv$isotopes)
    }
    if (is.null(add)) # UNDONE: doc that this is the default
    {
        data(adducts, package = "enviPat", envir = epEnv)
        add <- epEnv$adducts
    }
    
    indsToGNames <- function(inds, gNames)
    {
        inds <- strsplit(inds, "/")
        return(sapply(inds, \(i) paste0(gNames[as.integer(i)], collapse = "/")))
    }
    
    componList <- lapply(componList, function(comp)
    {
        comp <- copy(comp)
        compS <- comp[, c("mz", "intensity", "ret"), with = FALSE]

        # UNDONE: configurable args
        # UNDONE: store objects in slots
        ps <- nontarget::pattern.search(compS, iso = iso, ppm = FALSE, mztol = 0.02)
        # NOTE: nontarget::adduct.search() calls stop() when there are no results ...
        as <- tryCatch(nontarget::adduct.search(compS, adducts = add, ppm = FALSE, mztol = 0.05, use_adducts = c("M+H", "M+K", "M+Na", "M+NH4")), error = \(...) NULL)
        
        comp[, ID := .I]
        
        # Parsing the information from pattern.search() is quite a journey... To summarize:
        # - The annotations are in ps$Patterns
        # - This table is divided into two parts:
        #     1. collapsed isotope grouping and interaction information (=distance from monoisotope) per charge level
        #     2. the collapsed annotations and metadata for each of the peak ID in the "to ID" column
        #        (so _not_ of the peak ID of the row!)
        # - To figure out the actual charge levels of the isotope groups, we need ps[["Peaks in pattern groups"]],
        #   which contains the collapsed charge levels for each group.
        #
        # To make the data a bit easier to parse and make things more consistent, we will make sure that the final table
        # contains the isotope grouping, interaction and annotations all in one row and for the peak ID of that row.
        # Thus, the isotope grouping+interatcion is repeated for each annotation, i.e. like other metadata. Furthermore,
        # instead of pointing to a peak with a higher interaction level, we point to the origin peak. Also, we add
        # "mono" to the peaks that are monoisotopes. Finally, we select the "best" isotope grouping in case there are
        # multiple charge levels.
        #
        # To get there, all the collapsed information is first converted to long format, then merged and finally
        # collapsed again to one row per peak ID.
        
        rmCols <- c(names(compS), "int", "m/z")
        
        isoTab <- as.data.table(ps$Patterns[, setdiff(names(ps$Patterns), rmCols)])
        setnames(isoTab,
                 c("peak ID", "group ID", "interaction level", "to ID", "isotope(s)", "mass tolerance", "charge level"),
                 c("ID", "isogroup", "iso_interaction", "iso_to", "isotope", "iso_mz_tol", "charge"))
        isoTab <- isoTab[isogroup != 0]
        
        isoPeaks <- isoTab[, c("ID", "isogroup", "iso_interaction"), with = FALSE]
        if (nrow(isoPeaks) > 0)
        {
            isoPeaks <- rbindlist(lapply(seq_len(nrow(isoPeaks)), function(row)
            {
                data.table(ID = isoPeaks$ID[row],
                           isogroup = as.integer(strsplit(isoPeaks$isogroup[row], "/")[[1]]),
                           iso_interaction = as.integer(strsplit(isoPeaks$iso_interaction[row], "/")[[1]]))
            }))
            isoGroupsCharges <- rbindlist(lapply(seq_len(nrow(ps[["Peaks in pattern groups"]])), function(row)
            {
                # NOTE: charge is only a character if multiple were collapsed
                data.table(groupID = as.integer(strsplit(ps[["Peaks in pattern groups"]][["group ID"]][row], "/")[[1]][-1]),
                           charge = as.integer(strsplit(as.character(ps[["Peaks in pattern groups"]][["charge level"]][row]), "/")[[1]]))
            }), idcol = "isoCluster")
            
            isoPeaks[isoGroupsCharges, c("charge", "isoCluster") := .(i.charge, i.isoCluster), on = c("isogroup" = "groupID")]
            
            isoCands <- isoTab[, c("ID", "iso_to", "isotope", "iso_mz_tol", "charge"), with = FALSE][!is.na(isotope) & iso_to != "0"]
            isoCands <- rbindlist(lapply(seq_len(nrow(isoCands)), function(row)
            {
                data.table(ID = as.integer(strsplit(isoCands$iso_to[row], "/")[[1]]),
                           isotope = strsplit(isoCands$isotope[row], "/")[[1]],
                           iso_mz_tol = strsplit(isoCands$iso_mz_tol[row], "/")[[1]],
                           charge = as.integer(strsplit(isoCands$charge[row], "/")[[1]]),
                           iso_link = isoCands$ID[row])
            }))
            
            isoCands <- merge(isoPeaks, isoCands, by = c("ID", "charge"), all.x = TRUE, sort = FALSE)
            isoCands[is.na(isotope), isotope := "mono"]
            setorderv(isoCands, c("ID", "isogroup"))
            
            isoGroups <- isoCands[, .(has13C = any(isotope == "13C"), size = uniqueN(ID),
                                      isoCluster = unique(isoCluster), charge = unique(charge)), by = isogroup]
            
            # keep if
            # cluster size == 1 OR
            # cluster size > 1 AND (has 13C AND no other cluster has 13C) OR
            # size is largest
            
            isoGroups[, keep := {
                wh13C <- which(has13C); whSzMax <- which.max(size); whChMin <- which.min(charge)
                if (.N == 1)
                    TRUE
                else if (length(wh13C) == 1)
                    seq_len(.N) == wh13C
                else if (length(whSzMax) == 1)
                    seq_len(.N) == whSzMax
                else
                    seq_len(.N) == whChMin
            }, by = "isoCluster"]
            
            isoCands <- isoCands[isogroup %in% isoGroups[keep == TRUE]$isogroup]
            # NOTE: isogroup and charge should now be a single value for each cluster due to above filtering
            
            isoCands <- isoCands[, .(isogroup = as.integer(isogroup)[1],
                                     iso_interaction = paste0(iso_interaction, collapse = "/"),
                                     isotope = paste0(isotope, collapse = "/"),
                                     iso_mz_tol = paste0(iso_mz_tol, collapse = "/"),
                                     charge = as.integer(charge)[1],
                                     iso_link = paste0(comp[match(iso_link, ID)]$group, collapse = "/")), by = ID]

            isoTab[isoCands, c("isogroup", "iso_interaction", "isotope", "iso_mz_tol", "charge", "iso_link") :=
                       .(i.isogroup, i.iso_interaction, i.isotope, i.iso_mz_tol, i.charge, iso_link), on = "ID"]
        }
        isoTab[, iso_to := NULL]
        
        addTab <- NULL
        if (!is.null(as))
        {
            # Similarly to isotope information, we apply some data transformation to make things easier to parse. The
            # adduct from/to information is split over two columns. Then a grouping is made based on neutral mass, and
            # the "best" adduct group is selected in case of conflicts.

            addTab <- as.data.table(as$adducts[, setdiff(names(as$adducts), rmCols)])
            setnames(addTab, c("peak ID", "group ID", "to ID", "adduct(s)", "mass tolerance"),
                     c("ID", "addgroup", "add_to", "adduct", "add_mz_tol"))
            addTab <- addTab[addgroup != 0]
            addTabLong <- rbindlist(lapply(seq_len(nrow(addTab)), function(row)
            {
                ret <- data.table(ID = addTab$ID[row],
                                  adduct = strsplit(addTab$adduct[row], "//")[[1]],
                                  add_mz_tol = strsplit(addTab$add_mz_tol[row], "/")[[1]],
                                  add_to = strsplit(addTab$add_to[row], "/")[[1]])
                ret[, addgroup := addTab$addgroup[match(ID, addTab$ID)]]
                ret[, adduct_other := sub(".*>", "", adduct)]
                ret[, adduct := sub("<.*", "", adduct)]
                return(ret)
            }))
            
            # add adduct groups per neutral mass: for each ID, assign unique IDs per adduct and assign the same ID to
            # the IDs of corresponding add_to/adduct_to pairs.
            
            addTabLong[, addgroup2 := .GRP, by = .(ID, adduct)]
            for (row in seq_len(nrow(addTabLong)))
            {
                if (addTabLong$ID[row] < addTabLong$add_to[row])
                {
                    wh <- which(addTabLong$ID == addTabLong$add_to[row] &
                                    addTabLong$adduct == addTabLong$adduct_other[row])
                    set(addTabLong, i = wh, j = "addgroup2", value = addTabLong$addgroup2[row])
                }
            }
            
            # convert adduct and calculate neutral masses
            addObjs <- lapply(addTabLong$adduct, as.adduct, format = "nontarget", adductInfo = add)
            addTabLong[, adduct := sapply(addObjs, as.character)]
            addTabLong[, adduct_other := sapply(adduct_other, \(ao) as.character(as.adduct(ao, format = "nontarget", adductInfo = add)))]
            addTabLong[, neutralMass := calculateMasses(comp$mz[match(ID, comp$ID)], addObjs, type = "neutral")]
            
            prefAdducts <- c("[M+H]+", "[M-H]-") # UNDONE: make configurable
            addGroups <- addTabLong[, .(size = .N,
                                        prefMatch = min(match(adduct, prefAdducts, nomatch = length(prefAdducts) + 1))),
                                    by = "addgroup2"]
            # select 'best' adduct group in case there are neutral mass conflicts.
            addTabLong[, sel := {
                ag <- addgroup2
                grps <- copy(addGroups[addgroup2 %in% ag])
                if (nrow(grps) == 1)
                    TRUE
                else
                {
                    grps[, keep := {
                        if (any(prefMatch <= length(prefAdducts)))
                            seq_len(.N) == which.min(prefMatch)
                        else if (!allSame(size))
                            seq_len(.N) == which(size == max(size))
                        else
                            seq_len(.N) == which(addgroup2 == grps[1])
                    }]
                    addgroup2 %in% grps[keep == TRUE]$addgroup2
                }
            }, by = "ID"]
            
            # make sure both in an adduct pair are (de)selected
            addTabLong[sel == TRUE, sel := {
                addTabLong[adduct_other == .SD$adduct & ID == .SD$add_to & adduct == .SD$adduct_other & addgroup2 == .SD$addgroup2]$sel
            }, .SDcols = c("adduct", "add_to", "adduct_other", "addgroup2"), by = .I]
            addTabLong <- addTabLong[sel == TRUE]
            setorderv(addTabLong, c("ID", "addgroup2", "add_to"))

            addTab <- addTabLong[, .(addgroup = as.integer(unique(addgroup2)),
                                     adduct = unique(adduct),
                                     neutralMass = unique(neutralMass),
                                     add_link = paste0(comp[match(add_to, ID)]$group, collapse = "/"),
                                     add_link_adduct = paste0(adduct_other, collapse = "/"),
                                     add_link_mz_tol = paste0(add_mz_tol, collapse = "/")), by = ID]
        }
        
        comp <- merge(comp, isoTab, by = "ID", all.x = TRUE, sort = FALSE)
        if (!is.null(addTab))
            comp <- merge(comp, addTab, by = "ID", all.x = TRUE, sort = FALSE)
        comp[, ID := NULL]
        
        return(comp)
    })
    
    return(componList)
}

#' @rdname components-class
#' @export
componentsNet <- setClass("componentsNet", slots = c(featureComponents = "list", featureGraphs = "list"),
                          contains = "components")

setMethod("initialize", "componentsNet",
          function(.Object, ...) callNextMethod(.Object, algorithm = "compnet", ...))

#' @rdname components-class
#' @export
setMethod("expandForIMS", "componentsNet", function(obj, ...) cannotExpandComponIMS(obj))

setMethod("generateComponentsNet", "featureGroups", function(fGroups, ionization = NULL, minSize = 2)
{
    # checkPackage("cliqueMS", "rickhelmus/cliqueMS") # UNDONE
    
    # UNDONE: asserts
    checkmate::assertCount(minSize, positive = TRUE)
    
    fTable <- featureTable(fGroups)
    
    # EICs: get complete chromatograms so these can be compared, however, only keep the feature signal so other peaks
    # will not interfere correlation calculations.
    EICs <- getFeatureEIXs(fGroups, "EIC", EIXParams = getDefEICParams(window = Inf))
    EICs <- Map(names(EICs), EICs, f = function(ana, anaEICs)
    {
        at <- attr(anaEICs, "allXValues")
        return(Map(names(anaEICs), anaEICs, f = function(fg, eic)
        {
            m <- cbind(time = at, intensity = doFillEIXIntensities(at, eic[, "time"], eic[, "intensity"]))
            ft <- fTable[[ana]][group == fg]
            # UNDONE: limit retmin/retmax?
            m[m[, "time"] < ft$retmin | m[, "time"] > ft$retmax, "intensity"] <- 0
            # m[, "intensity"] <- m[, "intensity"] / max(m[, "intensity"]) * 100 # UNDONE: need normalization?
            return(m)
        }))
    })

    compsFeats <- Map(fTable, EICs, f = makeCompNetFeatures)
    compsFeatsTabs <- sapply(compsFeats, "[[", "components", simplify = FALSE)
    
    # generate consensus components: calculate pairwise grouping of features across analyses
    compsGroupsList <- sapply(compsFeatsTabs, \(x) lapply(x, "[[", "group"), simplify = FALSE)
    featGroupsList <- sapply(fTable, "[[", "group", simplify = FALSE)
    coCount <- getComponNetCoMatrix(compsGroupsList, featGroupsList, names(fGroups))
    
    distm <- as.dist(1 - coCount)
    hc <- fastcluster::hclust(distm, method = "complete") # UNDONE: make configurable
    ct <- cutree(hc, h = 0.5) # UNDONE: make configurable
    
    gInfo <- groupInfo(fGroups)
    allFeatCompTab <- rbindlist(lapply(compsFeatsTabs, \(x) rbindlist(x, idcol = "compID")), idcol = "analysis")
    componList <- lapply(sort(unique(ct)), function(id)
    {
        tab <- data.table(group = names(ct)[ct == id])
        tab[, c("ret", "mz") := .(gInfo$ret[match(group, gInfo$group)], gInfo$mz[match(group, gInfo$group)])]
        tab[, c("degreeMin", "degreeMax", "degreeMean", "corMin", "corMax", "corMean") := {
            fct <- allFeatCompTab[group %chin% grp, env = I(list(grp = group))]
            cors <- if (all(is.na(fct$corMin)))
                .(NA_real_, NA_real_, NA_real_)
            else
                .(min(fct$corMin, na.rm = TRUE), max(fct$corMax, na.rm = TRUE), mean(fct$corMean, na.rm = TRUE))
            c(.(min(fct$degree), max(fct$degree), mean(fct$degree)), cors)
        }, by = "group"]
        
        return(tab)
        # UNDONE: add more metadata?
    })
    
    if (length(componList) > 0)
    {
        # UNDONE: this will take mean feature intensities of _all_ features in the group, including those not in any
        # feature component.
        componList <- calculateComponentIntensities(componList, fGroups)
        names(componList) <- paste0("CMP", seq_along(componList))
    }
    
    # UNDONE: also filter feature components by size? Then also need to update graphs for plotting
    componList <- componList[sapply(componList, nrow) >= minSize]
    
    # UNDONE
    # componList <- annotateCompNetNontarget(componList, iso = NULL, add = NULL)
    componList <- annotateCompNetFM(componList, ionization = ionization)
    
    cInfo <- data.table(name = names(componList), cmp_ret = sapply(componList, function(cmp) mean(cmp$ret)),
                        cmp_retsd = sapply(componList, function(cmp) sd(cmp$ret)),
                        # neutral_mass = sapply(componList, function(cmp) mean(cmp$neutralMass)),
                        size = sapply(componList, nrow))
    
    return(componentsNet(featureComponents = compsFeatsTabs,
                         featureGraphs = sapply(compsFeats, "[[", "graph", simplify = FALSE),
                         componentInfo = cInfo, components = componList))
})

setMethod("plotGraph", "componentsNet", function(obj, analysis)
{
    checkmate::assertChoice(analysis, names(obj@featureGraphs))
    
    data <- visNetwork::toVisNetworkData(obj@featureGraphs[[analysis]])
    nodes <- as.data.table(data$nodes)
    nodes[, group := sapply(id, \(x) which(sapply(obj@featureComponents[[analysis]], \(y) x %chin% y$group))[1])]
    edges <- data$edges
    edges$value <- edges$weight; edges$title <- round(edges$weight, 2)
    nodes <- nodes[id %in% c(edges$from, edges$to)] # UNDONE: remove singletons during componentization
    visNetwork::visNetwork(nodes, edges) |> visNetwork::visIgraphLayout(physics = TRUE)
})
