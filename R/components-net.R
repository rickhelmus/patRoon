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

annotateCompNetNontarget <- function(componList, iso, add, ...)
{
    # UNDONE: check deps
    
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
    
    componList <- lapply(componList, function(comp)
    {
        comp <- copy(comp)

        # UNDONE: configurable args        
        ps <- nontarget::pattern.search(comp[, .(mz, intensity, ret)], iso = iso)
        # NOTE: nontarget::adduct.search() calls stop() when there are no results ...
        as <- tryCatch(nontarget::adduct.search(comp[, .(mz, intensity, ret)], adducts = add),
                       error = function(...) NULL)
        cb <- nontarget::combine(ps, as)
        browser()
        # comp[, c("isogroup", "isonr", "charge") := {
        #     ps <- nontarget::pattern.search(featureTable(fGroups)[[1]][group %in% group]$mz, iso = iso)
        #     if (nrow(ps) == 0)
        #         .(NA_real_, NA_real_, NA_real_)
        #     else
        #         .(ps$isogroup, ps$isonr, ps$charge)
        # }, by = "group"]
        
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
    
    gCount <- length(fGroups)
    gNames <- names(fGroups)
    coCount <- matrix(0, nrow = gCount, ncol = gCount, dimnames = list(gNames, gNames))
    coPossible <- matrix(0, nrow = gCount, ncol = gCount, dimnames = list(gNames, gNames))
    
    for (ana in analyses(fGroups))
    {
        for (cmp in compsFeats[[ana]]$components)
        {
            for (fgi in cmp$group)
            {
                for (fgj in cmp$group)
                {
                    coCount[fgi, fgj] <- coCount[fgi, fgj] + 1
                    if (fgi != fgj)
                        coCount[fgj, fgi] <- coCount[fgj, fgi] + 1
                }
            }
        }
        
        ft <- fTable[[ana]]
        for (fgi in ft$group)
        {
            for (fgj in ft$group)
            {
                coPossible[fgi, fgj] <- coPossible[fgi, fgj] + 1
                if (fgi != fgj)
                    coPossible[fgj, fgi] <- coPossible[fgj, fgi] + 1
            }
        }
    }
    
    coFrac <- coCount / coPossible
    distm <- as.dist(1 - coFrac)
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
    componList <- annotateCompNetNontarget(componList, iso = NULL, add = NULL)
    
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
