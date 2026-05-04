# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
NULL

getNetCompHClust <- function(rmat)
{
    distm <- as.dist(1 - abs(rmat))
    hc <- hclust(distm, method = "complete") # UNDONE: make method configurable?
    ct <- cutree(hc, h = 0.9)
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
    
    rmat <- corr$r
    rmat[rmat < 0.9 | corr$P >= 0.05] <- 0 # UNDONE: make configurable
    diag(rmat) <- 0
    
    graph <- igraph::graph_from_adjacency_matrix(rmat, mode = "undirected", weighted = TRUE, diag = FALSE)
    comps <- getNetCompCliques(graph) #getNetCompCommunity(graph) #getNetCompHCS(graph) #getNetCompHClust(rmat)
    
    return(list(graph = graph, components = comps))
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

setMethod("generateComponentsNet", "featureGroups", function(fGroups, ionization = NULL)
{
    # checkPackage("cliqueMS", "rickhelmus/cliqueMS") # UNDONE
    
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
    return(componentsNet(featureComponents = sapply(compsFeats, "[[", "components", simplify = FALSE),
                         featureGraphs = sapply(compsFeats, "[[", "graph", simplify = FALSE),
                         componentInfo = data.table(), components = list()))
})

setMethod("plotGraph", "componentsNet", function(obj, analysis)
{
    checkmate::assertChoice(analysis, names(obj@featureGraphs))
    
    data <- visNetwork::toVisNetworkData(obj@featureGraphs[[analysis]])
    nodes <- as.data.table(data$nodes)
    nodes[, group := sapply(id, \(x) which(sapply(obj@featureComponents[[analysis]], \(y) x %in% y))[1])]
    edges <- data$edges
    edges$value <- edges$weight; edges$title <- round(edges$weight, 2)
    nodes <- nodes[id %in% c(edges$from, edges$to)] # UNDONE: remove singletons during componentization
    visNetwork::visNetwork(nodes, edges) |> visNetwork::visIgraphLayout(physics = TRUE)
})
