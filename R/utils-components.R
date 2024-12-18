# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
NULL

# from https://stackoverflow.com/a/17531678
orderComponentsNames <- function(n) order(nchar(n), n)

mergeComponents <- function(compList, compNames, nameColumn)
{
    notEmpty <- lengths(compList) > 0
    compList <- compList[notEmpty]; compNames <- compNames[notEmpty]
    
    if (length(compList) == 0)
        return(list(components = list(), componentInfo = data.table()))
    
    retCInfo <- copy(componentInfo(compList[[1]]))
    retCInfo[, (nameColumn) := compNames[[1]]]
    retCTable <- componentTable(compList[[1]])

    if (length(compList) > 1)
    {
        for (mi in seq(2, length(compList)))
        {
            if (length(compList[[mi]]) == 0)
                next
            rci <- copy(componentInfo(compList[[mi]]))
            rci[, (nameColumn) := compNames[[mi]]]
            retCInfo <- rbind(retCInfo, rci, fill = TRUE)
            retCTable <- c(retCTable, componentTable(compList[[mi]]))
        }
    }               
    
    retCInfo[, name := paste0(name, "-", get(nameColumn))][]
    names(retCTable) <- retCInfo[["name"]]
    
    return(list(components = retCTable, componentInfo = retCInfo))
}

calculateComponentIntensities <- function(comps, fGroups)
{
    getGroupInt <- function(grp)
    {
        ints <- fGroups[[grp]]
        return(mean(ints[ints != 0]))
    }
    return(lapply(comps, function(cmp)
    {
        cmp <- copy(cmp)
        cmp[, intensity := sapply(group, getGroupInt)]
        cmp[, intensity_rel := intensity / max(intensity)]
        return(cmp[])
    }))
}

makeGraph <- function(components, onlyLinked, titles, width, height)
{
    cInfo <- copy(components@componentInfo)
    
    # convert link IDs to numeric indices
    cInfo[, id := .I]
    cInfo[, linksIDs := lapply(links, match, table = names(components))]
    allLinks <- unique(unlist(cInfo$linksIDs))
    cInfo <- cInfo[lengths(links) > 0 | id %in% allLinks]
    
    if (nrow(cInfo) == 0)
    {
        nodes <- data.table(id = character(), label = character(), group = numeric())
        edges <- data.table(from = character(), to = character())
    }
    else
    {
        edges <- rbindlist(mapply(cInfo$name, cInfo$linksIDs, FUN = function(n, l)
        {
            data.table::data.table(from = n, to = cInfo$name[match(unlist(l), cInfo$id)])
        }, SIMPLIFY = FALSE))
        
        graph <- igraph::simplify(igraph::graph_from_data_frame(edges, directed = FALSE))
        fc <- igraph::fastgreedy.community(graph)
        
        data <- visNetwork::toVisNetworkData(graph)
        nodes <- as.data.table(data$nodes)
        nodes[, group := fc$membership]
        edges <- data$edges
    }

    if (!onlyLinked)
    {
        unNodes <- data.table(id = setdiff(names(components), cInfo$name), group = 0)
        unNodes[, label := id]
        nodes <- rbind(nodes, unNodes)
    }
    
    nodes[, shape := "circle"]
    nodes[, title := titles[match(id, names(components))]]
    
    visNetwork::visNetwork(nodes = nodes, edges = edges, width = width, height = height)
}

printComponentsFiltered <- function(old, new)
{
    oldn <- length(old); newn <- length(new)
    oldresn <- if (oldn > 0) sum(sapply(old@components, nrow)) else 0
    newresn <- if (newn > 0) sum(sapply(new@components, nrow)) else 0
    printf("Done! Filtered %d (%.2f%%) components and %d (%.2f%%) entries. Remaining: %d components with %d entries.\n",
           oldn - newn, if (oldn == 0) 0 else (1 - (newn / oldn)) * 100,
           oldresn - newresn, if (oldresn == 0) 0 else (1 - (newresn / oldresn)) * 100, newn, newresn)
}

getIonizationFromAnnTable <- function(annTable)
{
    # UNDONE: default to positive OK?
    
    if (nrow(annTable) == 0)
        return("positive")
    
    allAdducts <- lapply(unique(annTable$adduct), as.adduct)
    allCharges <- sapply(allAdducts, slot, "charge")
    chRange <- range(allCharges)
    if (all(chRange < 0))
        return("negative")
    else if (all(chRange >= 0))
        return("positive")
    
    warning("Cannot determine ionization: both positive and negative adducts annotations are present. Defaulting to positive.",
            call. = FALSE)
    
    return("positive")
}
