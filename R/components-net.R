# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
NULL

getNetCompHClust <- function(rmat, method = "complete", h = 0.95)
{
    distm <- as.dist(1 - abs(rmat))
    hc <- fastcluster::hclust(distm, method = method)
    ct <- cutree(hc, h = h)
    clIDs <- sort(unique(ct))
    return(lapply(clIDs, \(id) names(ct[ct == id])))
}

getNetCompHCS <- function(graph, ...)
{
    return(RBGL::highlyConnSG(igraph::as_graphnel(graph), ...)$clusters)
}

getNetCompCommunity <- function(graph, func = igraph::cluster_walktrap, ...)
{
    return(unname(igraph::communities(func(graph, ...))))
}

getNetCompCliques <- function(graph, ...)
{
    cliques <- igraph::max_cliques(graph, ...)
    
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

# NOTE: this function is called by a withProg() block, so handles progression updates
makeCompNetFeatures <- function(fTable, EICs, sim, minSim, maxP, method, ...)
{
    eicm <- do.call(cbind, lapply(EICs, \(eic) eic[, "intensity"]))
    
    if (sim == "pearson")
    {
        corr <- Hmisc::rcorr(eicm, type = "pearson")
        rmat <- corr$r
        rmat[rmat < minSim | corr$P >= maxP] <- 0 # UNDONE: Doc that p threshold is applied.
    }
    else
    {
        rmat <- proxy::simil(eicm, method = "cosine", by_rows = FALSE) |> as.matrix()
        rmat[rmat < minSim] <- 0
    }
    diag(rmat) <- 0
    
    graph <- igraph::graph_from_adjacency_matrix(rmat, mode = "undirected", weighted = TRUE, diag = FALSE)
    compList <- switch(method,
                       community = getNetCompCommunity(graph, ...),
                       cliques = getNetCompCliques(graph, ...),
                       hcs = getNetCompHCS(graph, ...),
                       hclust = getNetCompHClust(rmat, ...))
    
    compTabs <- lapply(compList, function(grps)
    {
        if (length(grps) == 1)
            return(data.table(group = grps, degree = 0, corMin = NA_real_, corMax = NA_real_, corMean = NA_real_))
        
        subg <- igraph::induced_subgraph(graph, grps)
        
        tab <- data.table(group = grps, degree = igraph::degree(subg, normalized = TRUE))
        weights <- sapply(grps, \(g) igraph::E(subg)[igraph::incident(subg, g)]$weight, simplify = FALSE)
        weights <- lapply(weights, \(w) w[w > 0])
        tab[, c("corMin", "corMax", "corMean") := .(sapply(weights, min), sapply(weights, max), sapply(weights, mean))]
        return(tab)
    })
    
    doProgress()
    
    return(list(graph = graph, components = compTabs))
}

annotateCompNetFM <- function(componList, mzWindow, ionization, adducts, ...)
{
    adducts <- sapply(adducts, as.character)
    
    objects <- lapply(componList, function(comp)
    {
        list(fm = InterpretMSSpectrum::findMAIN(comp[, c("mz", "intensity"), with = FALSE], mzabs = mzWindow,
                                                ionmode = ionization, rules = adducts, ...))
    })
    
    componList <- lapply(seq_along(componList), function(i)
    {
        comp <- copy(componList[[i]])
        fmtab <- as.data.table(objects[[i]]$fm[[1]]) # HACK: [[1]] is how the print() method gets the table
        comp[, c("isogroup", "isonr", "charge", "adduct", "ppm") := .(fmtab$isogr, fmtab$iso, fmtab$charge, fmtab$adduct, fmtab$ppm)]
        comp[!is.na(adduct), neutralMass := mapply(mz, adduct, FUN = \(m, a) calculateMasses(m, as.adduct(a), type = "neutral"))]
        return(comp)
    })
    
    return(list(components = componList, objects = objects))
}

annotateCompNetNontarget <- function(componList, mzWindow, adducts, prefAdducts, patArgs, addArgs)
{
    # UNDONE: ignore ret
    
    epEnv <- new.env()
    if (is.null(patArgs[["iso"]])) # UNDONE: doc that this is the default
    {
        data(isotopes, package = "enviPat", envir = epEnv)
        patArgs$iso <- nontarget::make.isos(epEnv$isotopes)
    }
    if (is.null(addArgs[["adducts"]])) # UNDONE: doc that this is the default
    {
        data(adducts, package = "enviPat", envir = epEnv)
        addArgs$adducts <- epEnv$adducts
    }

    adducts <- sapply(adducts, as.character, format = "nontarget", adductInfo = addArgs$adduct, USE.NAMES = FALSE)
    
    indsToGNames <- function(inds, gNames)
    {
        inds <- strsplit(inds, "/")
        return(sapply(inds, \(i) paste0(gNames[as.integer(i)], collapse = "/")))
    }
    
    objects <- lapply(componList, function(comp)
    {
        compS <- comp[, c("mz", "intensity", "ret"), with = FALSE]
        utils::capture.output(ps <- do.call(nontarget::pattern.search, c(list(compS, ppm = FALSE, mztol = mzWindow), patArgs)))
        # NOTE: nontarget::adduct.search() calls stop() when there are no results ...
        utils::capture.output(
            as <- tryCatch(do.call(nontarget::adduct.search, c(list(compS, ppm = FALSE, mztol = mzWindow,
                                                                use_adducts = adducts), addArgs)), error = \(e) e)
        )
        if (inherits(as, "error"))
        {
            if (!grepl("No matches found", as$message))
                stop(as)
            as <- NULL
        }
        return(list(ps = ps, as = as))
    })
    
    componList <- lapply(seq_along(componList), function(i)
    {
        comp <- copy(componList[[i]])
        compS <- comp[, c("mz", "intensity", "ret"), with = FALSE]
        ps <- objects[[i]]$ps
        as <- objects[[i]]$as
        
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
        else
            isoTab[, iso_link := NA_character_]
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
            addObjs <- lapply(addTabLong$adduct, as.adduct, format = "nontarget", adductInfo = addArgs$adducts)
            addTabLong[, adduct := sapply(addObjs, as.character)]
            addTabLong[, adduct_other := sapply(adduct_other, \(ao) as.character(as.adduct(ao, format = "nontarget", adductInfo = addArgs$adducts)))]
            addTabLong[, neutralMass := calculateMasses(comp$mz[match(ID, comp$ID)], addObjs, type = "neutral")]
            
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
        else
        {
            addTab <- data.table(ID = integer(), addgroup = integer(), adduct = character(),
                                 neutralMass = numeric(), add_link = character(),
                                 add_link_adduct = character(), add_link_mz_tol = character())
        }
        
        comp <- merge(comp, isoTab, by = "ID", all.x = TRUE, sort = FALSE)
        comp <- merge(comp, addTab, by = "ID", all.x = TRUE, sort = FALSE)
        comp[, ID := NULL]
        
        return(comp)
    })
    names(componList) <- names(objects)
    
    return(list(components = componList, objects = objects))
}

#' Components class for network-based componentization.
#'
#' This class is derived from \code{\link{components}} and is used to store results from network-based componentization.
#'
#' Objects from this class are generated by \code{\link{generateComponentsNet}}.
#'
#' @slot featureComponents A \code{list} with feature components for each analysis.
#' @slot featureGraphs A \code{list} with feature network graphs (\CRANpkg{igraph}) for each analysis.
#' @slot annotationObjects A \code{list} with annotation objects for each component (from \code{\link{nontarget}} or
#'   \code{\link{InterpretMSSpectrum}}).
#'
#' @seealso \code{\link{components}}, \code{\link{generateComponents}}, \code{\link{generateComponentsNet}}
#'
#' @export
componentsNet <- setClass("componentsNet", slots = c(featureComponents = "list", featureGraphs = "list",
                                                     annotationObjects = "list"),
                          contains = "components")

setMethod("initialize", "componentsNet",
          function(.Object, ...) callNextMethod(.Object, algorithm = "compnet", ...))

#' @rdname components-class
#' @export
setMethod("expandForIMS", "componentsNet", function(obj, ...) cannotExpandComponIMS(obj))

#' Network-based componentization
#'
#' Uses networking to generate components by grouping features based on similarity of their chromatographic elution
#' profiles.
#'
#' @templateVar algo network
#' @templateVar do generate components
#' @templateVar generic generateComponents
#' @templateVar algoParam network
#' @template algo_generator
#'
#' @details Features are first grouped per analysis into feature components by constructing a network where nodes
#'   represent features and edges represent similarity between their elution profiles. Similarity is calculated using
#'   either Pearson correlation or cosine similarity. The network is then clustered using one of several methods (see
#'   \code{componMethod} argument). The resulting feature components are then merged across analyses using a consensus
#'   clustering approach, based on the pairwise presence similarity across analyses. Finally, components are annotated
#'   using either the \pkg{InterpretMSSpectrum} or \pkg{nontarget} package to identify isotopes, in-source fragments and
#'   adducts.
#'
#' @section \code{nontarget} annotation: When \code{annotAlgo="nontarget"}, annotation is performed using
#'   \pkg{nontarget}'s \code{\link[nontarget:pattern.search]{pattern.search}} and
#'   \code{\link[nontarget:adduct.search]{adduct.search}} functions.
#'
#'   \strong{Isotope annotation}: \code{pattern.search} identifies isotope patterns. Peaks are grouped by isotope
#'   patterns and charge levels. When multiple charge levels are found for the same isotope group, the "best" grouping
#'   is selected using the following priority: (1) groups containing a 13C isotope, (2) the largest group,
#'   or (3) the lowest charge level.
#'
#'   \strong{Adduct annotation}: \code{adduct.search} identifies adduct relationships between peaks. Adducts are grouped
#'   by their calculated neutral mass. If there is a conflict in neutral mass assignment, the "best" adduct group is
#'   selected using the following priority: (1) groups containing preferred adducts (see \code{annotPrefAdducts}), (2)
#'   the largest group, or (3) the first group encountered.
#'
#' @param componSim The similarity metric to use: \code{"pearson"} for Pearson correlation or \code{"cosine"} for cosine
#'   similarity.
#' @param componMinSim The minimum similarity threshold (\samp{0-1}). Feature pairs with a lower similarity are not
#'   connected in the network.
#' @param componMaxP The maximum \emph{p}-value threshold (\samp{0-1}) for Pearson correlation. Only applicable when
#'   \code{componSim="pearson"}. Feature pairs with a higher \emph{p}-value are not connected.
#' @param componMethod The network componentization method: \code{"community"} for community detection with
#'   \CRANpkg{igraph}, \code{"cliques"} for maximal cliques (using
#'   \code{\link[igraph:max_cliques]{igraph::max_cliques}}), \code{"hcs"} for highly connected subgraphs (using
#'   \code{\link[RBGL:highlyConnSG]{RBGL::highlyConnSG}}), or \code{"hclust"} for hierarchical clustering (using
#'   \code{\link[fastcluster:hclust]{fastcluster::hclust}}).
#' @param componArgs A \code{list} with additional arguments passed to the network componentization function.
#'
#'   For \code{componMethod="community"}, arguments are passed to the \CRANpkg{igraph} clustering function. The default
#'   function is \code{\link[igraph:cluster_walktrap]{igraph::cluster_walktrap}} and can be changed by setting the
#'   \code{func} argument.
#'
#'   For \code{componMethod="cliques"}: arguments are passed to \code{\link[igraph:max_cliques]{igraph::max_cliques}}.
#'
#'   For \code{componMethod="hcs"}: arguments are passed to \code{\link[RBGL:highlyConnSG]{RBGL::highlyConnSG}}.
#'
#'   For \code{componMethod="hclust"}: arguments \code{method} to set the clustering method (see
#'   \code{\link[fastcluster:hclust]{fastcluster::hclust}}) (\code{"complete" by default}) and \code{h} to set the
#'   clustering height threshold (\samp{0.95} by default).
#' @param groupClust,groupClustH Clustering method (see \code{\link[fastcluster:hclust]{fastcluster::hclust}}) and
#'   height at which to cut the tree for the consensus clustering step across analyses.
#' @param annotAlgo Annotation algorithm: \code{"imss"} for annotation using
#'   \code{\link[InterpretMSSpectrum:findMAIN]{InterpretMSSpectrum::findMAIN}} or \code{"nontarget"} for annotation
#'   using \code{\link[nontarget:pattern.search]{nontarget::pattern.search}} and
#'   \code{\link[nontarget:adduct.search]{nontarget::adduct.search}}.
#' @param annotAdducts A \code{character} vector or \code{adduct} object vector with adducts to use for annotation. At
#'   least two adducts are required. Adducts incompatible with the ionization mode are automatically removed.
#' @param annotPrefAdducts A \code{character} vector with preferential adducts. Only relevant if
#'   \code{annotAlgo="nontarget"}, see the \verb{nontarget annotation} section.
#' @param annotArgs A \code{list} with additional arguments passed to the annotation function. For
#'   \code{annotAlgo="imss"}, arguments are passed to
#'   \code{\link[InterpretMSSpectrum:findMAIN]{InterpretMSSpectrum::findMAIN}}. For \code{annotAlgo="nontarget"}, the
#'   list should contain \code{pattern} and \code{adduct} elements, which are themselves lists with arguments passed to
#'   \code{\link[nontarget:pattern.search]{nontarget::pattern.search}} and
#'   \code{\link[nontarget:adduct.search]{nontarget::adduct.search}}, respectively.
#'
#' @templateVar ion TRUE
#' @templateVar minSize TRUE
#' @template compon_algo-args
#' @template parallel-arg
#'
#' @inheritParams generateComponents
#'
#' @return A \code{\link{componentsNet}} object (derived from \code{\link{components}}). The
#'   \code{\link{componentTable}} of each component contains the following columns:
#'   \describe{
#'   \item{\code{group}}{The feature group name.}
#'   \item{\code{ret}}{The retention time.}
#'   \item{\code{mz}}{The mass-to-charge ratio.}
#'   \item{\code{degreeMin}, \code{degreeMax}, \code{degreeMean}}{The minimum, maximum and mean normalized degree of
#'     the feature in the network across all analyses.}
#'   \item{\code{corMin}, \code{corMax}, \code{corMean}}{The minimum, maximum and mean correlation (or similarity) of
#'     the feature with other features in the component across all analyses.}
#'   \item{\code{intensity}}{The mean intensity of the feature group across analyses where it was detected.}
#'   \item{\code{intensity_rel}}{The relative intensity within the component (\samp{0-1}).}
#'   }
#'
#'   When \code{annotAlgo="imss"}, the following annotation columns are added:
#'   \describe{
#'   \item{\code{isogroup}}{The isotope group.}
#'   \item{\code{isonr}}{The isotope number within the group.}
#'   \item{\code{charge}}{The charge.}
#'   \item{\code{adduct}}{The assigned adduct.}
#'   \item{\code{ppm}}{The mass error in ppm.}
#'   \item{\code{neutralMass}}{The calculated neutral mass.}
#'   }
#'
#'   When \code{annotAlgo="nontarget"}, the output from \pkg{nontarget}'s \code{pattern.search} and \code{adduct.search}
#'   is significantly restructured to simplify interpretation. Isotope and adduct information is
#'   collapsed to one row per feature, the "best" grouping is selected in case of conflicts (see the \strong{nontarget
#'   annotation} section), and links to related features are expressed as feature group names. The following annotation
#'   columns are added:
#'   \describe{
#'   \item{\code{isogroup}}{The isotope group ID.}
#'   \item{\code{iso_interaction}}{The isotope interaction level(s), i.e. the distance from the monoisotope. Multiple
#'     values are slash-separated.}
#'   \item{\code{isotope}}{The isotope label(s) (e.g. \code{"13C"}). \code{"mono"} indicates a monoisotope peak.
#'     Multiple values are slash-separated.}
#'   \item{\code{iso_mz_tol}}{The mass tolerance(s) for the isotope assignment. Multiple values are slash-separated.}
#'   \item{\code{charge}}{The charge level.}
#'   \item{\code{iso_link}}{The feature group name(s) of the origin (monoisotope) peak(s) that this feature is an
#'     isotope of. \code{NA} for monoisotope peaks or features not part of an isotope pattern.}
#'   \item{\code{addgroup}}{The adduct group ID, grouping features that share the same neutral mass.}
#'   \item{\code{adduct}}{The assigned adduct.}
#'   \item{\code{neutralMass}}{The calculated neutral mass.}
#'   \item{\code{add_link}}{The feature group name(s) of the peak(s) linked to this feature via an adduct
#'     relationship.}
#'   \item{\code{add_link_adduct}}{The adduct(s) of the linked peak(s).}
#'   \item{\code{add_link_mz_tol}}{The mass tolerance(s) for the adduct assignment.}
#'   }
#'
#'   The \code{\link{componentInfo}} table contains the columns \code{name} (component name), \code{cmp_ret} (mean
#'   retention time), \code{cmp_retsd} (retention time standard deviation) and \code{size} (number of features).
#'
#' @template compon_ims_unsupported
#'
#' @templateVar class componentsSet
#' @template compon_gen-sets-merged
#'
#' @source The componentization approach was inspired by \pkg{CAMERA} and \pkg{cliqueMS}.
#'
#' @references \addCitations{CAMERA} \cr\cr \insertRef{Senan2019}{patRoon} \cr\cr \insertRef{Loos2017}{patRoon} \cr\cr
#'   \addCitations{enviPat} \cr\cr \insertRef{Jaeger2016}{patRoon} \cr\cr \insertRef{Jaeger2017}{patRoon} \cr\cr
#'   \addCitations{igraph} \cr\cr \addCitations{fastcluster}
#'
#' @seealso \code{\link{componentsNet}}
#'
#' @templateVar what generateComponentsNet
#' @template main-rd-method
#' @export
setMethod("generateComponentsNet", "featureGroups", function(fGroups, ionization = NULL, minSize = 2,
                                                             mzWindow = defaultLim("mz", "medium"),
                                                             componSim = "cosine", componMinSim = 0.95,
                                                             componMaxP = 0.05, componMethod = "community",
                                                             componArgs = list(), groupClust = "complete",
                                                             groupClustH = 0.5, annotAlgo = "imss",
                                                             annotAdducts = c("[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+",
                                                                              "[M-H]-", "[M-H2O-H]-"),
                                                             annotPrefAdducts = c("[M+H]+", "[M-H]-"),
                                                             annotArgs = list())
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    ionization <- checkAndGetIonization(ionization, fGroups, add = ac)
    aapply(checkmate::assertCount, . ~ minSize, positive = TRUE, fixed = list(add = ac))
    checkmate::assertNumber(mzWindow, finite = TRUE, lower = 0, add = ac)
    checkmate::assertChoice(componSim, c("pearson", "cosine"), add = ac)
    checkmate::assertNumber(componMinSim, finite = TRUE, lower = 0, upper = 1, add = ac)
    checkmate::assertNumber(componMaxP, finite = TRUE, lower = 0, upper = 1, add = ac)
    checkmate::assertChoice(componMethod, c("community", "cliques", "hcs", "hclust"), add = ac)
    checkmate::assertList(componArgs, any.missing = FALSE, names = "unique", null.ok = TRUE, add = ac)
    checkmate::assertNumber(groupClustH, finite = TRUE, lower = 0, upper = 1, add = ac)
    checkmate::assertChoice(annotAlgo, c("imss", "nontarget"), add = ac)
    checkmate::assert(
        checkmate::checkCharacter(annotAdducts, min.chars = 1, any.missing = FALSE, min.len = 2),
        checkmate::checkList(annotAdducts, types = "adduct", min.len = 2, any.missing = FALSE),
        .var.name = "annotAdducts", add = ac
    )
    checkmate::assertCharacter(annotAdducts, min.chars = 1, any.missing = FALSE, unique = TRUE, add = ac)
    checkmate::assertCharacter(annotPrefAdducts, min.chars = 1, any.missing = FALSE, unique = TRUE, add = ac)
    checkmate::assertList(annotArgs, any.missing = FALSE, names = "unique", null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    # Check optional dependencies
    if (componSim == "pearson")
        checkPackage("Hmisc")
    if (componSim == "cosine")
        checkPackage("proxy")
    if (componMethod == "hcs")
        checkPackage("RBGL")
    if (annotAlgo == "imss")
        checkPackage("InterpretMSSpectrum")
    if (annotAlgo == "nontarget")
        checkPackage("nontarget")
    
    if (!is.list(annotAdducts))
        annotAdducts <- lapply(annotAdducts, as.adduct)
    
    # only retain relevant adducts for ionization (NOTE: otherwise nontarget throws an error)
    annotAdducts <- annotAdducts[sapply(annotAdducts, \(a) a@charge > 0) == (ionization == "positive")]
    if (length(annotAdducts) < 2)
    {
        stop("Need at least two adducts for annotation, but only ", length(annotAdducts),
             " were provided for ionization '", ionization, "'", call. = FALSE)
    }
    
    hash <- makeHash(fGroups, ionization, minSize, mzWindow, componSim, componMinSim, componMaxP, componMethod,
                     componArgs, groupClust, groupClustH, annotAlgo, annotAdducts, annotPrefAdducts, annotArgs)
    cd <- loadCacheData("componentsNet", hash)
    if (!is.null(cd))
        return(cd)
    
    fTable <- featureTable(fGroups)
    
    printf("Getting EICs and correlation matrices for %d analyses...\n", length(fTable))
    
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

    printf("Generating feature components for %d analyses...\n", length(fTable))
    compsFeats <- withProg(length(fTable), FALSE, Map(fTable, EICs, f = makeCompNetFeatures,
                                                      MoreArgs = c(list(sim = componSim, minSim = componMinSim,
                                                                        maxP = componMaxP, method = componMethod),
                                                                   componArgs)))
    compsFeatsTabs <- sapply(compsFeats, "[[", "components", simplify = FALSE)
    
    # generate consensus components: calculate pairwise grouping of features across analyses
    printf("Generating consensus components across analyses...")
    compsGroupsList <- sapply(compsFeatsTabs, \(x) lapply(x, "[[", "group"), simplify = FALSE)
    featGroupsList <- sapply(fTable, "[[", "group", simplify = FALSE)
    coCount <- getComponNetCoMatrix(compsGroupsList, featGroupsList, names(fGroups))
    
    distm <- as.dist(1 - coCount)
    hc <- fastcluster::hclust(distm, method = groupClust)
    ct <- cutree(hc, h = groupClustH)
    
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
    
    printf("Done! Generated %d components\n", length(componList))
    
    printf("Annotating components using %s... ", annotAlgo)
    annotResult <- if (annotAlgo == "nontarget")
    {
        if (length(annotArgs) == 0)
            annotArgs <- list(pattern = list(), adduct = list())
        annotateCompNetNontarget(componList, mzWindow = mzWindow, adducts = annotAdducts,
                                 prefAdducts = annotPrefAdducts, patArgs = annotArgs$pattern, addArgs = annotArgs$adduct)
    }
    else
        do.call(annotateCompNetFM, c(list(componList, mzWindow = mzWindow, ionization = ionization,
                                          adducts = annotAdducts), annotArgs))
    printf("Done!\n")
    
    componList <- annotResult$components
    
    cInfo <- data.table(name = names(componList), cmp_ret = sapply(componList, function(cmp) mean(cmp$ret)),
                        cmp_retsd = sapply(componList, function(cmp) sd(cmp$ret)),
                        size = sapply(componList, nrow))
    
    ret <- componentsNet(featureComponents = compsFeatsTabs,
                         featureGraphs = sapply(compsFeats, "[[", "graph", simplify = FALSE),
                         annotationObjects = annotResult$objects,
                         componentInfo = cInfo, components = componList)
    saveCacheData("componentsNet", ret, hash)
    return(ret)
})

#' @describeIn componentsNet Plots an interactive network graph for the feature components of an analysis.
#'
#' @param obj The \code{componentsNet} object to plot.
#' @param analysis The name of the analysis to plot.
#'
#' @template plotGraph
#'
#' @references \addCitations{igraph}
#'
#' @export
setMethod("plotGraph", "componentsNet", function(obj, analysis, width = NULL, height = NULL)
{
    checkmate::assertChoice(analysis, names(obj@featureGraphs))
    
    data <- visNetwork::toVisNetworkData(obj@featureGraphs[[analysis]])
    nodes <- as.data.table(data$nodes)
    nodes[, group := sapply(id, \(x) which(sapply(obj@featureComponents[[analysis]], \(y) x %chin% y$group))[1])]
    edges <- data$edges
    edges$value <- edges$weight; edges$title <- round(edges$weight, 2)
    nodes <- nodes[id %in% c(edges$from, edges$to)] # UNDONE: remove singletons during componentization?
    visNetwork::visNetwork(nodes, edges, width = width, height = height) |> visNetwork::visIgraphLayout(physics = TRUE)
})
