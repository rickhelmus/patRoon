# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

getTPParents <- function(parents, skipInvalid, prefCalcChemProps, neutralChemProps, checkWhat = "SMILES")
{
    if ((is.data.frame(parents) && nrow(parents) == 0) || length(parents) == 0)
        return(data.table(name = character(), SMILES = character(), InChI = character(), InChIKey = character(),
                          formula = character(), neutralMass = numeric()))
    
    if (is.data.frame(parents))
        parents <- prepareSuspectList(parents, NULL, skipInvalid, checkDesc = TRUE,
                                      prefCalcChemProps = prefCalcChemProps, neutralChemProps = neutralChemProps,
                                      calcMZs = FALSE)
    else if (inherits(parents, "compounds"))
    {
        compTab <- as.data.table(parents)
        if (!is.null(compTab[["compoundName"]]))
            compTab[, name := ifelse(nzchar(compoundName), compoundName, identifier)]
        else
            setnames(compTab, "identifier", "name")
        parents <- compTab[, intersect(c("name", "SMILES", "InChI", "InChIKey", "neutral_formula", "neutralMass",
                                         "molNeutralized"), names(compTab)), with = FALSE]
        setnames(parents, "neutral_formula", "formula")
        parents <- prepareSuspectList(parents, NULL, skipInvalid, checkDesc = FALSE,
                                      prefCalcChemProps = prefCalcChemProps, neutralChemProps = neutralChemProps,
                                      calcMZs = FALSE)
    }
    else # suspect screening
    {
        parents <- copy(screenInfo(parents))
        parents <- unique(parents, by = "name")
        # only keep columns that are relevant for suspect lists (ie when convertToSuspects is called)
        keepCols <- c("name", "rt", "formula", "SMILES", "InChI", "InChIKey", "neutralMass", "molNeutralized",
                      "fragments_mz", "fragments_formula", "adduct")
        parents <- parents[, intersect(keepCols, names(parents)), with = FALSE]
    }
    
    if (is.null(parents[[checkWhat]]))
        stop(sprintf("No %s information available for parents. Please include either %s columns.", checkWhat,
                     if (checkWhat == "formula") "formula or SMILES/InChI" else "SMILES or InChI"), call. = FALSE)
    
    noData <- is.na(parents[[checkWhat]]) | !nzchar(parents[[checkWhat]])
    if (any(noData))
    {
        do.call(if (skipInvalid) warning else stop,
                list(sprintf("The following parents miss mandatory %s: %s", checkWhat,
                             paste0(parents$name[noData], collapse = ","))))
        parents <- parents[!noData]
    }
    
    return(parents)
}

doConvertToMFDB <- function(prodAll, parents, out, includeParents)
{
    # UNDONE: cache?

    if (nrow(prodAll) == 0)
        stop("Cannot create MetFrag database: no data", call. = FALSE)

    prodAll <- unique(prodAll, by = "name") # omit duplicates from the same parent
    
    # combine equal TPs from different parents
    prodAll[, c("name", "parent") := .(paste0(name, collapse = ","), paste0(parent, collapse = ",")), by = "InChIKey"]
    
    # ... and remove now duplicates
    prodAll <- unique(prodAll, by = "InChIKey")
    
    # set to MetFrag style names
    setnames(prodAll,
             c("name", "formula", "neutralMass"),
             c("Identifier", "MolecularFormula", "MonoisotopicMass"))
    
    if (includeParents)
    {
        pars <- copy(parents)
        setnames(pars,
                 c("name", "formula", "neutralMass"),
                 c("Identifier", "MolecularFormula", "MonoisotopicMass"))
        pars[, CompoundName := Identifier]
        prodAll <- rbind(pars, prodAll, fill = TRUE)
    }
    
    # Add required InChIKey1 column
    prodAll[, InChIKey1 := getIKBlock1(InChIKey)]
    
    # equalize identifiers and names
    prodAll[, CompoundName := Identifier]
    
    keepCols <- c("Identifier", "MolecularFormula", "MonoisotopicMass", "SMILES", "InChI", "InChIKey", "InChIKey1",
                  "molNeutralized", "ALogP", "LogP", "XLogP", "parent")
    
    fwrite(prodAll[, intersect(keepCols, names(prodAll)), with = FALSE], out)
}

doPlotTPGraph <- function(TPTab, parents, cmpTab, structuresMax, prune, onlyCompletePaths, width, height)
{
    # UNDONE: don't make name unique, but use IDs?
    
    if (nrow(TPTab) == 0)
        stop("No TPs to plot", call. = FALSE)
    
    TPTab <- copy(TPTab)
    TPTab[, c("name_orig", "name") := .(name, make.unique(name))]
    TPTab[, parent_name := fifelse(is.na(parent_ID), parent, name[match(parent_ID, ID)]), by = "parent"]
    
    if (!is.null(cmpTab))
    {
        TPTab <- TPTab[parent %chin% cmpTab$parent_name] # omit missing root parents
        TPTab[, present := name_orig %chin% cmpTab$TP_name]
        
        TPTab[, childPresent := FALSE]
        markChildPresent <- function(TPNames)
        {
            if (length(TPNames) == 0)
                return()
            TPTab[name %chin% TPNames, childPresent := TRUE]
            pars <- TPTab[name %chin% TPNames]$parent_name
            markChildPresent(pars[TPTab[name %chin% pars]$childPresent == FALSE])
        }
        markChildPresent(TPTab[present == TRUE]$parent_name)
        
        if (prune)
            TPTab <- TPTab[present == TRUE | childPresent == TRUE]
        if (onlyCompletePaths)
        {
            TPTab <- TPTab[present == TRUE]
            # keep removing TPs without parent until no change
            oldn <- nrow(TPTab)
            repeat
            {
                TPTab <- TPTab[parent_name == parent | parent_name %chin% name]
                newn <- nrow(TPTab)
                if (oldn == newn)
                    break
                oldn <- newn
            }
        }
    }
    
    TPTab[, parent_formula := fifelse(is.na(parent_ID),
                                      parents$formula[match(parent_name, parents$name)],
                                      formula[match(parent_name, name)])]
    TPTab[, formulaDiff := mapply(formula, parent_formula, FUN = getFormulaDiffText)]
    
    nodes <- data.table(id = union(TPTab$parent, TPTab$name))
    nodes[, isTP := id %chin% TPTab$name]
    nodes[isTP == TRUE, label := paste0("TP", TPTab$chem_ID[match(id, TPTab$name)])]
    nodes[isTP == FALSE, label := id]
    nodes[, group := if (.N > 1) label else "unique", by = "label"]
    nodes[, present := isTP == FALSE | TPTab$present[match(id, TPTab$name)]]
    nodes[present == TRUE, shapeProperties := list(list(list(useBorderWithImage = TRUE)))]
    nodes[present == FALSE, shapeProperties := list(list(list(useBorderWithImage = FALSE)))]
    nodes[, present := NULL]
    nodes[isTP == FALSE, level := 0]
    nodes[isTP == TRUE, level := TPTab$generation[match(id, TPTab$name)]]
    
    if (nrow(nodes) <= structuresMax && nrow(nodes) > 0)
    {
        # UNDONE: make util?
        imgf <- tempfile(fileext = ".svg") # temp file is re-used
        getURIFromSMILES <- function(SMILES)
        {
            mol <- getMoleculesFromSMILES(SMILES, emptyIfFails = TRUE)[[1]]
            saveRCDKStructure(mol, "svg", imgf, 500, 500, transparent = FALSE)
            return(knitr::image_uri(imgf))
        }
        nodes[, shape := "image"]
        nodes[, SMILES := fifelse(isTP, TPTab$SMILES[match(id, TPTab$name)], parents$SMILES[match(id, parents$name)])]
        nodes[, image := getURIFromSMILES(SMILES[1]), by = "SMILES"]
        nodes[, SMILES := NULL]
    }
    else
        nodes[, shape := "ellipse"]
    
    TPCols <- intersect(c("name", "name_lib", "SMILES", "formula", "generation", "accumulation", "production",
                          "globalAccumulation", "likelihood", "Lipinski_Violations", "Insecticide_Likeness_Violations",
                          "Post_Em_Herbicide_Likeness_Violations", "transformation", "transformation_ID", "enzyme",
                          "biosystem", "evidencedoi", "evidencedref", "sourcecomment", "datasetref", "similarity",
                          "mergedBy", "coverage"), names(TPTab))
    nodes[isTP == TRUE, title := sapply(id, function(TP)
    {
        TPTabSub <- TPTab[name == TP, TPCols, with = FALSE]
        return(paste0(names(TPTabSub), ": ", TPTabSub, collapse = "<br>"))
    })]
    
    edges <- data.table(from = TPTab$parent_name, to = TPTab$name, label = TPTab$formulaDiff)
    
    visNetwork::visNetwork(nodes = nodes, edges = edges, width = width, height = height) %>%
        visNetwork::visNodes(shapeProperties = list(useBorderWithImage = FALSE)) %>%
        visNetwork::visEdges(arrows = "to", font = list(align = "top", size = 12)) %>%
        visNetwork::visOptions(highlightNearest = list(enabled = TRUE, hover = TRUE, algorithm = "hierarchical"),
                               selectedBy = list(variable = "group", main = "Select duplicate TPs",
                                                 values = unique(nodes$group[nodes$group != "unique"]))) %>%
        visNetwork::visHierarchicalLayout(enabled = TRUE, sortMethod = "directed")
}

getProductsFromLib <- function(TPLibrary, TPLibrarySel, generations, matchGenerationsBy, matchIDBy)
{
    results <- split(TPLibrarySel, by = "parent_name")
    
    curTPIDs <- setNames(vector("integer", length = length(results)), names(results))
    prepTPs <- function(r, pn, pid, gen, prvLogPDiff)
    {
        # remove parent columns
        set(r, j = grep("^parent_", names(r), value = TRUE), value = NULL)
        
        # remove TP_ prefix
        cols <- grep("^TP_", names(r), value = TRUE)
        setnames(r, cols, sub("^TP_", "", cols))
        
        r[, retDir := 0] # may be changed below
        r[, generation := gen]
        
        r[, ID := curTPIDs[pn] + seq_len(nrow(r))]
        curTPIDs[pn] <<- curTPIDs[pn] + nrow(r)
        r[, parent_ID := pid]
        
        # make it additive so LogPDiff corresponds to the original parent
        if (!is.null(r[["LogPDiff"]]) && !is.null(prvLogPDiff))
            r[, LogPDiff := LogPDiff + prvLogPDiff]
        
        return(r)
    }
    results <- Map(results, names(results), f = prepTPs, MoreArgs = list(pid = NA_integer_, gen = 1, prvLogPDiff = NULL))
    
    if (generations > 1)
    {
        for (gen in seq(2, generations))
        {
            results <- Map(results, names(results), f = function(r, pn)
            {
                tps <- r[generation == (gen-1)]
                nexttps <- rbindlist(lapply(split(tps, seq_len(nrow(tps))), function(tpRow)
                {
                    nt <- copy(TPLibrary[get(paste0("parent_", matchGenerationsBy)) == tpRow[[matchGenerationsBy]]])
                    return(prepTPs(nt, pn, tpRow$ID, gen, tpRow$LogPDiff))
                }))
                return(rbind(r, nexttps))
            })
        }
    }
    
    # fill in chem IDs and names now that we sorted out all TPs
    results <- Map(results, names(results), f = function(r, pn)
    {
        # prune TPs that equal the main parent
        # UNDONE: this does not work if matchGenerationsBy=="name" (default for formulas) 
        r <- r[generation == 1 | get(matchGenerationsBy) != TPLibrary[parent_name == pn][[paste0("parent_", matchGenerationsBy)]][1]]
        set(r, j = "chem_ID", value = match(r[[matchIDBy]], unique(r[[matchIDBy]])))
        setnames(r, "name", "name_lib")
        set(r, j = "name", value = paste0(pn, "-TP", r$chem_ID))
    })

    if (!is.null(TPLibrarySel[["retDir"]]))
    {
        results <- Map(results, TPLibrarySel[match(names(results), parent_name)]$retDir,
                       f = data.table::set, MoreArgs = list(i = NULL, j = "retDir"))
    }
    else if (!is.null(TPLibrarySel[["parent_LogP"]]) && !is.null(TPLibrarySel[["TP_LogP"]]))
    {
        results <- Map(results, TPLibrarySel[match(names(results), parent_name)]$parent_LogP,
                       f = function(r, pLogP) set(r, j = "retDir", value = fifelse(r$LogP < pLogP, -1, 1)))
    }
    else if (!is.null(TPLibrarySel[["LogPDiff"]]))
    {
        results <- lapply(results, function(x) set(x, j = "retDir", value = fcase(x$LogPDiff < 0, -1,
                                                                                  x$LogPDiff > 0, 1,
                                                                                  default = 0)))
    }
    
    results <- pruneList(results, checkZeroRows = TRUE)

    return(results)
}

prepareDataForTPLibrary <- function(parents, TPLibrary, generations, matchParentsBy, matchGenerationsBy, matchIDBy,
                                    neutralizeTPs)
{
    TPLibrary <- copy(as.data.table(TPLibrary))
    
    # prepare chem infos
    for (wh in c("parent", "TP"))
    {
        cols <- c("SMILES", "InChI", "InChIKey", "formula", "neutralMass", "molNeutralized")
        
        # temporarily remove parent/TP prefix for prepareChemTable
        whcols <- intersect(paste0(wh, "_", cols), names(TPLibrary))
        setnames(TPLibrary, whcols, sub(paste0(wh, "_"), "", whcols))
        TPLibrary <- prepareChemTable(TPLibrary, FALSE, neutralizeTPs)
        
        # put back prefix
        regcols <- intersect(cols, names(TPLibrary))
        # clean out NA columns
        for (col in regcols)
        {
            if (all(is.na(TPLibrary[[col]])))
                TPLibrary[, (col) := NULL]
        }
        regcols <- intersect(regcols, names(TPLibrary))
        setnames(TPLibrary, regcols, paste0(wh, "_", regcols))
        
        if ("InChIKey1" %in% c(matchParentsBy, matchGenerationsBy))
        {
            whcol <- paste0(wh, "_InChIKey")
            if (is.null(TPLibrary[[whcol]]))
                stop(sprintf("Cannot match by InChIKey1: missing %s column in the library", whcol), call. = FALSE)
            TPLibrary[, (paste0(whcol, 1)) := getIKBlock1(get(whcol))]
        }
        
        for (mb in union(matchParentsBy, matchGenerationsBy))
        {
            whcol <- paste0(wh, "_", mb)
            if (is.null(TPLibrary[[whcol]]))
                stop(sprintf("Cannot match by %s: missing %s column in the library", mb, whcol), call. = FALSE)
        }
    }
    
    TPLibrarySel <- TPLibrary
    if (!is.null(parents))
    {
        # match with library
        
        dataLib <- TPLibrary[[paste0("parent_", matchParentsBy)]]
        dataSusp <- parents[[matchParentsBy]]
        
        if (matchParentsBy != "name")
        {
            # rename from suspect list
            TPLibrary[, parent_name_lib := parent_name] # store original
            TPLibrary[, parent_name := parents[match(dataLib, dataSusp)]$name]
        }
        
        # only take data in both
        dataInBoth <- intersect(dataLib, dataSusp)
        TPLibrarySel <- TPLibrary[dataLib %chin% dataInBoth]
        parents <- parents[dataSusp %chin% dataInBoth]
    }
    else
    {
        parents <- unique(TPLibrary[, grepl("^parent_", names(TPLibrary)), with = FALSE], by = "parent_name")
        setnames(parents, sub("^parent_", "", names(parents)))
    }
    
    products <- getProductsFromLib(TPLibrary, TPLibrarySel, generations, matchGenerationsBy, matchIDBy)
    parents <- parents[name %in% names(products)]
    products <- products[match(parents$name, names(products))] # sync order
    
    return(list(parents = parents, products = products))
}
