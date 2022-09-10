getTPParents <- function(parents, skipInvalid, prefCalcChemProps)
{
    if (is.data.frame(parents))
        parents <- prepareSuspectList(parents, NULL, skipInvalid, checkDesc = TRUE,
                                      prefCalcChemProps = prefCalcChemProps, calcMZs = FALSE)
    else if (inherits(parents, "compounds"))
    {
        compTab <- as.data.table(parents)
        if (!is.null(compTab[["compoundName"]]))
            compTab[, name := ifelse(nzchar(compoundName), compoundName, identifier)]
        else
            setnames(compTab, "identifier", "name")
        parents <- compTab[, c("name", "SMILES", "InChI", "InChIKey", "neutral_formula", "neutralMass"), with = FALSE]
        setnames(parents, "neutral_formula", "formula")
        parents <- prepareSuspectList(parents, NULL, skipInvalid, checkDesc = FALSE,
                                      prefCalcChemProps = prefCalcChemProps, calcMZs = FALSE)
    }
    else # suspect screening
    {
        parents <- copy(screenInfo(parents))
        parents <- unique(parents, by = "name")
        # only keep columns that are relevant for suspect lists (ie when convertToSuspects is called)
        keepCols <- c("name", "rt", "formula", "SMILES", "InChI", "InChIKey", "neutralMass",
                      "fragments_mz", "fragments_formula", "adduct")
        parents <- parents[, intersect(keepCols, names(parents)), with = FALSE]
    }
    
    if (is.null(parents[["SMILES"]]))
        stop("No SMILES information available for parents. Please include either SMILES or InChI columns.")
    
    noSM <- is.na(parents$SMILES) | !nzchar(parents$SMILES)
    if (any(noSM))
    {
        do.call(if (skipInvalid) warning else stop,
                list("The following parents miss mandatory SMILES: ", paste0(parents$name[noSM], collapse = ",")))
        parents <- parents[!noSM]
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
                  "ALogP", "LogP", "XLogP", "parent")
    
    fwrite(prodAll[, intersect(keepCols, names(prodAll)), with = FALSE], out)
}

doPlotTPGraph <- function(TPTab, parents, cmpTab, structuresMax, prune, onlyCompletePaths)
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
    TPTab[, formulaDiff := mapply(formula, parent_formula, FUN = function(f, pf)
    {
        sfl <- splitFormulaToList(subtractFormula(f, pf))
        ret <- ""
        subfl <- sfl[sfl < 0]
        if (length(subfl) > 0)
            ret <- paste0("-", formulaListToString(abs(subfl)))
        addfl <- sfl[sfl > 0]
        if (length(addfl) > 0)
            ret <- if (nzchar(ret)) paste0(ret, " +", formulaListToString(addfl)) else paste0("+", formulaListToString(addfl))
        return(ret)
    })]
    
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
        imgf <- tempfile(fileext = ".png") # temp file is re-used
        getURIFromSMILES <- function(SMILES)
        {
            mol <- getMoleculesFromSMILES(SMILES, emptyIfFails = TRUE)[[1]]
            withr::with_png(imgf, withr::with_par(list(mar = rep(0, 4)), plot(getRCDKStructurePlot(mol, 150, 150))))
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
    
    visNetwork::visNetwork(nodes = nodes, edges = edges) %>%
        visNetwork::visNodes(shapeProperties = list(useBorderWithImage = FALSE)) %>%
        visNetwork::visEdges(arrows = "to", font = list(align = "top", size = 12)) %>%
        visNetwork::visOptions(highlightNearest = list(enabled = TRUE, hover = TRUE, algorithm = "hierarchical"),
                               selectedBy = list(variable = "group", main = "Select duplicate TPs",
                                                 values = unique(nodes$group[nodes$group != "unique"]))) %>%
        visNetwork::visHierarchicalLayout(enabled = TRUE, sortMethod = "directed")
}
