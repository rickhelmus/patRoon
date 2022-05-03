getTPParents <- function(parents, skipInvalid)
{
    if (is.data.frame(parents))
        parents <- prepareSuspectList(parents, NULL, skipInvalid, calcMZs = FALSE)
    else if (inherits(parents, "compounds"))
    {
        compTab <- as.data.table(parents)
        if (!is.null(compTab[["compoundName"]]))
            compTab[, name := ifelse(nzchar(compoundName), compoundName, identifier)]
        else
            setnames(compTab, "identifier", "name")
        parents <- compTab[, c("name", "SMILES", "InChI", "InChIKey"), with = FALSE]
        parents <- prepareSuspectList(parents, NULL, skipInvalid, calcMZs = FALSE)
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
