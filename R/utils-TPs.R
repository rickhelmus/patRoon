getTPParents <- function(parents, adduct, skipInvalid)
{
    if (is.data.frame(parents))
        parents <- prepareSuspectList(parents, adduct, skipInvalid, calcMZs = FALSE)
    else if (inherits(parents, "compounds"))
    {
        compTab <- as.data.table(parents)
        if (!is.null(compTab[["compoundName"]]))
            compTab[, name := ifelse(nzchar(compoundName), compoundName, identifier)]
        else
            setnames(compTab, "Identifier", "name")
        parents <- compTab[, c("name", "SMILES", "InChI", "InChIKey"), with = FALSE]
    }
    else # suspect screening
        parents <- copy(screenInfo(parents)) # UNDONE: keep all columns?
    
    if (is.null(parents[["SMILES"]]))
        stop("No SMILES information available for parents. Please include either SMILES or InChI columns.")
    
    noSM <- is.na(parents$SMILES) | !nzchar(parents$SMILES)
    if (any(noSM))
    {
        do.call(if (skipInvalid) warning else stop,
                "The following parents miss mandatory SMILES: ", paste0(parents$name[noSM], collapse = ","))
        parents <- parents[!noSM]
    }
    
    return(parents)
}

doConvertToMFDB <- function(predAll, parents, out, includeParents)
{
    # UNDONE: cache?
    
    # set to MetFrag style names
    setnames(predAll,
             c("name", "formula", "neutralMass"),
             c("Identifier", "MolecularFormula", "MonoisotopicMass"))
    if (!is.null(predAll[["Precursor Major Isotope Mass"]])) # BT
        setnames(predAll, "Precursor Major Isotope Mass", "Precursor MonoisotopicMass")
    
    if (includeParents)
    {
        pars <- copy(parents)
        setnames(pars,
                 c("name", "formula", "neutralMass"),
                 c("Identifier", "MolecularFormula", "MonoisotopicMass"))
        pars[, CompoundName := Identifier]
        predAll <- rbind(pars, predAll, fill = TRUE)
    }
    
    # Add required InChIKey1 column
    predAll[, InChIKey1 := getIKBlock1(InChIKey)]
    
    # equalize identifiers and names
    predAll[, CompoundName := Identifier]
    
    keepCols <- c("Identifier", "MolecularFormula", "MonoisotopicMass", "Precursor MonoisotopicMass",
                  "SMILES", "InChI", "InChIKey", "InChIKey1", "ALogP") # UNDONE: more?
    
    fwrite(predAll[, intersect(keepCols, names(predAll)), with = FALSE], out)
}
