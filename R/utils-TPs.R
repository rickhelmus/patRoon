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
            setnames(compTab, "identifier", "name")
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
                list("The following parents miss mandatory SMILES: ", paste0(parents$name[noSM], collapse = ",")))
        parents <- parents[!noSM]
    }
    
    return(parents)
}

doConvertToMFDB <- function(prodAll, parents, out, includeParents)
{
    # UNDONE: cache?
    
    # set to MetFrag style names
    setnames(prodAll,
             c("name", "formula", "neutralMass"),
             c("Identifier", "MolecularFormula", "MonoisotopicMass"))
    if (!is.null(prodAll[["parent_Major Isotope Mass"]])) # BT
        setnames(prodAll, "parent_Major Isotope Mass", "Parent MonoisotopicMass")
    
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
    
    keepCols <- c("Identifier", "MolecularFormula", "MonoisotopicMass", "Precursor MonoisotopicMass",
                  "SMILES", "InChI", "InChIKey", "InChIKey1", "ALogP") # UNDONE: more?
    
    fwrite(prodAll[, intersect(keepCols, names(prodAll)), with = FALSE], out)
}
