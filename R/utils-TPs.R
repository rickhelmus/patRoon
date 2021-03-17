getTPSuspects <- function(suspects, adduct, skipInvalid)
{
    if (is.data.frame(suspects))
        suspects <- prepareSuspectList(suspects, adduct, skipInvalid, calcMZs = FALSE)
    else if (inherits(suspects, "compounds"))
    {
        compTab <- as.data.table(suspects)
        if (!is.null(compTab[["compoundName"]]))
            compTab[, name := ifelse(nzchar(compoundName), compoundName, identifier)]
        else
            setnames(compTab, "Identifier", "name")
        suspects <- compTab[, c("name", "SMILES", "InChI", "InChIKey"), with = FALSE]
    }
    else # suspect screening
        suspects <- copy(screenInfo(suspects)) # UNDONE: keep all columns?
    
    if (is.null(suspects[["SMILES"]]))
        stop("No SMILES information available for suspects. Please include either SMILES or InChI columns.")
    
    noSM <- is.na(suspects$SMILES) | !nzchar(suspects$SMILES)
    if (any(noSM))
    {
        do.call(if (skipInvalid) warning else stop,
                "The following suspects miss mandatory SMILES: ", paste0(suspects$name[noSM], collapse = ","))
        suspects <- suspects[!noSM]
    }
    
    return(suspects)
}

doConvertToMFDB <- function(predAll, suspects, out, includePrec)
{
    # UNDONE: cache?
    
    # set to MetFrag style names
    setnames(predAll,
             c("name", "formula", "neutralMass"),
             c("Identifier", "MolecularFormula", "MonoisotopicMass"))
    if (!is.null(predAll[["Precursor Major Isotope Mass"]])) # BT
        setnames(predAll, "Precursor Major Isotope Mass", "Precursor MonoisotopicMass")
    
    if (includePrec)
    {
        precs <- copy(suspects)
        setnames(precs,
                 c("name", "formula", "neutralMass"),
                 c("Identifier", "MolecularFormula", "MonoisotopicMass"))
        precs[, CompoundName := Identifier]
        predAll <- rbind(precs, predAll, fill = TRUE)
    }
    
    # Add required InChIKey1 column
    predAll[, InChIKey1 := getIKBlock1(InChIKey)]
    
    # equalize identifiers and names
    predAll[, CompoundName := Identifier]
    
    keepCols <- c("Identifier", "MolecularFormula", "MonoisotopicMass", "Precursor MonoisotopicMass",
                  "SMILES", "InChI", "InChIKey", "InChIKey1", "ALogP") # UNDONE: more?
    
    fwrite(predAll[, intersect(keepCols, names(predAll)), with = FALSE], out)
}
