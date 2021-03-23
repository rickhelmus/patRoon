#' @include main.R
#' @include TP.R
NULL

#' @export
TPPredictionsLibrary <- setClass("TPPredictionsLibrary", contains = "TPPredictions")

setMethod("initialize", "TPPredictionsLibrary",
          function(.Object, ...) callNextMethod(.Object, algorithm = "library", ...))


#' @export
predictTPsLibrary <- function(suspects = NULL, TPLibrary = NULL, adduct = NULL, skipInvalid = TRUE,
                              matchSuspectsBy = "InChIKey")
{
    # UNDONE: default match by IK or IK1?
    
    checkmate::assert(
        checkmate::checkNull(suspects),
        checkmate::checkClass(suspects, "data.frame"),
        checkmate::checkClass(suspects, "compounds"),
        checkmate::checkClass(suspects, "featureGroupsScreening"),
        checkmate::checkClass(suspects, "featureGroupsScreeningSet"),
        .var.name = "suspects"
    )
    
    ac <- checkmate::makeAssertCollection()
    if (!is.null(TPLibrary))
    {
        checkmate::assertDataFrame(TPLibrary, any.missing = FALSE, min.rows = 1, add = ac)
        libCols <- c("precursor_name", "precursor_SMILES", "TP_name", "TP_SMILES")
        assertHasNames(TPLibrary, libCols, add = ac)
        for (cl in libCols)
            assertListVal(TPLibrary, cl, checkmate::assertCharacter, min.chars = 1, any.missing = FALSE,
                          unique = cl == "precursor_name", add = ac)
    }
        
    if (is.data.frame(suspects))
        assertSuspectList(suspects, needsAdduct = FALSE, skipInvalid = TRUE, add = ac)
    checkmate::assertFlag(skipInvalid, add = ac)
    checkmate::assertChoice(matchSuspectsBy, c("InChIKey", "InChIKey1", "InChI", "SMILES"), null.ok = FALSE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (is.null(TPLibrary))
        TPLibrary <- copy(PubChemTransformations) # default to embedded PC transformations
    else
    {
        TPLibrary <- copy(as.data.table(TPLibrary))
        
        # add chem infos where necessary
        for (wh in c("precursor", "TP"))
        {
            for (col in c("formula", "InChI", "InChIKey"))
            {
                whcol <- paste0(wh, "_", col)
                if (is.null(TPLibrary[[whcol]]))
                {
                    whSMI <- paste0(wh, "_SMILES")
                    TPLibrary[, (whcol) := switch(col,
                                                  formula = convertToFormulaBabel(get(whSMI), "smi"),
                                                  InChI = babelConvert(get(whSMI), "smi", "inchi"),
                                                  InChIKey = babelConvert(get(whSMI), "smi", "inchikey"))]
                }
            }
            whmcol <- paste0(wh, "_neutralMass")
            if (is.null(TPLibrary[[whmcol]]))
                TPLibrary[, (whmcol) := sapply(get(paste0(wh, "_formula")), getFormulaMass)]
        }
    }

    if (!is.null(suspects))
    {
        suspects <- getTPSuspects(suspects, adduct, skipInvalid)
        
        # match with library
        if (matchSuspectsBy == "InChIKey1")
        {
            dataLib <- getIKBlock1(TPLibrary$precursor_InChIKey)
            dataSusp <- getIKBlock1(suspects$InChIKey)
        }
        else
        {
            dataLib <- TPLibrary[[paste0("precursor_", matchSuspectsBy)]]
            dataSusp <- suspects[[matchSuspectsBy]]
        }
        
        # only take data in both
        dataInBoth <- intersect(dataLib, dataSusp)
        dataInBoth <- dataInBoth[match(dataSusp, dataInBoth, nomatch = 0)] # align suspect order
        TPLibrary <- TPLibrary[match(dataInBoth, dataLib)] # match and align suspect order
        suspects <- suspects[dataSusp %chin% dataInBoth]
        
        # rename from suspect list
        TPLibrary[, precursor_name := suspects$name]
    }
    else
    {
        suspects <- unique(TPLibrary[, grepl("^precursor_", names(TPLibrary)), with = FALSE], by = "precursor_name")
        setnames(suspects, sub("^precursor_", "", names(suspects)))
    }
    
    results <- split(TPLibrary, by = "precursor_name")
    results <- Map(results, names(results), f = function(r, pn)
    {
        # remove precursor columns
        set(r, j = grep("^precursor_", names(r), value = TRUE), value = NULL)
        
        # remove TP_ prefix
        cols <- grep("^TP_", names(r), value = TRUE)
        setnames(r, cols, sub("^TP_", "", cols))
        
        # make TP names unique
        r[, name := paste0(pn, "-TP-", name)]
        
        r[, RTDir := 0] # may be changed below
        
        return(r)
    })
    
    if (!is.null(TPLibrary[["precursor_LogP"]]) && !is.null(TPLibrary[["TP_LogP"]]))
    {
        results <- Map(results, TPLibrary$precursor_LogP,
                       f = function(r, pLogP) set(r, j = "RTDir", value = fifelse(r$LogP < pLogP, -1, 1)))
    }
    
    results <- pruneList(results, checkZeroRows = TRUE)
    suspects <- suspects[name %in% names(results)]
    
    return(TPPredictionsLibrary(suspects = suspects, predictions = results))
}

#' @export
setMethod("convertToMFDB", "TPPredictionsLibrary", function(pred, out, includePrec)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertPathForOutput(out, overwrite = TRUE, add = ac) # NOTE: assert doesn't work on Windows...
    checkmate::assertFlag(includePrec, add = ac)
    checkmate::reportAssertions(ac)
    
    allTPs <- rbindlist(pred@predictions)
    
    doConvertToMFDB(allTPs, suspects(pred), out, includePrec)
})

setMethod("linkPrecursorsToFGroups", "TPPredictionsLibrary", function(pred, fGroups)
{
    return(screenInfo(fGroups)[name %in% names(pred), c("name", "group"), with = FALSE])
})
