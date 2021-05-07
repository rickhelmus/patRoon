#' @include main.R
#' @include TP.R
NULL

#' @export
transformationProductsLibrary <- setClass("transformationProductsLibrary", contains = "transformationProducts")

setMethod("initialize", "transformationProductsLibrary",
          function(.Object, ...) callNextMethod(.Object, algorithm = "library", ...))


#' @export
generateTPsLibrary <- function(parents = NULL, TPLibrary = NULL, adduct = NULL, skipInvalid = TRUE,
                               matchSuspectsBy = "InChIKey")
{
    # UNDONE: default match by IK or IK1?
    
    checkmate::assert(
        checkmate::checkNull(parents),
        checkmate::checkClass(parents, "data.frame"),
        checkmate::checkClass(parents, "compounds"),
        checkmate::checkClass(parents, "featureGroupsScreening"),
        checkmate::checkClass(parents, "featureGroupsScreeningSet"),
        .var.name = "parents"
    )
    
    ac <- checkmate::makeAssertCollection()
    if (!is.null(TPLibrary))
    {
        checkmate::assertDataFrame(TPLibrary, any.missing = FALSE, min.rows = 1, add = ac)
        libCols <- c("parent_name", "parent_SMILES", "TP_name", "TP_SMILES")
        assertHasNames(TPLibrary, libCols, add = ac)
        for (cl in libCols)
            assertListVal(TPLibrary, cl, checkmate::assertCharacter, min.chars = 1, any.missing = FALSE, add = ac)
    }
        
    if (is.data.frame(parents))
        assertSuspectList(parents, needsAdduct = FALSE, skipInvalid = TRUE, add = ac)
    checkmate::assertFlag(skipInvalid, add = ac)
    checkmate::assertChoice(matchSuspectsBy, c("InChIKey", "InChIKey1", "InChI", "SMILES"), null.ok = FALSE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (is.null(TPLibrary))
        TPLibrary <- copy(PubChemTransformations) # default to embedded PC transformations
    else
    {
        TPLibrary <- copy(as.data.table(TPLibrary))
        
        # add chem infos where necessary
        for (wh in c("parent", "TP"))
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

    if (!is.null(parents))
    {
        parents <- getTPParents(parents, adduct, skipInvalid)
        
        # match with library
        if (matchSuspectsBy == "InChIKey1")
        {
            dataLib <- getIKBlock1(TPLibrary$parent_InChIKey)
            dataSusp <- getIKBlock1(parents$InChIKey)
        }
        else
        {
            dataLib <- TPLibrary[[paste0("parent_", matchSuspectsBy)]]
            dataSusp <- parents[[matchSuspectsBy]]
        }
        
        # only take data in both
        dataInBoth <- intersect(dataLib, dataSusp)
        dataInBoth <- dataInBoth[match(dataSusp, dataInBoth, nomatch = 0)] # align suspect order
        TPLibrary <- TPLibrary[match(dataInBoth, dataLib)] # match and align suspect order
        parents <- parents[dataSusp %chin% dataInBoth]
        
        # rename from suspect list
        TPLibrary[, parent_name := parents$name]
    }
    else
    {
        parents <- unique(TPLibrary[, grepl("^parent_", names(TPLibrary)), with = FALSE], by = "parent_name")
        setnames(parents, sub("^parent_", "", names(parents)))
    }
    
    results <- split(TPLibrary, by = "parent_name")
    results <- Map(results, names(results), f = function(r, pn)
    {
        # remove parent columns
        set(r, j = grep("^parent_", names(r), value = TRUE), value = NULL)
        
        # remove TP_ prefix
        cols <- grep("^TP_", names(r), value = TRUE)
        setnames(r, cols, sub("^TP_", "", cols))
        
        # make TP names unique
        r[, name := paste0(pn, "-TP-", name)]
        
        r[, retDir := 0] # may be changed below
        
        return(r)
    })
    
    if (!is.null(TPLibrary[["parent_LogP"]]) && !is.null(TPLibrary[["TP_LogP"]]))
    {
        results <- Map(results, TPLibrary[match(names(results), parent_name)]$parent_LogP,
                       f = function(r, pLogP) set(r, j = "retDir", value = fifelse(r$LogP < pLogP, -1, 1)))
    }
    
    results <- pruneList(results, checkZeroRows = TRUE)
    parents <- parents[name %in% names(results)]
    
    return(transformationProductsLibrary(parents = parents, products = results))
}

#' @export
setMethod("convertToMFDB", "transformationProductsLibrary", function(TPs, out, includeParents)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertPathForOutput(out, overwrite = TRUE, add = ac) # NOTE: assert doesn't work on Windows...
    checkmate::assertFlag(includeParents, add = ac)
    checkmate::reportAssertions(ac)
    
    allTPs <- rbindlist(TPs@products)
    
    doConvertToMFDB(allTPs, parents(TPs), out, includeParents)
})

setMethod("linkParentsToFGroups", "transformationProductsLibrary", function(TPs, fGroups)
{
    return(screenInfo(fGroups)[name %in% names(TPs), c("name", "group"), with = FALSE])
})
