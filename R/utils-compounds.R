#' @include main.R
#' @include compounds.R
#' @include mspeaklists.R
#' @include formulas.R
NULL

getIKBlock1 <- function(IK) strtrim(IK, 14)

mergeFragInfo <- function(fiLeft, fiRight, leftName, rightName)
{
    # UNDONE: what about multiple formula candidates?

    fiLeft <- copy(fiLeft); fiRight <- copy(fiRight)

    if (is.null(fiLeft[["mergedBy"]]))
        fiLeft[, mergedBy := leftName]
    if (is.null(fiRight[["mergedBy"]]))
        fiRight[, mergedBy := rightName]

    if (nrow(fiLeft) == 0)
        fiLeft <- fiRight
    else if (nrow(fiRight) > 0)
    {
        # for overlap: just add label
        fiLeft <- merge(fiLeft, fiRight[, c("PLIndex", "mergedBy"), with = FALSE], all.x = TRUE, by = "PLIndex")
        fiLeft[is.na(mergedBy.y), mergedBy := mergedBy.x]
        fiLeft[is.na(mergedBy.x), mergedBy := mergedBy.y]
        fiLeft[!is.na(mergedBy.x) & !is.na(mergedBy.y), mergedBy := paste(mergedBy.x, mergedBy.y, sep = ",")]
        fiLeft[, c("mergedBy.x", "mergedBy.y") := NULL]

        # add unique
        fiUnique <- fiRight[!PLIndex %in% fiLeft$PLIndex]
        if (nrow(fiUnique) > 0)
            fiLeft <- rbind(fiLeft, fiUnique, fill = TRUE)
    }

    return(fiLeft)
}

#' @details \code{compoundScorings} displays an overview of scorings may be
#'   applied to rank candidate compounds (see \verb{Scorings} section below).
#'
#' @param includeSuspectLists,onlyDefault,includeNoDB A logical specifying
#'   whether scoring terms related to suspect lists, default scoring terms and
#'   non-database specific scoring terms should be included in the output,
#'   respectively.
#'
#' @return \code{compoundScorings} returns a \code{data.frame} with information
#'   on which scoring terms are used, what their algorithm specific name is and
#'   other information such as to which database they apply and short remarks.
#' @rdname compound-generation
#' @export
compoundScorings <- function(algorithm = NULL, database = NULL, includeSuspectLists = TRUE,
                             onlyDefault = FALSE, includeNoDB = TRUE)
{
    algos <- c("metfrag", "sirius")

    ret <- patRoon:::compScorings # stored inside R/sysdata.rda

    ac <- checkmate::makeAssertCollection()
    checkmate::assertChoice(algorithm, algos, null.ok = TRUE, add = ac)
    checkmate::assertString(database, na.ok = FALSE, null.ok = TRUE, add = ac)
    checkmate::assertFlag(includeSuspectLists, add = ac)
    checkmate::reportAssertions(ac)

    if (!is.null(algorithm))
        ret <- ret[nzchar(ret[[algorithm]]), names(ret) != setdiff(algos, algorithm)]
    if (!is.null(database))
    {
        if (includeNoDB)
            ret <- ret[!nzchar(ret$database) | ret$database == database, ]
        else
            ret <- ret[ret$database == database, ]
    }
    if (!includeSuspectLists)
        ret <- ret[!ret$suspect_list, ]
    if (onlyDefault)
        ret <- ret[ret$default, ]

    return(ret)
}

compScoreNames <- function(onlyNums) unique(compoundScorings(includeSuspectLists = !onlyNums)$name)

makeDBIdentLink <- function(db, ident)
{
    ident <- as.character(ident)
    
    if (is.na(ident) || !nzchar(ident))
        return(character())
    
    # CSI:FingerID/PubChemLite might return multiple identifiers, separated by ; or a space
    # set consensus results can also merge multiple identifiers
    idlist <- unlist(strsplit(ident, ";| "))
    
    if (grepl("pubchem", tolower(db)))
        fmt <- "<a target=\"_blank\" href=\"https://pubchem.ncbi.nlm.nih.gov/compound/%s\">%s</a>"
    else if (tolower(db) == "chemspider")
        fmt <- "<a target=\"_blank\" href=\"http://www.chemspider.com/Search.aspx?q=%s\">%s</a>"
    else if (startsWith(idlist[1], "DTX"))
        fmt <- "<a target=\"_blank\" href=\"https://comptox.epa.gov/dashboard/dsstoxdb/results?search=%s\">%s</a>"
    else
        fmt <- "%s"
    
    return(paste0(sprintf(fmt, idlist, idlist), collapse = "; "))
}

getCompInfoList <- function(compResults, compIndex, mConsNames, addHTMLURL)
{
    columns <- names(compResults)

    resultRow <- compResults[compIndex, ]

    addValText <- function(curText, fmt, cols)
    {
        cols <- getAllMergedConsCols(cols, columns, mConsNames)
        ret <- character()
        for (cl in cols)
        {
            if (!is.null(resultRow[[cl]]) && !is.na(resultRow[[cl]]) &&
                (!is.character(resultRow[[cl]]) || nzchar(resultRow[[cl]])))
            {
                fm <- sprintf("%s: %s", cl, fmt)
                ret <- c(ret, sprintf(fm, resultRow[[cl]]))
            }
        }

        return(c(curText, ret))
    }

    ctext <- character()

    if (addHTMLURL)
    {
        addIdURL <- function(param, ident, db) return(sprintf("%s: %s", param, makeDBIdentLink(db, ident)))

        dbcols <- getAllMergedConsCols("database", columns, mConsNames)
        
        if (!is.null(resultRow[["identifier"]])) # compounds were not merged, can use 'regular' column
            ctext <- c(ctext, addIdURL("identifier", resultRow$identifier, resultRow$database))
        else
        {
            idcols <- getAllMergedConsCols("identifier", columns, mConsNames)

            if (allSame(resultRow[, idcols, with = FALSE])) # no need to show double ids
                ctext <- c(ctext, addIdURL("identifier", resultRow[[idcols[1]]], resultRow[[dbcols[1]]]))
            else
            {
                for (i in seq_along(idcols))
                    ctext <- c(ctext, addIdURL(idcols[i], resultRow[[idcols[i]]], resultRow[[dbcols[i]]]))
            }
        }
        
        relatedIDCols <- getAllMergedConsCols("relatedCIDs", columns, mConsNames)
        for (i in seq_along(relatedIDCols))
            ctext <- c(ctext, addIdURL(relatedIDCols[i], resultRow[[relatedIDCols[i]]], resultRow[[dbcols[i]]]))
    }
    else
    {
        ctext <- addValText(ctext, "%s", "identifier")
        ctext <- addValText(ctext, "%s", "relatedCIDs")
    }

    ctext <- addValText(ctext, "%s", c("compoundName", "neutral_formula", "SMILES"))

    if (length(getAllMergedConsCols("InChIKey", columns, mConsNames)) > 0)
        ctext <- addValText(ctext, "%s", "InChIKey")
    else # only add InChIKey1/2 if full isn't available
        ctext <- addValText(ctext, "%s", c("InChIKey1", "InChIKey2"))

    ctext <- addValText(ctext, "%.2f", c("XlogP", "AlogP"))

    # PubChemLite
    ctext <- addValText(ctext, "%s", c("FP", "compoundName2"))
    
    # Dashboard
    ctext <- addValText(ctext, "%s", c("CASRN", "QCLevel"))

    # FOR-IDENT
    ctext <- addValText(ctext, "%s", c("tonnage", "categories"))

    # TP prediction DB
    ctext <- addValText(ctext, "%s", c("precursor", "transformation", "enzymes"))

    return(ctext)
}

buildMFLandingURL <- function(mfSettings, peakList, precursorMz)
{
    # Via personal communication Steffen/Emma, see https://github.com/Treutler/MetFamily/blob/22b9f46b2716b805c24c03d260045605c0da8b3e/ClusteringMS2SpectraGUI.R#L2433
    # Code adopted from MetFamily R package: https://github.com/Treutler/MetFamily

    if (is.null(mfSettings))
    {
        # no settings given, simply default to PubChem
        mfSettings <- list(MetFragDatabaseType = "PubChem")
    }

    mfSettings$IonizedPrecursorMass <- precursorMz
    mfSettings$NeutralPrecursorMass <- "" # make sure user needs to calculate it and remove default

    PL <- paste(peakList$mz, peakList$intensity, sep = " ", collapse = "; ")
    mfSettings$PeakList <- PL

    if (mfSettings$MetFragDatabaseType == "ExtendedPubChem")
        mfSettings$MetFragDatabaseType <- "PubChem" # user should tick box for now...
    else if (!mfSettings$MetFragDatabaseType %in% c("KEGG", "PubChem", "ChemSpider", "LipidMaps", "MetaCyc", "LocalInChI", "LocalSDF"))
        mfSettings$MetFragDatabaseType <- NULL # not all databases are supported yet.

    # Allowed parameters: list taken from error page when unsupported parameter is given
    mfSettings <- mfSettings[names(mfSettings) %in%
                                 c("FragmentPeakMatchAbsoluteMassDeviation", "FragmentPeakMatchRelativeMassDeviation",
                                   "DatabaseSearchRelativeMassDeviation", "PrecursorCompoundIDs", "IonizedPrecursorMass",
                                   "NeutralPrecursorMass", "NeutralPrecursorMolecularFormula", "PrecursorIonMode",
                                   "PeakList", "MetFragDatabaseType")]

    setstr <- paste0(paste0(names(mfSettings), "=", mfSettings), collapse = "&")
    ret <- paste0("https://msbi.ipb-halle.de/MetFragBeta/landing.xhtml?", setstr)
    #ret <- sprintf("<a target=\"_blank\" href=\"%s\">MetFragWeb</a>", ret)

    return(ret)
}

getCompoundsSpecPlotTitle <- function(compoundName, formula, compoundName2 = NULL, formula2 = NULL)
{
    hasCName <- !is.null(compoundName) && !is.na(compoundName) && nzchar(compoundName)
    hasCName2 <- !is.null(compoundName2) && !is.na(compoundName2) && nzchar(compoundName2)
    if (hasCName && hasCName2)
        compoundName <- paste0(compoundName, "/", compoundName2)
    return(subscriptFormula(formula, over = if (hasCName) compoundName else NULL, formulas2 = formula2))
}
