# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

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
        fiLeft <- merge(fiLeft, fiRight[, c("PLID", "mergedBy", "ion_formula", "neutral_loss"), with = FALSE], all.x = TRUE, by = "PLID")
        fiLeft[is.na(mergedBy.y), mergedBy := mergedBy.x]
        fiLeft[is.na(mergedBy.x), mergedBy := mergedBy.y]
        fiLeft[!is.na(mergedBy.x) & !is.na(mergedBy.y), mergedBy := paste(mergedBy.x, mergedBy.y, sep = ",")]
        fiLeft[, c("mergedBy.x", "mergedBy.y") := NULL]
        fiLeft[, ion_formula := fifelse(!is.na(ion_formula.x), ion_formula.x, ion_formula.y)]
        fiLeft[, neutral_loss := fifelse(!is.na(neutral_loss.x), neutral_loss.x, neutral_loss.y)]
        fiLeft[, c("ion_formula.x", "ion_formula.y", "neutral_loss.x", "neutral_loss.y") := NULL]

        # add unique
        fiUnique <- fiRight[!PLID %in% fiLeft$PLID]
        if (nrow(fiUnique) > 0)
            fiLeft <- rbind(fiLeft, fiUnique, fill = TRUE)
        
        setorderv(fiLeft, "PLID")
    }

    return(fiLeft)
}

#' Scorings terms for compound candidates
#'
#' Returns an overview of scorings may be applied to rank candidate compounds.
#'
#' @param algorithm The algorithm: \code{"metfrag"} or \code{"sirius"}. Set to \code{NULL} to return all scorings.
#' @param database The database for which results should be returned (\emph{e.g.} \code{"pubchem"}). Set to \code{NULL}
#'   to return all scorings.
#' @param includeSuspectLists,onlyDefault,includeNoDB A logical specifying whether scoring terms related to suspect
#'   lists, default scoring terms and non-database specific scoring terms should be included in the output,
#'   respectively.
#'
#' @return A \code{data.frame} with information on which scoring terms are used, what their algorithm specific name is
#'   and other information such as to which database they apply and short remarks.
#'
#' @seealso generateCompounds
#'
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

addCompoundScore <- function(compounds, scoreName, updateScore, scoreWeight)
{
    if (updateScore)
    {
        compounds@groupAnnotations <- Map(groupNames(compounds), annotations(compounds), f = function(grp, ann)
        {
            if (all(is.na(ann[[scoreName]])) || max(ann[[scoreName]], na.rm = TRUE) == 0)
                return(ann)
            ann <- copy(ann)
            norm <- ann[[scoreName]] / max(ann[[scoreName]])
            ann[, score := score + (scoreWeight * norm)]
            return(ann)
        })
    }
    
    compounds@scoreRanges <- Map(compounds@scoreRanges, annotations(compounds), f = function(sc, ann)
    {
        rn <- if (all(is.na(ann[[scoreName]]))) c(NA_real_, NA_real_) else range(ann[[scoreName]], na.rm = TRUE)
        ret <- c(sc, setNames(list(rn), scoreName))
        # extend score range if necessary
        if (updateScore)
            ret$score <- c(min(ret$score, ann$score, na.rm = TRUE), max(ret$score, ann$score, na.rm = TRUE))
        return(ret)
    })
    compounds@scoreTypes <- union(compounds@scoreTypes, scoreName)
    
    return(compounds)
}

makeDBIdentLink <- function(db, ident)
{
    ident <- as.character(ident)
    
    if (length(ident) == 0)
        return(character())
    if (is.na(db) || is.na(ident))
        return("NA")
    
    # CSI:FingerID/PubChemLite might return multiple identifiers, separated by ; or a space
    # set consensus results can also merge multiple identifiers
    idlist <- strsplit(ident, ";| ")
    
    if (grepl("pubchem", tolower(db)))
        fmt <- "<a target=\"_blank\" href=\"https://pubchem.ncbi.nlm.nih.gov/compound/{ id }\">{ id }</a>"
    else if (tolower(db) == "chemspider")
        fmt <- "<a target=\"_blank\" href=\"http://www.chemspider.com/Search.aspx?q={ id }\">{ id }</a>"
    else if (startsWith(idlist[[1]][1], "DTX"))
        fmt <- "<a target=\"_blank\" href=\"https://comptox.epa.gov/dashboard/dsstoxdb/results?search={ id }\">{ id }</a>"
    else if (tolower(db) == "library")
        fmt <- paste0("{ id } (",
                      "<a target=\"_blank\" href=\"https://massbank.eu/MassBank/RecordDisplay?id={ id }\">MB.eu</a>",
                      " | ",
                      "<a target=\"_blank\" href=\"https://mona.fiehnlab.ucdavis.edu/spectra/display/{ id }\">MoNA</a>)")
    else
        fmt <- NULL

    ret <- character(length(ident))
    NAIdent <- is.na(ident) | !nzchar(ident)
    ret[NAIdent] <- NA_character_
    ret[!NAIdent] <- sapply(idlist[!NAIdent], function(id) {
        if (!is.null(fmt))
            id <- glue::glue(fmt, id = id)
        return(paste0(id, collapse = "; "))
    })
    
    return(ret)
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
    
    ctext <- addValText(ctext, "%s", "molNeutralized")

    ctext <- addValText(ctext, "%.2f", c("XlogP", "AlogP", "LogP"))

    # PubChemLite
    ctext <- addValText(ctext, "%s", c("FP", "compoundName2"))
    
    # Dashboard
    ctext <- addValText(ctext, "%s", c("CASRN", "QCLevel"))

    # FOR-IDENT
    ctext <- addValText(ctext, "%s", c("tonnage", "categories"))

    # TP prediction DB
    ctext <- addValText(ctext, "%s", c("parent", "transformation", "enzyme", "evidencedoi"))

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
