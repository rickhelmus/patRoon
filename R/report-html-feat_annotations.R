# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include report-html.R
NULL

makeAnnTitleTab <- function(title, tab)
{
    return(htmltools::div(
        style = list(margin = "10px 20px"),
        htmltools::div(style = list("text-align" = "center", "font-weight" = "bold"), title),
        tab
    ))
}

makeAnnKable <- function(tab, mark = NULL, ...)
{
    withr::local_options(list(knitr.kable.NA = ""))
    kab <- knitr::kable(tab, format = "html", escape = FALSE, ...) %>%
        kableExtra::kable_styling(font_size = 13, bootstrap_options = c("striped", "condensed"))
    # BUG: row_spec needs to be before scroll_box
    if (!is.null(mark))
        kab <- kableExtra::row_spec(kab, mark, italic = TRUE)
    kab <- kableExtra::scroll_box(kab, box_css = "border: 1px solid #ddd;",
                                  extra_css = "overflow-x: hidden; overflow-y: auto; height: 200px;")
    return(htmltools::span(dangerouslySetInnerHTML = list("__html" = kab)))
}

makeAnnDetailsReact <- function(title, tab, sets)
{
    if (nrow(tab) > 0)
    {
        tab <- copy(tab)
        cols <- copy(names(tab))
        for (col in cols)
        {
            if (is.na(tab[[col]]) | (is.character(tab[[col]]) & !nzchar(tab[[col]])))
                set(tab, j = col, value = NULL)
        }
    }
    ptab <- makePropTab(tab, sets, FALSE)
    return(makeAnnTitleTab(title, makeAnnKable(ptab, col.names = sub("common", "", names(ptab), fixed = TRUE))))
}

makeAnnPLReact <- function(apl)
{
    rtab <- if (is.null(apl) || nrow(apl) == 0)
        htmltools::div(align = "center", "No annotation available.")
    else
    {
        isPrec <- apl$precursor
        
        apl[, c("ID", "annotated", "precursor") := NULL]
        apl[, c("mz", "intensity") := .(round(mz, 5), round(intensity))]
        apl[, ion_formula := subscriptFormulaHTML(ion_formula)]
        apl[, neutral_loss := subscriptFormulaHTML(neutral_loss)]
        if (!is.null(apl[["ion_formula_MF"]]))
            apl[, ion_formula_MF := subscriptFormulaHTML(ion_formula_MF)]
        
        colDefs <- pruneList(list(
            ion_formula = reactable::colDef(html = TRUE, minWidth = 125),
            neutral_loss = reactable::colDef(html = TRUE, minWidth = 125),
            ion_formula_MF = if (!is.null(apl[["ion_formula_MF"]])) reactable::colDef(html = TRUE, minWidth = 150) else NULL
        ))
        
        makeAnnKable(apl, mark = if (any(isPrec)) which(isPrec) else NULL)
    }
    
    return(makeAnnTitleTab("Peak list annotations", rtab))
}

makeAnnScoreReact <- function(annRow, sets)
{
    annRow[, (names(annRow)) := lapply(.SD, round, 2), .SDcols = names(annRow)]
    ptab <- makePropTab(annRow, sets, FALSE)
    return(makeAnnTitleTab("Scorings", makeAnnKable(ptab)))
}

getAnnReactImgCell <- function(value) htmltools::img(src = value, style = list("max-height" = "300px"))

makeAnnReactable <- function(tab, id, detailsTabFunc = NULL, annPLTabFunc = NULL, scoreTabFunc = NULL, ...)
{
    details <- if (is.null(detailsTabFunc))
        NULL
    else
    {
        function(index)
        {
            htmltools::div(style = list(margin = "12px 45px", display = "flex", "flex-wrap" = "no-wrap",
                                        background = "#FCFCFC", border = "dashed 1px",
                                        "justify-content" = "space-between", "overflow-x" = "auto"),
                           detailsTabFunc(index), annPLTabFunc(index), scoreTabFunc(index))
        }
    }
    
    return(makeReactable(tab, id, compact = TRUE, details = details, pagination = TRUE,
                         language = reactable::reactableLang(noData = "No annotations available"), ...))
}

genHTMLReportPlotsStructs <- function(fGroups, compounds, settings, outPath, parallel)
{
    scrStructInfo <- if (isScreening(fGroups)) screenInfo(fGroups)[, c("SMILES", "InChIKey"), with = FALSE] else NULL
    compStructInfo <- NULL
    if (!is.null(compounds) && length(compounds) != 0)
    {
        compStructInfo <- subsetDTColumnsIfPresent(as.data.table(compounds), c("group", "SMILES", "InChIKey",
                                                                               "InChIKey1"))
        if (is.null(compStructInfo[["InChIKey"]])) # not present with eg SIRIUS
        {
            # fall back to first block. UNDONE: calculate InChIKey for SIRIUS results? Needs OpenBabel...
            setnames(compStructInfo, "InChIKey1", "InChIKey")
        }

        compStructInfo[, index := seq_len(.N), by = "group"]
        compStructInfo <- compStructInfo[index <= settings$compounds$topMost]
        compStructInfo <- removeDTColumnsIfPresent(compStructInfo, c("group", "index", "InChIKey1"))
    }
    
    structInfo <- rbindlist(list(scrStructInfo, compStructInfo))
    if (nrow(structInfo) > 0 && any(!is.na(structInfo$SMILES)))
    {
        structInfo <- unique(structInfo[!is.na(SMILES)], by = "InChIKey")
        cat("Generate structures...\n")
        return(doApply("Map", parallel, structInfo$InChIKey, structInfo$SMILES, f = function(ik, smi)
        {
            # NOTE: we use the InChIKey here instead of makeHash()
            pf <- file.path(getHTMLReportPlotPath(outPath), paste0("struct-", ik, ".svg"))
            if (!file.exists(pf))
                saveRCDKStructure(getMoleculesFromSMILES(smi)[[1]], "svg", pf, 100, 100)
            doProgress()
            return(pf)
        }))
    }
    return(list())
}

genHTMLReportPlotsMSPeakLists <- function(MSPeakLists, settings, outPath, parallel)
{
    if (!settings$MSPeakLists$spectra)
        return(list())
    
    cat("Generate MS spectra...\n")
    
    if (length(MSPeakLists) == 0)
        return(list())
    
    return(doApply("sapply", parallel, groupNames(MSPeakLists), function(grp)
    {
        ret <- list()
        
        args <- list(MSPeakLists, groupName = grp, title = "")
        pp <- list(mar = c(4.1, 4.1, 0.2, 0.2))
        
        ret$MS <- makeHTMLReportPlot("spec-MS", outPath, "plotSpectrum", c(list(MSLevel = 1), args), parParams = pp,
                                     width = 7, height = 4)
        
        ret$MSMS <- if (!is.null(MSPeakLists[[grp]][["MSMS"]]))
        {
            makeHTMLReportPlot("spec-MSMS", outPath, "plotSpectrum", c(list(MSLevel = 2), args), parParams = pp,
                               width = 7, height = 4)
        }
        else
            ""
        
        doProgress()
        
        return(ret)
    }, simplify = FALSE))
}

genHTMLReportPlotsFormulas <- function(formulas, MSPeakLists, settings, outPath, parallel)
{
    if (!settings$formulas$include)
        return(list())
    
    cat("Generate formula annotation plots...\n")
    
    if (length(formulas) == 0)
        return(list())
    
    return(doApply("Map", parallel, groupNames(formulas), annotations(formulas), f = function(grp, ann)
    {
        ret <- list()
        
        if (nrow(ann) > settings$formulas$topMost)
            ann <- ann[seq_len(settings$formulas$topMost)]
        
        ret$spectra <- sapply(seq_len(nrow(ann)), function(index)
        {
            if (is.null(MSPeakLists[[grp]][["MSMS"]]))
                return("")
            makeHTMLReportPlot("form-spec", outPath, "plotSpectrum", list(formulas, index, grp,
                                                                          MSPeakLists = MSPeakLists),
                               width = 7, height = 4)
        })
        
        ret$scores <- sapply(seq_len(nrow(ann)), function(index)
        {
            makeHTMLReportPlot("form-scores", outPath, "plotScores",
                               list(formulas, index, grp, normalizeScores = settings$formulas$normalizeScores,
                                    excludeNormScores = settings$formulas$exclNormScores), width = 6, height = 5)
        })
        
        doProgress()
        
        return(ret)
    }))
}

genHTMLReportPlotsCompounds <- function(compounds, MSPeakLists, formulas, settings, outPath, parallel)
{
    cat("Generate compound annotation plots...\n")
    
    if (length(compounds) == 0)
        return(list())
    
    return(doApply("Map", parallel, groupNames(compounds), annotations(compounds), f = function(grp, ann)
    {
        ret <- list()
        
        if (nrow(ann) > settings$compounds$topMost)
            ann <- ann[seq_len(settings$compounds$topMost)]
        
        ret$spectra <- sapply(seq_len(nrow(ann)), function(index)
        {
            if (is.null(MSPeakLists[[grp]][["MSMS"]]))
                return("")
            makeHTMLReportPlot("comp-spec", outPath, "plotSpectrum",
                               list(compounds, index, grp, MSPeakLists, formulas, FALSE, title = ""),
                               parParams = list(mar = c(4.1, 4.1, 0.2, 0.2)), width = 7, height = 4, pointsize = 16)
        })
        
        ret$scores <- sapply(seq_len(nrow(ann)), function(index)
        {
            makeHTMLReportPlot("comp-scores", outPath, "plotScores",
                               list(compounds, index, grp, normalizeScores = settings$compounds$normalizeScores,
                                    excludeNormScores = settings$compounds$exclNormScores),
                               parParams = list(mar = c(4.1, 4.1, 0.4, 0.2)), width = 7, height = 4, pointsize = 16)
        })
        
        doProgress()
        
        return(ret)
    }))
}

genHTMLReportPlotsCompsCluster <- function(compsCluster, settings, outPath, parallel)
{
    cat("Generate compound cluster plots...\n")
    
    if (length(compsCluster) == 0)
        return(list())
    
    return(doApply("Map", parallel, groupNames(compsCluster), cutClusters(compsCluster), f = function(grp, ct)
    {
        ret <- list()
        
        ret$dendro <- makeHTMLReportPlot("comp-clust-dendro", outPath, "plot", list(compsCluster, groupName = grp),
                                         width = 12, height = 4, pointsize = 16)
        
        ret$mcs <- sapply(sort(unique(ct)), function(cli)
        {
            makeHTMLReportPlot("comp-clust-mcs", outPath, "plotStructure",
                               list(compsCluster, groupName = grp, cluster = cli, 100, 100), width = 5, height = 4)
        })
        
        doProgress()
        
        return(ret)
    }))
}


reportHTMLUtils$methods(
    genMSPLTable = function(MSLevel)
    {
        if (MSLevel == 1)
            mdprintf("MS peak lists... ")
        
        MSPeakLists <- objects$MSPeakLists[, names(objects$fGroups)]
        
        id <- if (MSLevel == 2) "MSMSPLTab" else "MSPLTab"
        theight <- 350
        
        if (length(MSPeakLists) == 0)
        {
            # dummy table to show empty results
            return(makeAnnReactable(data.table(ID = integer()), id, height = theight))
        }
        
        tab <- as.data.table(MSPeakLists)
        tab <- tab[type == if (MSLevel == 2) "MSMS" else "MS"][, type := NULL]
        isPrec <- tab$precursor
        tab[, precursor := NULL]
        
        colDefs <- list(
            group = reactable::colDef(show = FALSE, filterMethod = reactExactFilter()),
            mz = reactable::colDef(format = reactable::colFormat(digits = 5)),
            intensity = reactable::colDef(format = reactable::colFormat(digits = 0))
        )
        
        colDefs <- setReactNumRangeFilters(id, tab, colDefs)
        if (!is.null(tab[["set"]]))
            colDefs$set <- reactable::colDef(filterInput = function(values, name) reactSelectFilter(id, values, name))
        
        return(makeAnnReactable(tab, id, columns = colDefs, height = theight, filterable = TRUE,
                                rowClass = function(index) if (isPrec[index]) "fw-bold" else ""))
    },
    
    genFormulasTable = function()
    {
        mdprintf("Formulas... ")
        
        formulas <- objects$formulas[names(objects$fGroups)]
        
        mcn <- mergedConsensusNames(formulas)
        
        if (length(formulas) == 0)
        {
            # dummy table to show empty results
            return(makeAnnReactable(data.table(formula = character()), "formulasTab"))
        }
        
        hash <- makeHash(formulas, objects$MSPeakLists[, groupNames(formulas)], settings$formulas)
        # set fixDTs to FALSE, since it considerably slows down cache loading and we're not dealing with data.tables here
        cd <- loadCacheData("reportHTMLFormulas", hash, fixDTs = FALSE)
        if (!is.null(cd))
            return(cd)
        
        tab <- as.data.table(formulas)
        
        # NOTE: for consensus results, duplicate algo columns (eg explainedPeaks) are only shown in details
        tab <- subsetDTColumnsIfPresent(tab, c("group", "neutral_formula", "neutralMass", "explainedPeaks",
                                               "explainedIntensity", "error"))
        
        tab[, candidate := seq_len(.N), by = "group"]
        
        tab <- tab[candidate <= settings$formulas$topMost]
        
        for (col in names(tab))
        {
            if (is.numeric(tab[[col]]))
                tab[, (col) := round(get(col), if (col == "neutralMass") 5 else 2)]
        }
        
        tab[, spectrum := plots$formulas[[group]]$spectra, by = "group"]
        tab[, scorings := plots$formulas[[group]]$scores, by = "group"]
        
        if (hasSuspects())
        {
            scr <- screenInfo(objects$fGroups)
            if (!is.null(scr[["formula"]]))
            {
                wh <- which(tab$neutral_formula %chin% scr$formula)
                if (length(wh) == 0)
                    tab[, suspect := ""]
                else
                {
                    tab[wh, suspect := {
                        paste0(unique(scr[formula == neutral_formula]$name), collapse = ", ")
                    }, by = "neutral_formula"]
                }
            }
            else
                tab[, suspect := ""]
        }
        
        getFormDetails <- function(index)
        {
            ft <- formulas[[tab$group[index]]][tab$candidate[index]]
            ft <- ft[, setdiff(names(ft), names(tab)), with = FALSE]
            takeCols <- getAllMergedConsCols(c("neutral_formula", "ion_formula", "neutralMass", "ion_formula_mz",
                                               "error", "error_frag_median", "error_frag_median_abs",
                                               "explainedPeaks", "explainedIntensity"), names(ft), mcn)
            ft <- ft[, takeCols, with = FALSE]
            for (col in getAllMergedConsCols(c("neutral_formula", "ion_formula"), names(ft), mcn))
                set(ft, j = col, value = subscriptFormulaHTML(ft[[col]]))
            massCols <- getAllMergedConsCols(c("neutralMass", "ion_formula_mz"), names(ft), mcn)
            numCols <- names(ft)[sapply(ft, is.numeric)]
            for (col in numCols)
                set(ft, j = col, value = round(ft[[col]], if (col %chin% massCols) 5 else 2))
            return(makeAnnDetailsReact("Formula properties", ft,
                                       if (isFGSet(objects$fGroups)) sets(objects$fGroups) else NULL))
        }
        
        getAnnPLDetails <- function(index)
        {
            apl <- annotatedPeakList(formulas, index = tab$candidate[index], groupName = tab$group[index],
                                     MSPeakLists = objects$MSPeakLists, onlyAnnotated = TRUE)
            return(makeAnnPLReact(apl))
        }
        
        getScoreDetails <- function(index)
        {
            fRow <- formulas[[tab$group[index]]][tab$candidate[index]]
            cols <- getAllMergedConsCols(annScoreNames(formulas, FALSE), names(fRow), mcn)
            return(makeAnnScoreReact(fRow[, cols, with = FALSE],
                                     if (isFGSet(objects$fGroups)) sets(objects$fGroups) else NULL))
        }
        
        setcolorder(tab, intersect(c("candidate", "suspect"), names(tab)))
        
        colDefs <- pruneList(list(
            group = reactable::colDef(show = FALSE, filterMethod = reactExactFilter()),
            candidate = reactable::colDef("#", minWidth = 15),
            suspect = if (!is.null(tab[["suspect"]])) reactable::colDef("suspect(s)", filterMethod = reactSuspectFilter()) else NULL,
            neutral_formula = reactable::colDef("formula",
                                                cell = function(value) htmltools::span(dangerouslySetInnerHTML = list("__html" = subscriptFormulaHTML(value)))),
            neutralMass = reactable::colDef("neutral mass"),
            spectrum = reactable::colDef(cell = getAnnReactImgCell, minWidth = 200),
            scorings = reactable::colDef(cell = getAnnReactImgCell, minWidth = 200)
        ))
        
        colDefs <- setReactNumRangeFilters("formulasTab", tab, colDefs)
        
        CSVCols <- setdiff(names(tab), c("spectrum", "scorings"))
        ret <- makeAnnReactable(tab, "formulasTab", columns = colDefs, getFormDetails, getAnnPLDetails,
                                getScoreDetails, meta = list(CSVCols = CSVCols))
        saveCacheData("reportHTMLFormulas", ret, hash)
        return(ret)
    },
    
    genCompoundsTable = function()
    {
        mdprintf("Compounds... ")
        
        compounds <- objects$compounds[names(objects$fGroups)]
        
        mcn <- mergedConsensusNames(compounds)
        
        if (length(compounds) == 0)
        {
            # dummy table to show empty results
            return(makeAnnReactable(data.table(compound = character()), "compoundsTab"))
        }
        
        hash <- makeHash(compounds, objects$formulas[groupNames(compounds)], objects$MSPeakLists[, groupNames(compounds)],
                         settings$compounds)
        cd <- loadCacheData("reportHTMLCompounds", hash, fixDTs = FALSE)
        if (!is.null(cd))
            return(cd)
        
        tab <- as.data.table(compounds)
        
        if (!is.null(tab[["database"]]))
            databases <- tab$database
        
        # NOTE: for consensus results, duplicate algo columns (eg identifier) are only shown in details
        tab <- subsetDTColumnsIfPresent(tab, c("group", "compoundName", "compoundName2", "identifier",
                                               "neutral_formula", "neutralMass", "explainedPeaks", "score", "InChIKey",
                                               "UID"))
        
        tab[, candidate := seq_len(.N), by = "group"]
        
        tab <- tab[candidate <= settings$compounds$topMost]
        
        cmpNames2 <- tab[["compoundName2"]]
        if (!is.null(cmpNames2))
            tab[, compoundName2 := NULL]
        
        if (!is.null(tab[["neutralMass"]]))
            tab[, neutralMass := round(neutralMass, 5)]
        if (!is.null(tab[["score"]]))
            tab[, score := round(score, 2)]
        
        if (is.null(tab[["InChIKey"]])) # SIRIUS
            tab[, structure := plots$structs[UID]]
        else
            tab[, structure := plots$structs[InChIKey]]
        
        tab[, spectrum := plots$compounds[[group]]$spectra, by = "group"]
        tab[, scorings := plots$compounds[[group]]$scores, by = "group"]
        
        if (hasSuspects())
        {
            scr <- screenInfo(objects$fGroups)
            if (!is.null(scr[["InChIKey"]]) && any(!is.na(scr$InChIKey)))
            {
                scr <- data.table::copy(scr)
                scr[, IK1 := getIKBlock1(InChIKey)]
                wh <- which(tab$UID %chin% scr$IK1)
                if (length(wh) == 0)
                    tab[, suspect := ""]
                tab[wh, suspect := {
                    paste0(unique(scr[IK1 == UID[1]]$name), collapse = ", ")
                }, by = "UID"]
            }
            else
                tab[, suspect := ""]
        }
        
        tab <- removeDTColumnsIfPresent(tab, c("InChIKey", "UID"))
        
        getCompCell <- function(value, index)
        {
            cn2 <- cmpNames2[index]
            if (!is.null(cn2) && !is.na(cn2))
            {
                if (is.na(value))
                    value <- cn2
                value <- htmltools::tagList(htmltools::strong(value), htmltools::br(), paste0("(", cn2, ")"))
            }
            else
                value <- htmltools::strong(value)
            return(value)
        }
        
        getCompDetails <- function(index)
        {
            ct <- data.table::copy(compounds[[tab$group[index]]][tab$candidate[index]])
            ct <- ct[, setdiff(names(ct), names(tab)), with = FALSE]
            takeCols <- c(
                "compoundName",
                "identifier",
                "relatedCIDs",
                "database",
                "explainedPeaks",
                "neutral_formula",
                "SMILES",
                "InChI",
                "InChIKey",
                "XlogP", "AlogP", "LogP",
                
                # PubChemLite
                "FP", "compoundName2", # UNDONE: update?
                
                # CompTox
                "CASRN", "QCLevel",
                
                # FOR-IDENT
                "tonnage", "categories",
                
                # TP DB
                "parent", "transformation", "enzyme", "evidencedoi"
            )
            ct <- ct[, getAllMergedConsCols(takeCols, names(ct), mcn), with = FALSE]
            
            if (!is.null(tab[["compoundName"]]) && !is.null(tab[["compoundName2"]]))
                set(ct, j = "compoundName2", value = NULL) # already in main table compoundName column
            
            set(ct, j = "neutral_formula", value = subscriptFormulaHTML(ct$neutral_formula))
            
            numCols <- names(ct)[sapply(ct, is.numeric)]
            for (col in numCols)
                set(ct, j = col, value = round(ct[[col]], 2))
            
            for (col in getAllMergedConsCols(c("identifier", "relatedCIDs"), names(ct), mcn))
            {
                # get db: handle consensus results
                cn <- sub("identifier|relatedCIDs", "", col)
                db <- if (nzchar(cn)) ct[[paste0("database", cn)]] else ct$database
                set(ct, j = col, value = makeDBIdentLink(db, ct[[col]]))
            }
            
            return(makeAnnDetailsReact("Compound properties", ct,
                                       if (isFGSet(objects$fGroups)) sets(objects$fGroups) else NULL))
        }
        
        getAnnPLDetails <- function(index)
        {
            apl <- annotatedPeakList(compounds, index = tab$candidate[index], groupName = tab$group[index],
                                     MSPeakLists = objects$MSPeakLists, formulas = objects$formulas,
                                     onlyAnnotated = TRUE)
            return(makeAnnPLReact(apl))
        }
        
        getScoreDetails <- function(index)
        {
            cRow <- compounds[[tab$group[index]]][tab$candidate[index]]
            cols <- getAllMergedConsCols(annScoreNames(compounds, FALSE), names(cRow), mcn)
            return(makeAnnScoreReact(cRow[, cols, with = FALSE],
                                     if (isFGSet(objects$fGroups)) sets(objects$fGroups) else NULL))
        }
        
        setcolorder(tab, intersect(c("candidate", "compoundName", "suspect", "structure"), names(tab)))
        
        colDefs <- pruneList(list(
            group = reactable::colDef(show = FALSE, filterMethod = reactExactFilter()),
            candidate = reactable::colDef("#", minWidth = 15),
            compoundName = if (!is.null(tab[["compoundName"]])) reactable::colDef("compound", cell = getCompCell) else NULL,
            suspect = if (!is.null(tab[["suspect"]])) reactable::colDef("suspect(s)", filterMethod = reactSuspectFilter()) else NULL,
            identifier = if (!is.null(tab[["identifier"]])) reactable::colDef(cell = function(value, index) htmltools::span(dangerouslySetInnerHTML = list("__html" = makeDBIdentLink(databases[index], value)))) else NULL,
            neutral_formula = reactable::colDef("formula",
                                                cell = function(value) htmltools::span(dangerouslySetInnerHTML = list("__html" = subscriptFormulaHTML(value)))),
            neutralMass = if (!is.null(tab[["neutralMass"]])) reactable::colDef("neutral mass") else NULL,
            structure = reactable::colDef(cell = getAnnReactImgCell, minWidth = 125),
            spectrum = reactable::colDef(cell = getAnnReactImgCell, minWidth = 200),
            scorings = reactable::colDef(cell = getAnnReactImgCell, minWidth = 200)
        ))
        
        colDefs <- setReactNumRangeFilters("compoundsTab", tab, colDefs)
        
        mfWebLinks <- sapply(names(objects$fGroups), function(grp)
        {
            if (is(compounds, "compoundsMF") && grp %chin% names(compounds))
            {
                # make link with used MF settings
                set <- settings(compounds)
            }
            else
                set <- NULL
            
            return(buildMFLandingURL(set, objects$MSPeakLists[[grp]][["MSMS"]],
                                     groupInfo(objects$fGroups)[grp, "mzs"]))
        })
        
        CSVCols <- setdiff(names(tab), c("structure", "spectrum", "scorings"))
        
        ret <- makeAnnReactable(tab, "compoundsTab", columns = colDefs, getCompDetails, getAnnPLDetails,
                                getScoreDetails, meta = list(mfWebLinks = mfWebLinks, CSVCols = CSVCols))
        saveCacheData("reportHTMLCompounds", ret, hash)
        return(ret)
    },
    
    genCompClustsImgs = function()
    {
        elements <- unlist(Map(names(plots$compsCluster), plots$compsCluster, f = function(grp, pl)
        {
            lapply(seq_along(pl$mcs), function(i)
            {
                img <- htmltools::img(src = pl$mcs[[i]], class = paste0("mcs mcs-", grp),
                                      style = "max-height: 250px; display: none;")
            })
        }), recursive = FALSE)
        elements <- do.call(htmltools::tagList, elements)
        return(htmltools::div(style = list(display = "flex"), elements))
    }
)
