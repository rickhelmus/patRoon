#' @include main.R
#' @include report-html.R
NULL

getFGReactTab <- function(objects, settings, ...)
{
    tab <- as.data.table(objects$fGroups, qualities = "score", average = TRUE,
                         concAggrParams = settings$concAggrParams, toxAggrParams = settings$toxAggrParams, ...)
    
    if (settings$features$chromatograms$small)
        tab[, chrom_small := group]
    if (settings$features$chromatograms$large)
        tab[, chrom_large := group]
    
    if (!is.null(objects$MSPeakLists) || !is.null(objects$formulas) || !is.null(objects$compounds))
    {
        tab[, annotations := {
            grp <- group[1]
            hasMSMS <- !is.null(objects$MSPeakLists) && !is.null(objects$MSPeakLists[[grp]]) &&
                !is.null(objects$MSPeakLists[[grp]][["MSMS"]])
            hasForm <- !is.null(objects$formulas) && !is.null(objects$formulas[[grp]])
            hasComp <- !is.null(objects$compounds) && !is.null(objects$compounds[[grp]])
            ann <- c(if (hasMSMS) "MS2", if (hasForm) "Form", if (hasComp) "Comp")
            paste0(ann, collapse = ";")
        }, by = "group"]
    }

    if (isFGSet(objects$fGroups))
    {
        # collapse adduct and ion_mz columns
        combCols <- function(x) { x <- x[!is.na(x)]; return(if (length(x) == 0) "" else paste0(x, collapse = ";")) }
        tab[, adduct := combCols(unlist(.SD)), .SDcols = paste0("adduct-", sets(objects$fGroups)), by = seq_len(nrow(tab))]
        tab[, ion_mz := combCols(unlist(.SD)), .SDcols = paste0("ion_mz-", sets(objects$fGroups)), by = seq_len(nrow(tab))]
        tab[, (grep("^(adduct|ion_mz)\\-", names(tab), value = TRUE)) := NULL]
    }
    
    return(tab)
}

getFGScreeningReactTab <- function(tab, plots)
{
    tab <- copy(tab)
    scols <- setdiff(names(tab), "group")
    setnames(tab, scols, paste0("susp_", scols))
    if (!is.null(tab[["susp_InChIKey"]]) && !all(is.na(tab$susp_InChIKey)))
        tab[, susp_structure := plots$structs[getIKBlock1(susp_InChIKey)]]
    return(tab)
}


genHTMLReportPlotsChromsLarge <- function(fGroups, settings, outPath, EICs, EICParams, parallel)
{
    if (!settings$features$chromatograms$large)
        return(list())
    
    cat("Generate large chromatograms...\n")
    doApply("sapply", parallel, names(fGroups), function(grp)
    {
        doProgress()
        makeHTMLReportPlot("chrom_large-", outPath, "plotChroms",
                           list(fGroups, groupName = grp, retMin = settings$features$retMin,
                                intMax = settings$features$chromatograms$intMax, EICs = EICs,
                                EICParams = EICParams, groupBy = "rGroups", title = "", bty = "l"),
                           parParams = list(mar = c(4.1, 4.1, 0.2, 0.2)),
                           width = 6, height = 4, bg = "transparent", pointsize = 16)
    }, simplify = FALSE)
}

genHTMLReportPlotsChromsSmall <- function(fGroups, settings, outPath, EICs, EICParams, parallel)
{
    if (!settings$features$chromatograms$small)
        return(list())
    
    cat("Generate small chromatograms...\n")
    doApply("sapply", parallel, names(fGroups), function(grp)
    {
        doProgress()
        makeHTMLReportPlot("chrom_small", outPath, "plotChroms",
                           list(fGroups, groupName = grp, retMin = settings$features$retMin, EICs = EICs,
                                EICParams = modifyList(EICParams, list(topMost = 1, topMostByRGroup = FALSE,
                                                                       onlyPresent = TRUE)),
                                showFGroupRect = FALSE, showPeakArea = TRUE, title = "",
                                intMax = settings$features$chromatograms$intMax, bty = "n"),
                           parParams = list(mai = c(0, 0, 0, 0), lwd = 10), width = 12, height = 4, bg = "transparent",
                           pointsize = 16)
    }, simplify = FALSE)
}

genHTMLReportPlotsChromsFeatures <- function(fGroups, settings, outPath, EICs, EICParams, parallel)
{
    if (isFALSE(settings$features$chromatograms$features))
        return(list())
    
    anas <- analyses(fGroups)
    cat("Generate individual feature chromatograms...\n")
    doApply("sapply", parallel, names(fGroups), function(grp)
    {
        doProgress()
        Map(anas, seq_along(anas), f = function(ana, anai)
        {
            if (settings$features$chromatograms$features != "all" && fGroups[[grp]][anai] == 0)
                return("")
            makeHTMLReportPlot("chrom_feat", outPath, "plotChroms",
                               list(fGroups, analysis = ana, groupName = grp, retMin = settings$features$retMin,
                                    EICs = EICs, EICParams = modifyList(EICParams, list(topMost = NULL,
                                                                                        onlyPresent = settings$features$chromatograms$features != "all"),
                                                                        keep.null = TRUE), showFGroupRect = FALSE,
                                    showPeakArea = TRUE, title = "", intMax = settings$features$chromatograms$intMax,
                                    bty = "l"),
                               parParams = list(mar = c(4.1, 4.1, 0.2, 0.2)), width = 6, height = 4, bg = "transparent",
                               pointsize = 20, scaling = 1)
        })
    }, simplify = FALSE)
}

genHTMLReportPlotsIntPlots <- function(fGroups, settings, outPath, parallel)
{
    if (!settings$features$intensityPlots)
        return(list())
    
    cat("Generate intensity plots...\n")
    
    mainArgs <- list(average = TRUE, normalize = TRUE, plotArgs = list(bty = "l"))
    if (isFGSet(fGroups))
        mainArgs <- c(mainArgs, list(groupBy = "set", showLegend = TRUE))
    
    doApply("sapply", parallel, names(fGroups), function(grp)
    {
        doProgress()
        makeHTMLReportPlot("int_plot", outPath, "plotInt", c(list(fGroups[, grp]), mainArgs),
                           parParams = list(mar = c(4.1, 4.1, 1, 0.1)), width = 8, height = 4, bg = "transparent",
                           pointsize = 16)
    }, simplify = FALSE)
}


reportHTMLUtils$methods(
    genMainTablePlain = function()
    {
        mdprintf("Feature groups... ")
        makeMainResultsReactable(getFGReactTab(objects, settings), "Plain", settings$features$retMin, plots,
                                 initView = "Plain")
    },
    
    genMainTableSusByGroup = function()
    {
        makeMainResultsReactable(getFGReactTab(objects, settings, onlyHits = TRUE), "SusByGroup",
                                 settings$features$retMin, plots, initView = "SuspectsByGroup")
    },
    genMainTableSusCandSuspect = function()
    {
        makeMainResultsReactable(getFGScreeningReactTab(screenInfo(objects$fGroups), plots), "SusCandSuspect",
                                 settings$features$retMin, plots)
    },
    genMainTableSusBySuspect = function()
    {
        tab <- getFGScreeningReactTab(screenInfo(objects$fGroups), plots)
        tab[, susp_groups := paste0(group, collapse = ", "), by = "susp_name"][, group := NULL]
        tab <- unique(tab, by = "susp_name")
        makeMainResultsReactable(tab, "SusBySuspect", settings$features$retMin, plots, initView = "SuspectsBySuspect")
    },
    genMainTableSusCandGroup = function()
    {
        tab <- getFGReactTab(objects, settings, collapseSuspects = NULL, onlyHits = TRUE)
        # HACK: use a different name (and col definition) so that we get a hidden column used for filtering
        setnames(tab, "susp_name", "susp_ID")
        makeMainResultsReactable(tab, "SusCandGroup", settings$features$retMin, plots)
    },

    # HACK: the ISTD table is essentially the same as what screenInfo() returns for suspects. So the following functions
    # do largely the same.
    
    genMainTableISTDsByGroup = function()
    {
        tab <- getFGReactTab(objects, settings)
        tab <- tab[group %chin% internalStandards(objects$fGroups)$group]
        makeMainResultsReactable(tab, "ISTDsByGroup", settings$features$retMin, plots, initView = "ISTDsByGroup")
    },
    genMainTableISTDsCandISTD = function()
    {
        makeMainResultsReactable(getFGScreeningReactTab(internalStandards(objects$fGroups), plots), "ISTDsCandISTD",
                                 settings$features$retMin, plots)
    },
    genMainTableISTDsByISTD = function()
    {
        tab <- getFGScreeningReactTab(internalStandards(objects$fGroups), plots)
        # UNDONE: give other name
        tab[, susp_groups := paste0(group, collapse = ", "), by = "susp_name"][, group := NULL]
        tab <- unique(tab, by = "susp_name")
        makeMainResultsReactable(tab, "ISTDsByISTD", settings$features$retMin, plots, initView = "ISTDsByISTD")
    },
    genMainTableISTDsCandGroup = function()
    {
        tab <- getFGScreeningReactTab(internalStandards(objects$fGroups), plots)
        # HACK: use a different name (and col definition) so that we get a hidden column used for filtering
        setnames(tab, "susp_name", "susp_ID")

        ftab <- getFGReactTab(objects, settings)
        ftab <- ftab[group %chin% internalStandards(objects$fGroups)$group]
                
        ftab <- merge(tab, ftab, by = "group")

        makeMainResultsReactable(ftab, "ISTDsCandGroup", settings$features$retMin, plots)
    },

    genSuspInfoTable = function(id)
    {
        tab <- as.data.table(objects$fGroups, collapseSuspects = NULL)
        tab <- tab[, getAllSuspCols(paste0("susp_", suspMetaDataCols()), names(tab),
                                    mergedConsensusNames(objects$fGroups)), with = FALSE]
        
        setnames(tab, sub("^susp_", "", names(tab)))
        
        tab <- unique(tab, by = "name")
        
        for (col in intersect(c("neutralMass", "mz"), names(tab)))
            set(tab, j = col, value = round(tab[[col]], 5))
        if (!is.null(tab[["rt"]]))
        {
            if (settings$features$retMin)
                set(tab, j = "rt", value = tab$rt / 60)
            set(tab, j = "rt", value = round(tab$rt, 2))
        }
        if (!is.null(tab[["formula"]]))
            set(tab, j = "formula", value = subscriptFormulaHTML(tab$formula))
        
        ptab <- makePropTab(tab, NULL, "name")
        makePropReactable(ptab, id, "name", minPropWidth = 120, minValWidth = 150)
    },
    
    genISTDInfoTable = function()
    {
        tab <- data.table::copy(internalStandards(objects$fGroups))
        tab <- subsetDTColumnsIfPresent(tab, suspMetaDataCols())
        tab <- unique(tab, by = "name")
        
        for (col in intersect(c("neutralMass", "mz"), names(tab)))
            set(tab, j = col, value = round(tab[[col]], 5))
        if (!is.null(tab[["rt"]]))
        {
            if (settings$features$retMin)
                set(tab, j = "rt", value = tab$rt / 60)
            set(tab, j = "rt", value = round(tab$rt, 2))
        }
        if (!is.null(tab[["formula"]]))
            set(tab, j = "formula", value = subscriptFormulaHTML(tab$formula))
        
        ptab <- makePropTab(tab, NULL, "name")
        makePropReactable(ptab, "ISTDInfoTab", "name", minPropWidth = 120, minValWidth = 150)
    },
    
    genFeaturesTable = function()
    {
        mdprintf("Features... ")
        
        tab <- as.data.table(getFeatures(objects$fGroups))
        tab <- removeDTColumnsIfPresent(tab, "adduct") # can already be seen in group table
        tab <- removeDTColumnsIfPresent(tab, featureQualityNames(group = FALSE, scores = FALSE)) # only show scores
        
        if (settings$features$retMin)
            tab[, c("ret", "retmin", "retmax") := .(ret / 60, retmin / 60, retmax / 60)]
        
        for (col in names(tab)[sapply(tab, is.numeric)])
            set(tab, j = col, value = round(tab[[col]], if (col %in% c("mz", "mzmin", "mzmax", "ion_mz")) 5 else 2))
        
        anaInfo <- analysisInfo(objects$fGroups)
        
        if (settings$features$chromatograms$features == "all")
        {
            # add data for 'missing' analyses
            missingTab <- rbindlist(sapply(unique(tab$group), function(grp)
            {
                mt <- data.table(analysis = setdiff(analyses(objects$fGroups), tab[group == grp]$analysis))
                if (!is.null(tab[["set"]]))
                    mt[, set := anaInfo$set[match(analysis, anaInfo$analysis)]]
                return(mt)
            }, simplify = FALSE), idcol = "group")
            if (nrow(missingTab) > 0)
            {
                tab <- rbind(tab, missingTab, fill = TRUE)
                # make sure analyses retain order of anaInfo
                tab[, anaInd := match(analysis, anaInfo$analysis)]
                setorderv(tab, "anaInd")
                tab[, anaInd := NULL]
            }
        }
        
        tab[, rGroup := anaInfo[match(tab$analysis, analysis)]$group]

        colDefs <- list(
            group = reactable::colDef(show = FALSE, filterMethod = getReactFilterMethodExact()),
            rGroup = reactable::colDef("replicate group")
        )
        if (!isFALSE(settings$features$chromatograms$features))
        {
            # add EICs
            tab[, chromatogram := ""] # dummy value, not needed
            
            colDefs$chromatogram <- reactable::colDef(minWidth = 175, cell = function(value, index)
            {
                htmltools::img(src = plots$chromsFeatures[[tab$group[index]]][[tab$analysis[index]]])
            })
        }
        if (!is.null(tab[["set"]]))
            colDefs$set <- reactable::colDef(filterInput = makeReactFilterInputSelect("featuresTab"))
        
        fqn <- featureQualityNames(group = FALSE, scores = TRUE)
        for (col in fqn)
        {
            if (!is.null(tab[[col]]))
                colDefs[[col]] <- reactable::colDef(show = FALSE)
        }
        
        colDefs <- setReactNumRangeFilters("featuresTab", tab, colDefs)
        
        colDefs$analysis <- setReactSelRangeFilter("featuresTab", reactable::colDef())
        colDefs$rGroup <- setReactSelRangeFilter("featuresTab", colDefs$rGroup)
        
        setcolorder(tab, intersect(c("analysis", "rGroup", "ID", "chromatogram"), names(tab)))
        
        CSVCols <- setdiff(names(tab), "chromatogram")
        
        makeReactable(tab, "featuresTab", compact = TRUE, defaultExpanded = TRUE, columns = colDefs, filterable = FALSE,
                      meta = list(colToggles = list(qualities = fqn), CSVCols = CSVCols, internFilterable = "group",
                                  neverFilterable = "chromatogram"), pagination = TRUE)
    },
    
    genConcsTable = function()
    {
        concs <- data.table::copy(concentrations(objects$fGroups))
        
        imgTag <- function(IK1, cand) sprintf("<img src='%s' alt='%s' style='max-height: 300px;')></img>", plots$structs[IK1], cand)
        
        concs[type == "SIRIUS_FP", candidate := subscriptFormulaHTML(candidate)]
        concs[type == "suspect", candidate := imgTag(getIKBlock1(screenInfo(objects$fGroups)[match(candidate, SMILES)]$InChIKey),
                                                     candidate)]
        if (!is.null(objects[["compounds"]]))
        {
            cTab <- as.data.table(objects$compounds)
            concs[type == "compound", candidate := imgTag(cTab[match(candidate, SMILES)]$UID, candidate)]
        }
        
        colDefs <- list(
            group = reactable::colDef(show = FALSE),
            candidate = reactable::colDef(html = TRUE)
        )
        if (!is.null(concs[["candidate_name"]]))
            colDefs$candidate_name <- reactable::colDef("candidate name")
        colDefs[analyses(objects$fGroups)] <- list(reactable::colDef(format = reactable::colFormat(digits = 2)))
        colDefs$type <- setReactSelRangeFilter("concsTab", reactable::colDef())
        colDefs <- setReactNumRangeFilters("concsTab", concs, colDefs)
        
        makeReactable(concs, "concsTab", columns = colDefs, pagination = TRUE, filterable = FALSE,
                      meta = list(internFilterable = "group", neverFilterable = "group"))
    },

    genToxTable = function()
    {
        tox <- data.table::copy(toxicities(objects$fGroups))
        
        imgTag <- function(IK1, cand) sprintf("<img src='%s' alt='%s' style='max-height: 300px;')></img>", plots$structs[IK1], cand)
        
        tox[type == "SIRIUS_FP", candidate := subscriptFormulaHTML(candidate)]
        tox[type == "suspect", candidate := imgTag(getIKBlock1(screenInfo(objects$fGroups)[match(candidate, SMILES)]$InChIKey),
                                                     candidate)]
        if (!is.null(objects[["compounds"]]))
        {
            cTab <- as.data.table(objects$compounds)
            tox[type == "compound", candidate := imgTag(cTab[match(candidate, SMILES)]$UID, candidate)]
        }
        
        colDefs <- list(
            group = reactable::colDef(show = FALSE),
            candidate = reactable::colDef(html = TRUE),
            LC50 = reactable::colDef(format = reactable::colFormat(digits = 2))
        )
        if (!is.null(tox[["candidate_name"]]))
            colDefs$candidate_name <- reactable::colDef("candidate name")
        colDefs$type <- setReactSelRangeFilter("toxTab", reactable::colDef())
        colDefs <- setReactNumRangeFilters("toxTab", tox, colDefs)
        
        makeReactable(tox, "toxTab", columns = colDefs, pagination = TRUE, filterable = FALSE,
                      meta = list(internFilterable = "group", neverFilterable = NULL))
    },
    
    genSuspAnnTable = function(id)
    {
        tab <- as.data.table(objects$fGroups, collapseSuspects = NULL)
        mcn <- mergedConsensusNames(objects$fGroups)
        tab <- tab[, unique(c("group", "susp_name", getAllSuspCols(paste0("susp_", suspAnnCols()), names(tab), mcn))),
                   with = FALSE]
        tab[, suspID := paste0(susp_name, "-", group)][, c("susp_name", "group") := NULL]
        setnames(tab, sub("^susp_", "", names(tab)))
        
        rcols <- getAllSuspCols(c("annSimForm", "annSimComp", "annSimBoth", "maxFragMatchesRel"), names(tab), mcn)
        if (length(rcols) > 0)
            tab[, (rcols) := lapply(.SD, round, 2), .SDcols = rcols]
        
        ptab <- makePropTab(tab, if (isFGSet(objects$fGroups)) sets(objects$fGroups) else NULL, "suspID")
        makePropReactable(ptab, id, "suspID", minPropWidth = 150, minValWidth = 100)
    },
    
    genISTDGraph = function(set = NULL)
    {
        args <- list(objects$fGroups, width = "100%", height = "90%") # NOTE: bit less height to avoid vertical scrollbar
        if (!is.null(set))
            args <- c(args, set = set)
        do.call(plotGraph, args)
    },
    
    makeFGToolbar = function(tableID)
    {
        makeToolbar(tableID, columnToggles = list(
            list(value = "group", name = "Group info", checked = TRUE),
            list(value = "intensities", name = "Intensities", checked = TRUE),
            maybeInclUI(hasConcs(), list(value = "concentrations", name = "Concentrations",
                                       checked = TRUE)),
            maybeInclUI(hasFQualities(), list(value = "qualities", name = "Quality scores")),
            maybeInclUI(settings$features$chromatograms$large,
                      list(value = "chrom_large", name = "Large chromatograms"))
        ))
    },
    
    mainTabToClass = function(main) if (main) "detailsMainTable" else "detailsCandTable",
    
    makeFGTableCard = function(tab, main, ...)
    {
        makeMainTableCard(tab, ..., cl = mainTabToClass(main), hd = "Feature groups",
                          toolbar = makeFGToolbar(tab$elementId))
    },
    makeCandTableCard = function(tab, main, ..., groupBy = NULL)
    {
        tb <- makeToolbar(tab$elementId, groupBy = groupBy, toggleExpand = !is.null(groupBy))
        makeMainTableCard(tab, ..., cl = mainTabToClass(main), toolbar = tb)
    },
    
    genDetailsPlainUI = function()
    {
        list(
            makeFGTableCard(genMainTablePlain(), TRUE, "Plain")
        )
    },
    
    genDetailsSuspectsUI = function()
    {
        hasForms <- !is.null(screenInfo(objects$fGroups)[["formula"]])
        groupBy <- if (hasForms)
        {
            list(
                list(value = "", name = "None"),
                list(value = "susp_formula", name = "Formula")
            )
        }
        else
            NULL
        
        list(
            makeFGTableCard(genMainTableSusByGroup(), TRUE, "SuspectsByGroup"),
            makeCandTableCard(genMainTableSusCandSuspect(), FALSE, "SuspectsByGroup", "Candidates",
                              groupBy = groupBy),
            makeCandTableCard(genMainTableSusBySuspect(), TRUE, "SuspectsBySuspect", "Suspects",
                              groupBy = groupBy),
            makeFGTableCard(genMainTableSusCandGroup(), FALSE, "SuspectsBySuspect")
        )
    },
    
    genDetailsISTDsUI = function()
    {
        hasForms <- !is.null(internalStandards(objects$fGroups)[["formula"]])
        groupBy <- if (hasForms)
        {
            list(
                list(value = "", name = "None"),
                list(value = "susp_formula", name = "Formula")
            )
        }
        else
            NULL
        
        list(
            makeFGTableCard(genMainTableISTDsByGroup(), TRUE, "ISTDsByGroup"),
            makeCandTableCard(genMainTableISTDsCandISTD(), FALSE, "ISTDsByGroup", "Internal standards",
                              groupBy = groupBy),
            makeCandTableCard(genMainTableISTDsByISTD(), TRUE, "ISTDsByISTD", "Internal standards",
                              groupBy = groupBy),
            makeFGTableCard(genMainTableISTDsCandGroup(), FALSE, "ISTDsByISTD")
        )
    },
    
    genISTDsGraphUI = function()
    {
        layoutArgs <- list(width = 1, height = "100%", style = "padding-bottom: 10px; padding-right: 10px;")
        if (!hasSets())
        {
            do.call(bslib::layout_column_wrap, c(layoutArgs, list(
                bslib::card(
                    bslib::card_header("Internal standard assignments"),
                    bslib::card_body(genISTDGraph())
                )        
            )))
        }
        else
        {
            navs <- lapply(getFGSets(), function(s)
            {
                bslib::nav(
                    s,
                    bslib::card_body(genISTDGraph(s))
                )
            })
            do.call(bslib::layout_column_wrap, c(layoutArgs, list(
                do.call(bslib::navs_tab_card, c(list(title = "Internal standard assignments"), navs))
            )))
        }
    }
)

