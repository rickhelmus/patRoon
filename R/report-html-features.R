#' @include main.R
#' @include report-html.R
NULL

roundFGTab <- function(ftab, fGroups)
{
    ftab <- copy(ftab)
    for (col in names(ftab)[(sapply(ftab, is.numeric))])
        set(ftab, j = col, value = round(ftab[[col]],
                                         if (col %in% c("mz", "neutralMass", "susp_d_mz")) 5 else 2))
    ftab[, (replicateGroups(fGroups)) := lapply(.SD, round, 0), .SDcols = replicateGroups(fGroups)]
    return(ftab)    
}

getFGTable <- function(fGroups, colSusp, retMin, concAggrParams, toxAggrParams)
{
    adtArgs <- list(fGroups, qualities = "score", average = TRUE, concAggrParams = concAggrParams,
                    toxAggrParams = toxAggrParams)
    if (isScreening(fGroups))
        adtArgs <- c(adtArgs, list(collapseSuspects = colSusp))
    tab <- do.call(as.data.table, adtArgs)
    tab <- subsetDTColumnsIfPresent(tab, c("group", "ret", "mz", replicateGroups(fGroups), "adduct", "neutralMass",
                                           paste0(replicateGroups(fGroups), "_conc"), "conc_types", "LC50", "LC50_types",
                                           paste0("susp_", c("name", "estIDLevel", "d_rt", "d_mz", "sets", "InChIKey")),
                                           featureQualityNames(scores = TRUE),
                                           grep("^ISTD_assigned", names(tab), value = TRUE)))
    
    if (retMin)
    {
        tab[, ret := ret / 60]
        if (!is.null(tab[["susp_d_rt"]]))
            tab[, susp_d_rt := susp_d_rt / 60]
    }
    
    tab <- roundFGTab(tab, fGroups)
    
    if (nrow(internalStandards(fGroups)) > 0)
    {
        wrapISTDs <- function(s) wrapStr(gsub(",", ", ", s, fixed = TRUE), 50)
        if (isFGSet(fGroups))
        {
            for (s in sets(fGroups))
            {
                cn <- paste0("ISTD_assigned-", s)
                tab[, (cn) := sapply(get(cn), wrapISTDs)]
            }
        }
        else if (!is.null(tab[["ISTD_assigned"]]))
            tab[, ISTD_assigned := sapply(ISTD_assigned, wrapISTDs)]
    }
    
    return(tab)
}

getFGColSepStyle <- function() list(borderLeft = "1px solid DarkGrey")

getFGColGrpStartCols <- function(groupDefs) sapply(groupDefs[-1], function(col) col$columns[1])

featGroupTabHasSusps <- function(tab) !is.null(tab[["susp_d_mz"]]) # HACK: this column should always be there if there are (non-collapsed) suspect results

getFeatGroupColDefs <- function(tab)
{
    colDefs <- list()
    
    setCD <- function(col, field, value)
    {
        if (col %chin% names(tab))
        {
            if (is.null(colDefs[[col]]))
                colDefs[[col]] <<- do.call(reactable::colDef, setNames(list(value), field))
            else
                colDefs[[col]][[field]] <<- value
        }
    }
    
    setCD("mz", "name", "m/z")
    if (featGroupTabHasSusps(tab))
        setCD("susp_name", "name", "suspect")
    else
        setCD("susp_name", "name", "name(s)")
    setCD("susp_d_rt", "name", "\U0394 ret")
    setCD("susp_d_mz", "name", "\U0394 mz")
    setCD("susp_estIDLevel", "name", "estIDLevel") # UNDONE: tool-tip?
    setCD("susp_estIDLevel", "align", "right")
    IDLCols <- getBrewerPal(5, "YlOrBr")
    setCD("susp_estIDLevel", "cell", function(value)
    {
        nidl <- numericIDLevel(value)
        htmltools::span(style = list("border-radius" = "30px", display = "inline-block", width = "2em",
                                     "background-color" = IDLCols[nidl], color = if (nidl <= 2) "black" else "white",
                                     "text-align" = "center", "font-weight" = "bold"),
                        value)
    })
    setCD("susp_sets", "name", "sets")
    # InChIKeys are only there for internal usage
    setCD("susp_InChIKey", "show", FALSE)
    
    featScoreNames <- intersect(featureQualityNames(scores = TRUE), names(tab))
    for (col in featScoreNames)
    {
        setCD(col, "name", sub("Score$", "", col))
        setCD(col, "show", FALSE) # hidden by default, controlled by checkbox
    }
    
    return(colDefs)
}

getFGGroupDefs <- function(tab, groupBy, rgs)
{
    colSepStyle <- getFGColSepStyle()
    hasSusp <- featGroupTabHasSusps(tab)
    isGrouped <- !is.null(groupBy)
    featScoreNames <- intersect(featureQualityNames(scores = TRUE), names(tab))
    concCols <- intersect(c(paste0(rgs, "_conc"), "conc_types"), names(tab))
    
    return(pruneList(list(
        # NOTE: the first group doesn't have a headerStyle (left border)
        # workaround for stickies: https://github.com/glin/reactable/issues/236#issuecomment-1107911895
        if (isGrouped) reactable::colGroup("", columns = groupBy, sticky = "left") else NULL,
        reactable::colGroup("feature", columns = intersect(c("group", "ret", "mz", "adduct", "neutralMass",
                                                             grep("^ISTD_assigned", names(tab), value = TRUE)),
                                                           names(tab)),
                            headerStyle = if (isGrouped) colSepStyle else NULL),
        if (hasSusp) reactable::colGroup("screening",
                                         columns = intersect(c("susp_d_rt", "susp_d_mz", "susp_estIDLevel", "susp_sets"),
                                                             names(tab)),
                                         headerStyle = colSepStyle)
        # may still be suspects, but collapsed
        else if (!is.null(tab[["susp_name"]])) reactable::colGroup("suspect", "susp_name", headerStyle = colSepStyle) else NULL,
        reactable::colGroup("intensity", columns = rgs, headerStyle = colSepStyle),
        if (length(concCols) > 0) reactable::colGroup("concentrations", columns = concCols, headerStyle = colSepStyle) else NULL,
        if (!is.null(tab[["LC50"]])) reactable::colGroup("toxicity", columns = c("LC50", "LC50_types"),
                                                         headerStyle = colSepStyle) else NULL,
        if (length(featScoreNames) > 0) reactable::colGroup("feature quality scores", columns = featScoreNames,
                                                            headerStyle = colSepStyle) else NULL
    )))
}

makeFGReactable <- function(tab, id, colDefs, groupDefs, visible, plots, settings, objects, ...)
{
    addToGroupDefs <- function(gd, grp, col)
    {
        for (i in seq_along(gd))
        {
            if (gd[[i]]$name == grp)
            {
                gd[[i]]$columns <- c(gd[[i]]$columns, col)
                break
            }
        }
        return(gd)
    }

    # sync column order
    tab <- copy(tab)
    setcolorder(tab, unlist(lapply(groupDefs, "[[", "columns")))
    
    if (settings$features$chromatograms$small || settings$features$chromatograms$large)
    {
        cols <- character()
        if (settings$features$chromatograms$small)
            cols <- "chrom_small"
        if (settings$features$chromatograms$large)
            cols <- c(cols, "chrom_large")
        
        # add EICs
        tab[, (cols) := group]
        tabn <- names(tab)
        setcolorder(tab, c(tabn[seq_len(match("group", tabn))], cols)) # move after group column
        
        if (settings$features$chromatograms$small)
        {
            cell <- function(value, index)
            {
                tag <- htmltools::img(src = plots$chromsSmall[[value]], style = list("max-height" = "20px"))
                if (settings$features$chromatograms$large)
                    tag <- htmltools::tagAppendAttributes(tag, "data-srcZoom" = plots$chromsLarge[[value]])
                else
                    tag <- htmltools::tagAppendAttributes(tag, class = "noZoomImg")
                return(tag)
            }
            colDefs$chrom_small <- reactable::colDef("chromatogram", filterable = FALSE, searchable = FALSE,
                                                     cell = cell)
        }
        if (settings$features$chromatograms$large)
        {
            colDefs$chrom_large <- reactable::colDef("chromatogram", minWidth = 400, show = FALSE, filterable = FALSE,
                                                     searchable = FALSE,
                                                     cell = function(value, index) htmltools::img(src = plots$chromsLarge[[value]]))
        }
        
        groupDefs <- addToGroupDefs(groupDefs, "feature", cols)
    }
    
    if (any(c("MSPeakLists", "formulas", "compounds") %in% names(objects)))
    {
        tab[, annotations := group]
        tab[, hasMSMS := sapply(group, function(g) !is.null(objects$MSPeakLists) &&
                                    !is.null(objects$MSPeakLists[[g]]) && !is.null(objects$MSPeakLists[[g]][["MSMS"]]))]
        tab[, hasFormulas := sapply(group, function(g) !is.null(objects$formulas) && !is.null(objects$formulas[[g]]))]
        tab[, hasCompounds := sapply(group, function(g) !is.null(objects$compounds) && !is.null(objects$compounds[[g]]))]
        
        annTag <- function(bg, fg, ...) htmltools::span(style = list("border-radius" = "30px", padding = "4px 5px",
                                                                     "margin-right" = "2px", "font-size" = "smaller",
                                                                     "background-color" = bg, color = fg), ...)
        colDefs$annotations <- reactable::colDef(width = 150, cell = function(value, index)
        {
            tags <- list()
            if (tab$hasMSMS[index])
                tags <- c(tags, list(annTag("#033c73", "white", "MS", htmltools::tags$sup("2"))))
            if (tab$hasFormulas[index])
                tags <- c(tags, list(annTag("#2d9ddd", "white", "Form")))
            if (tab$hasCompounds[index])
                tags <- c(tags, list(annTag("grey", "white", "Comp")))
            return(do.call(htmltools::tagList, tags))
        }, filterInput = function(values, name) reactSelFilterButton(id, name, "#filterSelModal",
                                                                     "filtFeatAnnSelModalInit", "Select"),
        filterMethod = htmlwidgets::JS("function(rows, columnId, filterValue)
{
    const MSMS = filterValue.has('MS/MS'), forms = filterValue.has('Formulas'), comps = filterValue.has('Compounds'),
                 none = filterValue.has('None');
    return rows.filter(function(row)
    {
        const rMSMS = row.values.hasMSMS, rforms = row.values.hasFormulas, rcomps = row.values.hasCompounds;
        return((none && !MSMS && !rforms && !rcomps) || MSMS && rMSMS || forms && rforms || comps && rcomps);
    })                                       
}"))
        groupDefs <- addToGroupDefs(groupDefs, "feature", "annotations")
        
        colDefs[c("hasMSMS", "hasFormulas", "hasCompounds")] <- list(reactable::colDef(show = FALSE))
    }
    
    
    colSepStyle <- getFGColSepStyle()
    grpStartCols <- getFGColGrpStartCols(groupDefs)
    
    bgstyle <- htmlwidgets::JS(sprintf("function(rowInfo, column, state)
{
    let ret = { }
    const gby = state.groupBy;
    if (gby.length !== 0)
    {
        if (rowInfo.level === 0)
        {
            ret.background = 'lightcyan';
            //ret.color = 'white';
        }
        else if (gby.length > 1 && rowInfo.level === 1 && column.id !== gby[0])
            ret.background = '#f7ffff';
    }
    
    if ([ %s ].includes(column.id))
        ret.borderLeft = '%s';

    return ret;
}", paste0("'", grpStartCols, "'", collapse = ","), colSepStyle))
    
    for (col in grpStartCols)
        colDefs[[col]]$headerStyle <- colSepStyle
    
    colDefs <- lapply(colDefs, function(cd)
    {
        cd$style <- bgstyle
        return(cd)
    })
    
    groupCols <- lapply(groupDefs, "[[", "columns")
    names(groupCols) <- sapply(groupDefs, "[[", "name")
    colToggles <- list(
        group = setdiff(as.vector(groupCols$feature), c("group", "chrom_small", "chrom_large")), # don't toggle group and chroms
        intensities = as.vector(groupCols$intensity),
        concentrations = as.vector(groupCols[["concentrations"]]),
        qualities = as.vector(groupCols[["feature quality scores"]]),
        chrom_large = "chrom_large"
    )
    
    colDefs <- setReactNumRangeFilters(id, tab, colDefs)
    
    onClick = htmlwidgets::JS("function(rowInfo, column)
{
    updateFeatTabRowSel(rowInfo.values, rowInfo.index);
}")
    
    CSVCols <- setdiff(names(tab), c("chrom_small", "chrom_large", "annotations", "hasMSMS", "hasFormulas",
                                     "hasCompounds"))
    
    headThemeStyle <- list(padding = "2px 4px")
    rt <- makeReactable(tab, id, highlight = TRUE, onClick = onClick, defaultExpanded = TRUE, columns = colDefs,
                        defaultColDef = reactable::colDef(style = bgstyle), columnGroups = groupDefs,
                        filterable = FALSE, pagination = TRUE,
                        theme = reactable::reactableTheme(headerStyle = headThemeStyle,
                                                          groupHeaderStyle = headThemeStyle,
                                                          cellPadding = "2px 4px"),
                        meta = list(selectedRow = 0, colToggles = colToggles, CSVCols = CSVCols),
                        rowStyle = htmlwidgets::JS("function(rowInfo, state)
{
    const sel = state.meta.selectedRow;
    let ret = { cursor: 'pointer' };
    if (sel != null && rowInfo != undefined && rowInfo.index === sel)
    {
        ret.background = '#eee';
        ret.fontWeight = 'bold';
    }
    return ret;
}"), ...)
    
    if (!visible)
        rt <- htmlwidgets::onRender(rt, htmlwidgets::JS("function(el, x) { el.style.display = 'none'; }"))
    
    return(rt)
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
                                EICParams = EICParams, colourBy = "rGroups", title = "", bty = "l"),
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
        mainArgs <- c(mainArgs, list(sets = TRUE))
    else
        mainArgs <- c(mainArgs, list(col = "black"))
    
    doApply("sapply", parallel, names(fGroups), function(grp)
    {
        doProgress()
        makeHTMLReportPlot("int_plot", outPath, "plotInt", c(list(fGroups[, grp]), mainArgs),
                           parParams = list(mar = c(4.1, 4.1, 1, 0.1)), width = 8, height = 4, bg = "transparent",
                           pointsize = 16)
    }, simplify = FALSE)
}


reportHTMLUtils$methods(
    genFGTablePlain = function()
    {
        mdprintf("Feature groups... ")
        tab <- getFGTable(objects$fGroups, ",", settings$features$retMin, settings$features$aggregateConcs,
                          settings$features$aggregateTox)
        groupDefs <- getFGGroupDefs(tab, NULL, replicateGroups(objects$fGroups))
        colDefs <- getFeatGroupColDefs(tab)
        makeFGReactable(tab, "detailsTabPlain", colDefs = colDefs, groupDefs = groupDefs, visible = TRUE, plots = plots,
                        settings = settings, objects = objects)
    },
    genFGTableSuspects = function()
    {
        tab <- getFGTable(objects$fGroups, NULL, settings$features$retMin, settings$features$aggregateConcs,
                          settings$features$aggregateTox)[!is.na(susp_name)]
        groupDefs <- getFGGroupDefs(tab, "susp_name", replicateGroups(objects$fGroups))
        colDefs <- getFeatGroupColDefs(tab)
        makeFGReactable(tab, "detailsTabSuspects", colDefs = colDefs, groupDefs = groupDefs, visible = FALSE,
                        plots = plots, settings = settings, objects = objects, groupBy = "susp_name")
    },
    genFGTableISTDs = function()
    {
        istds <- data.table::copy(internalStandards(objects$fGroups))
        istds <- subsetDTColumnsIfPresent(istds, c("name", "group", "InChIKey", "d_rt", "d_mz", "sets"))
        if (settings$features$retMin)
            istds[, d_rt := d_rt / 60]
        
        # HACK: the ISTD table is essentially the same as what screenInfo() returns for suspects. Rename the columns
        # here, so that groupDefs etc are set like suspects.
        rncols <- setdiff(names(istds), "group")
        setnames(istds, rncols, paste0("susp_", rncols))
        
        ftab <- getFGTable(objects$fGroups, ",", settings$features$retMin, settings$features$aggregateConcs,
                           settings$features$aggregateTox)
        ftab <- removeDTColumnsIfPresent(ftab, c("susp_name", "susp_sets", grep("^ISTD_assigned",
                                                                                names(ftab), value = TRUE)))
        tab <- merge(ftab, istds, by = "group")
        tab <- roundFGTab(tab, objects$fGroups)
        
        groupDefs <- getFGGroupDefs(tab, "susp_name", replicateGroups(objects$fGroups))
        colDefs <- getFeatGroupColDefs(tab)
        colDefs$susp_name$name <- "Internal standard" # HACK
        makeFGReactable(tab, "detailsTabISTDs", colDefs = colDefs, groupDefs = groupDefs, visible = FALSE,
                        plots = plots, settings = settings, objects = objects, groupBy = "susp_name")
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
            set(tab, j = col, value = round(tab[[col]], if (col %in% c("mz", "mzmin", "mzmax")) 5 else 2))
        
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
        
        tab[, rGroup := anaInfo[match(analysis, anaInfo$analysis), "group"]]
        
        colDefs <- list(
            group = reactable::colDef(show = FALSE, filterMethod = reactExactFilter()),
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
        {
            colDefs$set <- reactable::colDef(filterInput = function(values, name) reactSelectFilter("featuresTab",
                                                                                                    values, name))
        }
        
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
                      meta = list(featQualCols = fqn, CSVCols = CSVCols), pagination = TRUE)
    },
    
    genConcsTable = function()
    {
        concs <- data.table::copy(concentrations(objects$fGroups))
        
        imgTag <- function(IK, cand) sprintf("<img src='%s' alt='%s' style='max-height: 300px;')></img>", plots$structs[IK], cand)
        
        concs[type == "SIRIUS_FP", candidate := subscriptFormulaHTML(candidate)]
        concs[type == "suspect", candidate := imgTag(screenInfo(objects$fGroups)[match(candidate, SMILES)]$InChIKey,
                                                     candidate)]
        if (!is.null(objects[["compounds"]]))
        {
            cTab <- as.data.table(objects$compounds)
            concs[type == "compound", candidate := imgTag(cTab[match(candidate, SMILES)]$InChIKey, candidate)]
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
        
        makeReactable(concs, "concsTab", columns = colDefs, pagination = TRUE, filterable = FALSE)
    },

    genToxTable = function()
    {
        tox <- data.table::copy(toxicities(objects$fGroups))
        
        imgTag <- function(IK, cand) sprintf("<img src='%s' alt='%s' style='max-height: 300px;')></img>", plots$structs[IK], cand)
        
        tox[type == "SIRIUS_FP", candidate := subscriptFormulaHTML(candidate)]
        tox[type == "suspect", candidate := imgTag(screenInfo(objects$fGroups)[match(candidate, SMILES)]$InChIKey,
                                                     candidate)]
        if (!is.null(objects[["compounds"]]))
        {
            cTab <- as.data.table(objects$compounds)
            tox[type == "compound", candidate := imgTag(cTab[match(candidate, SMILES)]$InChIKey, candidate)]
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
        
        makeReactable(tox, "toxTab", columns = colDefs, pagination = TRUE, filterable = FALSE)
    },
    
    genSuspAnnTable = function()
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
        makePropReactable(ptab, "suspAnnTab", "suspID", minPropWidth = 150, minValWidth = 100)
    },
    
    genISTDGraph = function(set = NULL)
    {
        args <- list(objects$fGroups, width = "100%", height = "90%") # NOTE: bit less height to avoid vertical scrollbar
        if (!is.null(set))
            args <- c(args, set = set)
        do.call(plotGraph, args)
    }
)

