#' @include main.R
#' @include report-html.R
NULL

getFeatTable <- function(fGroups, colSusp)
{
    tab <- if (isScreening(fGroups))
        as.data.table(fGroups, qualities = "score", average = TRUE, collapseSuspects = colSusp)
    else
        as.data.table(fGroups, qualities = "score", average = TRUE)

    tab <- subsetDTColumnsIfPresent(tab, c("group", "ret", "mz", replicateGroups(fGroups), "adduct", "neutralMass",
                                           paste0("susp_", c("name", "estIDLevel", "d_rt", "d_mz", "sets", "InChIKey")),
                                           featureQualityNames(scores = TRUE),
                                           grep("^ISTD_assigned", names(tab), value = TRUE)))
    
    # if (rmdVars$retMin) UNDONE
    tab[, ret := ret / 60]
    if (!is.null(tab[["susp_d_rt"]]))
        tab[, susp_d_rt := susp_d_rt / 60]
    
    for (col in names(tab)[(sapply(tab, is.numeric))])
        set(tab, j = col, value = round(tab[[col]],
                                        if (col %in% c("mz", "neutralMass", "susp_d_mz")) 5 else 2))
    tab[, (replicateGroups(fGroups)) := lapply(.SD, round, 0), .SDcols = replicateGroups(fGroups)]
    
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

getFeatColSepStyle <- function() list(borderLeft = "1px solid DarkGrey")

getFeatColGrpStartCols <- function(groupDefs) sapply(groupDefs[-1], function(col) col$columns[1])

featTabHasSusps <- function(tab) !is.null(tab[["susp_d_mz"]]) # HACK: this column should always be there if there are (non-collapsed) suspect results

getFeatColDefs <- function(tab)
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
    if (featTabHasSusps(tab))
        setCD("susp_name", "name", "suspect")
    else
        setCD("susp_name", "name", "name(s)")
    setCD("susp_d_rt", "name", "\U0394 ret")
    setCD("susp_d_mz", "name", "\U0394 mz")
    setCD("susp_estIDLevel", "name", "estIDLevel") # UNDONE: tool-tip?
    setCD("susp_estIDLevel", "align", "right")
    # InChIKeys are only there for internal usage
    setCD("susp_InChIKey", "show", FALSE)
    
    featScoreNames <- intersect(featureQualityNames(scores = TRUE), names(tab))
    if (length(featScoreNames) > 0)
    {
        for (col in featScoreNames)
        {
            setCD(col, "name", sub("Score$", "", col))
            setCD(col, "show", FALSE) # hidden by default, controlled by checkbox
        }
    }
    
    return(colDefs)
}

getFeatGroupDefs <- function(tab, groupBy, rgs)
{
    colSepStyle <- getFeatColSepStyle()
    hasSusp <- featTabHasSusps(tab)
    isGrouped <- !is.null(groupBy)
    featScoreNames <- intersect(featureQualityNames(scores = TRUE), names(tab))
    
    return(pruneList(list(
        # NOTE: the first group doesn't have a headerStyle (left border)
        # workaround for stickies: https://github.com/glin/reactable/issues/236#issuecomment-1107911895
        if (isGrouped) reactable::colGroup("", columns = groupBy, sticky = "left") else NULL,
        reactable::colGroup("feature", columns = intersect(c("group", "ret", "mz", "adduct", "neutralMass",
                                                             grep("^ISTD_assigned", names(tab), value = TRUE)),
                                                           names(tab)),
                            headerStyle = if (isGrouped) colSepStyle else NULL),
        if (hasSusp) reactable::colGroup("screening",
                                         columns = intersect(c("susp_d_rt", "susp_d_mz", "susp_estIDLevel", "sets"),
                                                             names(tab)),
                                         headerStyle = colSepStyle)
        # may still be suspects, but collapsed
        else if (!is.null(tab[["susp_name"]])) reactable::colGroup("suspect", "susp_name", headerStyle = colSepStyle) else NULL,
        reactable::colGroup("intensity", columns = rgs, headerStyle = colSepStyle),
        if (length(featScoreNames) > 0) reactable::colGroup("feature quality scores", columns = featScoreNames,
                                                            headerStyle = colSepStyle) else NULL
    )))
}

makeFeatReactable <- function(tab, id, colDefs, groupDefs, visible, EICsTopMost, plots, ..., onClick = NULL)
{
    # sync column order
    tab <- copy(tab)
    setcolorder(tab, unlist(lapply(groupDefs, "[[", "columns")))
    
    oc <- htmlwidgets::JS(sprintf("function(rowInfo, column)
{
    const tabEl = '%s';
    Reactable.setMeta(tabEl, { selectedRow: rowInfo.index });
    
    if (rowInfo.values && document.getElementById('compoundsTab'))
        Reactable.setFilter('compoundsTab', 'group', rowInfo.values.group);
    %s;
}", id, if (!is.null(onClick)) paste0("(", onClick, ")(tabEl, rowInfo, column);") else ""))
 
    # add EICs
    tab[, c("chrom_small", "chrom_large") := group]
    tabn <- names(tab)
    setcolorder(tab, c(tabn[seq_len(match("group", tabn))], "chrom_small", "chrom_large")) # move after group column
    colDefs$chrom_small <- reactable::colDef("chromatogram", cell = function(value, index)
    {
        sparkline::sparkline(EICsTopMost[[value]]$intensity, xvalues = EICsTopMost[[value]]$time, type = "line")
    })
    colDefs$chrom_large <- reactable::colDef("chromatogram", minWidth = 400, show = FALSE, cell = function(value, index)
    {
        htmltools::img(src = plots$chroms[[value]])
    })
    
    for (i in seq_along(groupDefs))
    {
        if (groupDefs[[i]]$name == "feature")
        {
            groupDefs[[i]]$columns <- c(groupDefs[[i]]$columns, "chrom_small", "chrom_large")
            break
        }
    }
    
    colSepStyle <- getFeatColSepStyle()
    grpStartCols <- getFeatColGrpStartCols(groupDefs)
    
    bgstyle <- htmlwidgets::JS(sprintf("function(rowInfo, column, state)
{
    let ret = { }
    const gby = state.groupBy;
    if (gby.length !== 0)
    {
        if (rowInfo.level === 0)
        {
            ret.background = 'black';
            ret.color = 'white';
        }
        else if (gby.length > 1 && rowInfo.level === 1 && column.id !== gby[0])
            ret.background = 'LightGrey';
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
        qualities = as.vector(groupCols[["feature quality scores"]]),
        chrom_large = "chrom_large"
    )

    headThemeStyle <- list(padding = "2px 4px")
    rt <- reactable::reactable(tab, elementId = id, pagination = FALSE, wrap = FALSE, resizable = TRUE,
                               highlight = TRUE, bordered = TRUE, onClick = oc, defaultExpanded = TRUE,
                               columns = colDefs, defaultColDef = reactable::colDef(style = bgstyle),
                               columnGroups = groupDefs,
                               theme = reactable::reactableTheme(headerStyle = headThemeStyle,
                                                                 groupHeaderStyle = headThemeStyle,
                                                                 cellPadding = "2px 4px"),
                               meta = list(selectedRow = NULL, plots = plots, colToggles = colToggles),
                               rowStyle = htmlwidgets::JS("function(rowInfo, state)
{
    const sel = state.meta.selectedRow;
    let ret = { cursor: 'pointer' };
    if (sel != null && rowInfo.index === sel)
        ret.background = '#eee';
    return ret;
}"), ...)
    
    if (!visible)
        rt <- htmlwidgets::onRender(rt, htmlwidgets::JS("function(el, x) { el.style.display = 'none'; }"))
    
    return(rt)
}

makeAnnReactable <- function(tab, id, ...)
{
    return(reactable::reactable(tab, elementId = id, resizable = TRUE, bordered = TRUE, wrap = FALSE,
                                pagination = FALSE, compact = TRUE,
                                language = reactable::reactableLang(noData = "No annotations available"), ...))
}

reportHTMLUtils$methods(
    genFeatTablePlain = function()
    {
        tab <- getFeatTable(objects$fGroups, ",")
        groupDefs <- getFeatGroupDefs(tab, NULL, replicateGroups(objects$fGroups))
        colDefs <- getFeatColDefs(tab)
        makeFeatReactable(tab, "detailsTabPlain", colDefs = colDefs, groupDefs = groupDefs, visible = TRUE,
                          EICsTopMost, plots = plots)
    },
    genFeatTableSuspects = function()
    {
        tab <- getFeatTable(objects$fGroups, NULL)
        groupDefs <- getFeatGroupDefs(tab, "susp_name", replicateGroups(objects$fGroups))
        colDefs <- getFeatColDefs(tab)
        
        onClick <- "function(tabEl, rowInfo)
{
    const structEl = document.getElementById('struct_view-suspect');
    const rd = (rowInfo.level === 0) ? rowInfo.subRows[0] : rowInfo.values;
    structEl.src = Reactable.getState(tabEl).meta.plots.structs[rd.susp_InChIKey];
}"
        makeFeatReactable(tab, "detailsTabSuspects", colDefs = colDefs, groupDefs = groupDefs, visible = FALSE,
                          EICsTopMost, plots = plots, groupBy = "susp_name", onClick = onClick)
    },
    genFeatTableComponents = function()
    {
        tab <- getFeatTable(objects$fGroups, ",")
        groupDefs <- getFeatGroupDefs(tab, "component", replicateGroups(objects$fGroups))
        colDefs <- getFeatColDefs(tab)
        
        ctab <- as.data.table(objects$components)
        setnames(ctab, "name", "component")
        ctab <- removeDTColumnsIfPresent(ctab, c(
            # already present from features table
            "ret", "mz",
            
            # general component properties
            "mz_increment", "rt_increment", "ret_min", "ret_max", "ret_range",
            replicateGroups(objects$fGroups), # nontarget: presence in replicate groups
            "cmp_ret", "cmp_mz", "cmp_retsd", "neutral_mass", "cmp_ppm", "analysis", "size",
            
            # for graphs
            "links"
        ))
        if (all(ctab$intensity == 1))
            ctab[, intensity := NULL]
        tab <- merge(tab, ctab, by = "group", sort = FALSE)
        
        cmpGrpCols <- setdiff(names(ctab), c("component", "group"))
        if (length(cmpGrpCols) > 0)
            groupDefs <- c(groupDefs, list(reactable::colGroup("component", cmpGrpCols, headerStyle = getFeatColSepStyle())))
        
        onClick <- "function(tabEl, rowInfo)
{
    let chromEl = document.getElementById('chrom_view-component');
    let specEl = document.getElementById('spectrum_view-component');
    let profileRelEl = document.getElementById('profileRel_view-component');
    let profileAbsEl = document.getElementById('profileAbs_view-component');
    const rd = (rowInfo.level === 0) ? rowInfo.subRows[0] : rowInfo.values;
    const pl = Reactable.getState(tabEl).meta.plots.components[rd.component];
    chromEl.src = pl.chrom;
    specEl.src = pl.spec;
    if (profileRelEl != undefined)
    {
        profileRelEl.src = pl.profileRel;
        profileAbsEl.src = pl.profileAbs;
    }
}"
        
        makeFeatReactable(tab, "detailsTabComponents", colDefs = colDefs, groupDefs = groupDefs, visible = FALSE,
                          EICsTopMost, plots = plots, groupBy = "component", onClick = onClick)
    },
    genFeatTableTPs = function()
    {
        fromTPs <- objects$components@fromTPs
        
        tabTPsFeat <- getFeatTable(objects$fGroups, if (fromTPs) NULL else ",")
        
        tabCompon <- as.data.table(objects$components)
        tabCompon <- subsetDTColumnsIfPresent(tabCompon, c("name", "parent_name", "parent_group", "group", "TP_retDir",
                                                           "TP_name", "retDir", "retDiff", "mzDiff", "formulaDiff",
                                                           "specSimilarity", "mergedBy"))

        if (fromTPs)
        {
            tabTPs <- merge(tabCompon, tabTPsFeat, by.x = c("group", "TP_name"), by.y = c("group", "susp_name"))
            setnames(tabTPs, c("name", "TP_name", "parent_name"),
                     c("component", "susp_name", "parent_susp_name"), skip_absent = TRUE)
        }
        else
        {
            tabTPs <- merge(tabCompon, tabTPsFeat, by = "group")
            setnames(tabTPs, "name", "component")
            tabTPs <- removeDTColumnsIfPresent(tabTPs, "susp_name")
        }

        parAggr <- function(parentCol) htmlwidgets::JS(sprintf("function(values, rows, groupRows)
{
    const v = rows[0]['%s'];
    return (v == null) ? '' : '<i>' + v + '</i>';
}", parentCol))
        
        parAggred <- function(level) htmlwidgets::JS(sprintf("function(cellInfo, state)
{
    return (cellInfo.level === %d) ? cellInfo.value : '';
}", level))

        # NOTE: below values may be in components but then from suspect list
        # UNDONE: convert RT to minutes if needed
        tabTPs[, c("parent_ret", "parent_mz") := groupInfo(objects$fGroups)[parent_group, ]]
        
        for (col in intersect(c("parent_ret", "retDiff", "specSimilarity"), names(tabTPs)))
            tabTPs[, (col) := round(get(col), 2)]
        for (col in c("parent_mz", "mzDiff"))
            tabTPs[, (col) := round(get(col), 5)]

        # add parent intensities & screening info
        rgs <- replicateGroups(objects$fGroups)
        tabTPsPar <- unique(subsetDTColumnsIfPresent(tabTPsFeat,
                                                     c("group", rgs,
                                                       paste0("susp_", c("estIDLevel", "d_rt", "d_mz", "sets",
                                                                         "InChIKey")))),
                            by = "group")
        setnames(tabTPsPar, paste0("parent_", names(tabTPsPar)))
        tabTPs <- merge(tabTPs, tabTPsPar, by = "parent_group", sort = FALSE, all.x = TRUE)
        
        groupBy <- if (fromTPs) c("component", "susp_name") else "component"
        groupDefs <- getFeatGroupDefs(tabTPs, groupBy, rgs)
        # squeeze in TP column
        groupDefs <- c(groupDefs[1:2],
                       list(reactable::colGroup("TP", columns = intersect(c("TP_name", "retDiff", "mzDiff",
                                                                            "formulaDiff", "retDir",
                                                                            "specSimilarity", "mergedBy"),
                                                                          names(tabTPs)),
                                                headerStyle = getFeatColSepStyle())),
                       groupDefs[seq(3, length(groupDefs))])
        
        colDefs <- getFeatColDefs(tabTPs)
        
        # set parent 'aggregates': actual value of parent feature group
        for (col in grep("^parent_", names(tabTPs), value = TRUE))
        {
            colTP <- sub("^parent_", "", col)
            colDefs[[colTP]]$aggregate <- parAggr(col)
            colDefs[[colTP]]$aggregated <- parAggred(0)
            colDefs[[colTP]]$html <- TRUE
            colDefs[[col]] <- reactable::colDef(show = FALSE)
        }
        # similar for TP retDir
        colDefs$retDir <- reactable::colDef(aggregate = parAggr("TP_retDir"), aggregated = parAggred(1), html = TRUE)
        colDefs$TP_retDir <- reactable::colDef(show = FALSE)
        
        # these are grouped
        if (!is.null(tabTPs[["TP_name"]]))
            colDefs$TP_name <- reactable::colDef("name")
        colDefs$retDiff <- reactable::colDef(name = "\U0394 ret")
        colDefs$mzDiff <- reactable::colDef(name = "\U0394 mz")
        colDefs$formulaDiff <- reactable::colDef(name = "\U0394 formula")
        
        # InChIKeys are only there for internal usage
        if (!is.null(tabTPs[["parent_susp_InChIKey"]]))
            colDefs$parent_susp_InChIKey <- reactable::colDef(show = FALSE)
        
        onClick <- "function(tabEl, rowInfo)
{
    const chromEl = document.getElementById('chrom_view-tp');
    const structEl = document.getElementById('struct_view-tp');
    let rd;
    
    // check level: handle clicking on expandable rows
    if (rowInfo.level === 0)
        rd = rowInfo.subRows[0]._subRows[0];
    else if (rowInfo.level === 1)
        rd = rowInfo.subRows[0];
    else
        rd = rowInfo.values;
        
    chromEl.src = Reactable.getState(tabEl).meta.plots.chroms[rd.parent_group];
    structEl.src = Reactable.getState(tabEl).meta.plots.structs[rd.parent_susp_InChIKey];
    showTPGraph(rd.component);
}"
    
        makeFeatReactable(tabTPs, "detailsTabTPs", FALSE, EICsTopMost, plots, groupBy = groupBy, colDefs = colDefs,
                          groupDefs = groupDefs, onClick = onClick)
    },
    
    genCompoundTable = function()
    {
        compounds <- objects$compounds[names(objects$fGroups)]
        
        mcn <- mergedConsensusNames(compounds)
        
        if (length(compounds) == 0)
        {
            # dummy table to show empty results
            return(makeAnnReactable(data.table(compound = character()), "compoundsTab"))
        }

        tab <- as.data.table(compounds)
        
        # NOTE: for consensus results, duplicate algo columns (eg identifier) are only shown in details
        tab <- subsetDTColumnsIfPresent(tab, c("group", "compoundName", "compoundName2", "identifier", "database",
                                               "neutral_formula", "neutralMass", "explainedPeaks", "score", "InChIKey"))
        
        cmpIndices <- tab[, seq_len(.N), by = "group"][[2]]
        cmpNames2 <- tab[["compoundName2"]]
        if (!is.null(cmpNames2))
            tab[, compoundName2 := NULL]
        
        if (!is.null(tab[["identifier"]]))
            tab[, identifier := sapply(identifier, makeDBIdentLink, db = database[1])][, database := NULL]
        
        tab[, neutral_formula := subscriptFormulaHTML(neutral_formula)]
        tab[, neutralMass := round(neutralMass, 5)]
        if (!is.null(tab[["score"]]))
            tab[, score := round(score, 2)]
        
        tab[, structure := plots$structs[InChIKey]][, InChIKey := NULL]
        tab[, spectrum := plots$compounds[[group]]$spectra, by = "group"]
        tab[, scorings := plots$compounds[[group]]$scores, by = "group"]
        
        getImgCell <- function(value) htmltools::img(src = value, style = list("max-height" = "300px"))
        
        makeDetailsTable <- function(title, ...)
        {
            # Nested table: based on from reactable cookbook
            return(htmltools::div(style = list(margin = "10px 20px"),
                                  htmltools::div(style = list("text-align" = "center", "font-weight" = "bold"), title),
                                  reactable::reactable(pagination = FALSE, compact = TRUE, bordered = TRUE, wrap = FALSE,
                                                       fullWidth = FALSE, resizable = TRUE, striped = TRUE,
                                                       height = 200, ...)))
        }
        
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
            cit <- getCompInfoTable(compounds[[tab$group[index]]], cmpIndices[index], mcn, TRUE)
            cit <- cit[!property %in% names(tab)]
            if (!is.null(tab[["compoundName"]]))
                cit <- cit[property != "compoundName2"] # already in main table compoundName column
            return(makeDetailsTable("Compound properties", cit, columns = list(
                property = reactable::colDef(minWidth = 150),
                value = reactable::colDef(html = TRUE, minWidth = 250)
            )))
        }
        
        getAnnPLDetails <- function(index)
        {
            apl <- annotatedPeakList(compounds, index = cmpIndices[index], groupName = tab$group[index],
                                     MSPeakLists = objects$MSPeakLists, formulas = objects$formulas,
                                     onlyAnnotated = TRUE)
            
            if (is.null(apl) || nrow(apl) == 0)
                return(htmltools::div(align = "center", "No annotation available."))
            
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
            
            return(makeDetailsTable("Peak list annotations", apl, columns = colDefs,
                                    rowClass = function(index) if (isPrec[index]) "font-weight-bold" else ""))
        }
        
        getScoreDetails <- function(index)
        {
            cRow <- compounds[[tab$group[index]]][cmpIndices[index]]
            sc <- getAllMergedConsCols(annScoreNames(compounds, FALSE), names(cRow), mcn)
            
            scores <- cRow[, sc, with = FALSE]
            scores <- setnames(transpose(scores, keep.names = "score"), 2, "value")
            return(makeDetailsTable("Scorings", scores, columns = list(
                score = reactable::colDef(minWidth = 150),
                value = reactable::colDef(format = reactable::colFormat(digits = 2),
                                          minWidth = 100)
                )))
        }
        
        setcolorder(tab, intersect(c("compoundName", "structure"), names(tab)))
        
        colDefs <- pruneList(list(
            group = reactable::colDef(show = FALSE),
            compoundName = if (!is.null(tab[["compoundName"]])) reactable::colDef("compound", cell = getCompCell) else NULL,
            identifier = if (!is.null(tab[["identifier"]])) reactable::colDef(html = TRUE) else NULL,
            neutral_formula = reactable::colDef("formula", html = TRUE),
            neutralMass = reactable::colDef("neutral mass"),
            structure = reactable::colDef(cell = getImgCell, minWidth = 125),
            spectrum = reactable::colDef(cell = getImgCell, minWidth = 200),
            scorings = reactable::colDef(cell = getImgCell, minWidth = 200)
        ))
        
        return(makeAnnReactable(tab, "compoundsTab", columns = colDefs, details = function(index)
        {
            htmltools::div(style = list(margin = "12px 45px", display = "flex", "flex-wrap" = "no-wrap",
                                        background = "#FCFCFC", border = "dashed 1px",
                                        "justify-content" = "space-between", "overflow-x" = "auto"),
                           getCompDetails(index), getAnnPLDetails(index), getScoreDetails(index))
        }))
    }
)

