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
                                           paste0("susp_", c("name", "estIDLevel", "neutralMass", "formula", "d_rt",
                                                             "d_mz", "sets", "InChIKey"))))
    
    # if (rmdVars$retMin) UNDONE
    tab[, ret := ret / 60]
    if (!is.null(tab[["susp_d_rt"]]))
        tab[, susp_d_rt := susp_d_rt / 60]
    
    for (col in names(tab)[(sapply(tab, is.numeric))])
        set(tab, j = col, value = round(tab[[col]],
                                        if (col %in% c("mz", "neutralMass", "susp_neutralMass", "susp_d_mz")) 5 else 2))
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

getFeatColDefs <- function(tab, groupDefs, EICsTopMost)
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

    setCD("group", "cell", function(value, index)
    {
        htmltools::div(value,
                       htmltools::br(),
                       sparkline::sparkline(EICsTopMost[[value]]$intensity, xvalues = EICsTopMost[[value]]$time,
                                            type = "line"))
    })
    
    setCD("mz", "name", "m/z")
    if (featTabHasSusps(tab))
        setCD("susp_name", "name", "suspect")
    else
        setCD("susp_name", "name", "name(s)")
    setCD("susp_neutralMass", "name", "neutralMass")
    setCD("susp_formula", "name", "formula")
    setCD("susp_d_rt", "name", "\U0394 ret")
    setCD("susp_d_mz", "name", "\U0394 mz")
    setCD("susp_estIDLevel", "name", "estIDLevel") # UNDONE: tool-tip?
    setCD("susp_estIDLevel", "align", "right")
    # InChIKeys are only there for internal usage
    setCD("susp_InChIKey", "show", FALSE)
    
    colSepStyle <- getFeatColSepStyle()
    for (col in getFeatColGrpStartCols(groupDefs))
        colDefs[[col]]$headerStyle <- colSepStyle
    
    return(colDefs)
}

getFeatGroupDefs <- function(tab, groupBy, rgs)
{
    colSepStyle <- getFeatColSepStyle()
    hasSusp <- featTabHasSusps(tab)
    isGrouped <- !is.null(groupBy)
    
    return(pruneList(list(
        # NOTE: the first group doesn't have a headerStyle (left border)
        # workaround for stickies: https://github.com/glin/reactable/issues/236#issuecomment-1107911895
        if (isGrouped) reactable::colGroup("", columns = groupBy, sticky = "left") else NULL,
        reactable::colGroup("feature", columns = c("group", "ret", "mz"), headerStyle = if (isGrouped) colSepStyle else NULL),
        if (hasSusp) reactable::colGroup("screening",
                                         columns = intersect(c("susp_neutralMass", "susp_formula", "susp_d_rt",
                                                               "susp_d_mz", "susp_estIDLevel"),
                                                             names(tab)),
                                         headerStyle = colSepStyle)
        # may still be suspects, but collapsed
        else if (!is.null(tab[["susp_name"]])) reactable::colGroup("suspect", "susp_name", headerStyle = colSepStyle) else NULL,
        reactable::colGroup("intensity", columns = rgs, headerStyle = colSepStyle)
    )))
}

makeFeatReactable <- function(tab, id, colDefs, groupDefs, visible, plots, ..., onClick = NULL)
{
    oc <- htmlwidgets::JS(sprintf("function(rowInfo, column)
{
    const tabEl = '%s';
    Reactable.setMeta(tabEl, { selectedRow: rowInfo.index });
    
    if (rowInfo.values)
        Reactable.setFilter('compoundsTab', 'group', rowInfo.values.group);
    %s;
}", id, if (!is.null(onClick)) paste0("(", onClick, ")(tabEl, rowInfo, column);") else ""))
    
    bgstyle <- htmlwidgets::JS(sprintf("function(rowInfo, column, state)
{
    let ret = { }
    const gby = state.groupBy;
    debugger;
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
}", paste0("'", getFeatColGrpStartCols(groupDefs), "'", collapse = ","), getFeatColSepStyle()))
    
    colDefs <- lapply(colDefs, function(cd)
    {
        cd$style <- bgstyle
        return(cd)
    })

    # sync column order    
    tab <- copy(tab)
    setcolorder(tab, unlist(lapply(groupDefs, "[[", "columns")))
    
    rt <- reactable::reactable(tab, elementId = id, pagination = FALSE, wrap = FALSE, resizable = TRUE,
                               highlight = TRUE, onClick = oc, defaultExpanded = TRUE, columns = colDefs,
                               defaultColDef = reactable::colDef(style = bgstyle),
                               columnGroups = groupDefs, rowStyle = htmlwidgets::JS("function(rowInfo, state)
{
    const sel = state.meta.selectedRow;
    let ret = { cursor: 'pointer' };
    if (sel != null && rowInfo.index === sel)
        ret.background = '#eee';
    return ret;
}"), meta = list(selectedRow = NULL, plots = plots), ...)
    
    if (!visible)
        rt <- htmlwidgets::onRender(rt, htmlwidgets::JS("function(el, x) { el.style.display = 'none'; }"))
    
    return(rt)
}

reportHTMLGenerator$methods(
    genFeatTablePlain = function()
    {
        tab <- getFeatTable(objects$fGroups, ",")
        groupDefs <- getFeatGroupDefs(tab, NULL, replicateGroups(objects$fGroups))
        colDefs <- getFeatColDefs(tab, groupDefs, EICsTopMost)
        makeFeatReactable(tab, "detailsTabPlain", colDefs = colDefs, groupDefs = groupDefs, visible = TRUE,
                          plots = plots)
    },
    genFeatTableSuspects = function()
    {
        tab <- getFeatTable(objects$fGroups, NULL)
        groupDefs <- getFeatGroupDefs(tab, "susp_name", replicateGroups(objects$fGroups))
        colDefs <- getFeatColDefs(tab, groupDefs, EICsTopMost)
        makeFeatReactable(tab, "detailsTabSuspects", colDefs = colDefs, groupDefs = groupDefs, visible = TRUE,
                          plots = plots, groupBy = "susp_name")
    },
    genFeatTableComponents = function()
    {
        # tab <- getFeatTable(objects$fGroups, ",")
        # makeFeatReactable(tab, "detailsTabComponents", getFeatColDefs(tab, NULL), FALSE, plots)
        reactable::reactable(as.data.table(objects$fGroups), elementId = "detailsTabComponents")
    },
    genFeatTableTPs = function()
    {
        # UNDONE: put general suspect info (formula/neutralMass) not in group rows but in suspect aggregate (or simply omit?)
        
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
                                                       paste0("susp_", c("estIDLevel", "neutralMass", "formula", "d_rt",
                                                                         "d_mz", "sets", "InChIKey")))),
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
        
        colDefs <- getFeatColDefs(tabTPs, groupDefs, EICsTopMost)
        
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
    const chromEl = document.getElementById('chrom_view');
    const structEl = document.getElementById('struct_view');
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
    
        makeFeatReactable(tabTPs, "detailsTabTPs", FALSE, plots, groupBy = groupBy, colDefs = colDefs,
                          groupDefs = groupDefs, onClick = onClick)
    },
    
    genCompoundTable = function()
    {
        tab <- as.data.table(objects$compounds)[, c("group", "compoundName", "InChIKey")]
        tab[, structure := plotImg(plots$structs[InChIKey])][, InChIKey := NULL]
        tab[, spectrum := plotImg(plots$compounds[[group]]$spectra), by = "group"]
        tab[, scorings := plotImg(plots$compounds[[group]]$scores), by = "group"]
        return(reactable::reactable(tab, elementId = "compoundsTab", resizable = TRUE, bordered = TRUE,
                                    pagination = FALSE, columns = list(
            group = reactable::colDef(show = FALSE),
            structure = reactable::colDef(html = TRUE, minWidth = 150),
            spectrum = reactable::colDef(html = TRUE, minWidth = 600),
            scorings = reactable::colDef(html = TRUE, minWidth = 500)
        )))
    }
)

