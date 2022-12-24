#' @include main.R
#' @include report-html.R
NULL

getFeatTable <- function(fGroups, colSusp)
{
    tab <- if (isScreening(fGroups))
        as.data.table(fGroups, qualities = "score", average = TRUE, collapseSuspects = colSusp)
    else
        as.data.table(fGroups, qualities = "score", average = TRUE)
    
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

makeFeatReactable <- function(tab, id, visible, plots, ..., onClick = NULL)
{
    oc <- htmlwidgets::JS(sprintf("function(rowInfo, column)
{
    const tabEl = '%s';
    Reactable.setMeta(tabEl, { selectedRow: rowInfo.index });
    %s;
}", id, if (!is.null(onClick)) paste0("(", onClick, ")(tabEl, rowInfo, column);") else ""))
    
    rt <- reactable::reactable(tab, elementId = id, pagination = FALSE, wrap = FALSE, resizable = TRUE,
                               onClick = oc, defaultExpanded = TRUE,
                               rowStyle = htmlwidgets::JS("function(rowInfo, state)
{
    const sel = state.meta.selectedRow;
    let ret = { cursor: 'pointer' };
    if (sel != null && rowInfo.index === sel)
        ret.background = 'grey';
    else if (rowInfo.level === 0)
        ret.borderTop = '2px solid black';
    else if (rowInfo.level === 1)
        ret.borderBottom = '1px solid black';
    return ret;
}"), meta = list(selectedRow = NULL, plots = plots), ...)
    
    if (!visible)
        rt <- htmlwidgets::onRender(rt, htmlwidgets::JS("function(el, x) { el.style.display = 'none'; }"))
    
    return(rt)
}


reportHTMLGenerator$methods(
    genFeatTablePlain = function()
    {
        makeFeatReactable(getFeatTable(objects$fGroups, ","), "detailsTabPlain", TRUE, plots)
    },
    genFeatTableSuspects = function()
    {
        makeFeatReactable(getFeatTable(objects$fGroups, ","), "detailsTabSuspects", TRUE, plots)
    },
    genFeatTableComponents = function()
    {
        makeFeatReactable(getFeatTable(objects$fGroups, ","), "detailsTabComponents", TRUE, plots)
    },
    genFeatTableTPs = function()
    {
        # UNDONE: put general suspect info (formula/neutralMass) not in group rows but in suspect aggregate (or simply omit?)
        
        tabTPsFeat <- getFeatTable(objects$fGroups, NULL)
        tabTPsFeat <- subsetDTColumnsIfPresent(tabTPsFeat, c("group", "ret", "mz", replicateGroups(objects$fGroups),
                                                             "adduct", "neutralMass",
                                                             paste0("susp_", c("name", "estIDLevel", "neutralMass",
                                                                               "formula", "d_rt", "d_mz", "sets",
                                                                               "InChIKey"))))
        tabCompon <- as.data.table(objects$components)
        tabCompon <- subsetDTColumnsIfPresent(tabCompon, c("name", "parent_name", "parent_group", "group", "TP_retDir",
                                                           "TP_name", "retDir", "retDiff", "mzDiff", "formulaDiff",
                                                           "specSimilarity", "mergedBy"))

        tabTPs <- merge(tabCompon, tabTPsFeat, by.x = c("group", "TP_name"), by.y = c("group", "susp_name"))
        setnames(tabTPs, "name", "component")
        setnames(tabTPs, "TP_name", "suspect")
        setnames(tabTPs, "parent_name", "parent_suspect")

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
        
        colDefs <- list()
        # set parent 'aggregates': actual value of parent feature group
        for (col in grep("^parent_", names(tabTPs), value = TRUE))
        {
            colTP <- sub("^parent_", "", col)
            colDefs[[colTP]] <- reactable::colDef(aggregate = parAggr(col), aggregated = parAggred(0), html = TRUE)
            colDefs[[col]] <- reactable::colDef(show = FALSE)
        }
        # similar for TP retDir
        colDefs$retDir <- reactable::colDef(aggregate = parAggr("TP_retDir"), aggregated = parAggred(1), html = TRUE)
        colDefs$TP_retDir <- reactable::colDef(show = FALSE)
        
        # these are grouped
        colDefs$retDiff <- reactable::colDef(name = "\U0394 ret")
        colDefs$mzDiff <- reactable::colDef(name = "\U0394 mz")
        colDefs$formulaDiff <- reactable::colDef(name = "\U0394 formula")
        
        # and these...
        hasSusp <- !is.null(tabTPs[["susp_neutralMass"]])
        if (hasSusp)
        {
            colDefs$susp_neutralMass$name <- "neutralMass"
            if (!is.null(tabTPs[["susp_formula"]]))
                colDefs$susp_formula$name <- "formula"
            colDefs$susp_d_rt$name <- "\U0394 ret"
            colDefs$susp_d_mz$name <- "\U0394 mz"
            if (!is.null(tabTPs[["susp_estIDLevel"]]))
            {
                colDefs$susp_estIDLevel$name <- "estIDLevel" # UNDONE: tool-tip?
                colDefs$susp_estIDLevel$align <- "right"
            }
        }
        
        # InChIKeys are only there for internal usage
        if (!is.null(tabTPs[["susp_InChIKey"]]))
            colDefs$susp_InChIKey <- reactable::colDef(show = FALSE)
        if (!is.null(tabTPs[["parent_susp_InChIKey"]]))
            colDefs$parent_susp_InChIKey <- reactable::colDef(show = FALSE)
        
        colDefs$group$cell <- function(value, index)
        {
            htmltools::div(value,
                           htmltools::br(),
                           sparkline::sparkline(EICsTopMost[[value]]$intensity, xvalues = EICsTopMost[[value]]$time,
                                                type = "line"))
        }
        
        ststyle <- list(borderRight = "3px solid #eee") # from Reactable examples
        colDefs$suspect$style <- colDefs$suspect$headerStyle <- ststyle
        
        colGroups <- pruneList(list(
            # workaround for stickies: https://github.com/glin/reactable/issues/236#issuecomment-1107911895
            reactable::colGroup("", columns = c("component", "suspect"), sticky = "left",
                                headerStyle = ststyle),
            reactable::colGroup("feature", columns = c("group", "ret", "mz")),
            reactable::colGroup("TP", columns = intersect(c("retDiff", "mzDiff", "formulaDiff", "retDir",
                                                            "specSimilarity", "mergedBy"), names(tabTPs))),
            if (hasSusp) reactable::colGroup("screening",
                                             columns = intersect(c("susp_neutralMass", "susp_formula", "susp_d_rt",
                                                                   "susp_d_mz", "susp_estIDLevel"),
                                                                 names(tabTPs))) else NULL,
            reactable::colGroup("intensity", columns = rgs)
        ))
        
        setcolorder(tabTPs, unlist(lapply(colGroups, "[[", "columns")))
        
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
}"
        
        makeFeatReactable(tabTPs, "detailsTabTPs", FALSE, plots, groupBy = c("component", "suspect"), columns = colDefs,
                          columnGroups = colGroups, bordered = TRUE, onClick = onClick)
    }
)

