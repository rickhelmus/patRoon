#' @include main.R
#' @include report-html.R
NULL

getFeatTable <- function(fGroups, colSusp)
{
    tab <- as.data.table(fGroups, qualities = "score", average = TRUE, collapseSuspects = colSusp)
    
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

makeFeatReactable <- function(tab, id, visible, ...)
{
    rt <- reactable::reactable(tab, elementId = id, pagination = FALSE, wrap = FALSE, resizable = TRUE,
                               onClick = htmlwidgets::JS(sprintf("function(rowInfo, column)
{
    Reactable.setMeta('%s', { selectedRow: rowInfo.index });
}", id)), rowStyle = htmlwidgets::JS("function(rowInfo, state)
{
    const sel = state.meta.selectedRow;
    if (sel != null && rowInfo.index == sel)
        return { background: 'grey', cursor: 'pointer' };
    return { cursor: 'pointer' };
}"), meta = list(selectedRow = NULL), ...)
    
    if (!visible)
        rt <- htmlwidgets::onRender(rt, htmlwidgets::JS("function(el, x) { el.style.display = 'none'; }"))
    
    return(rt)
}


reportHTMLGenerator$methods(
    genFeatTablePlain = function()
    {
        makeFeatReactable(getFeatTable(objects$fGroups, ","), "detailsTabPlain", TRUE)
    },
    genFeatTableSuspects = function()
    {
        makeFeatReactable(getFeatTable(objects$fGroups, ","), "detailsTabSuspects", TRUE)
    },
    genFeatTableComponents = function()
    {
        makeFeatReactable(getFeatTable(objects$fGroups, ","), "detailsTabComponents", TRUE)
    },
    genFeatTableTPs = function()
    {
        tabTPsFeat <- getFeatTable(objects$fGroups, NULL)
        tabTPsFeat <- subsetDTColumnsIfPresent(tabTPsFeat, c("group", "ret", "mz", replicateGroups(objects$fGroups),
                                                             "adduct", "neutralMass",
                                                             paste0("susp_", c("name", "estIDLevel", "neutralMass",
                                                                               "formula", "d_rt", "d_mz", "sets"))))
        tabCompon <- as.data.table(objects$components)
        tabCompon <- subsetDTColumnsIfPresent(tabCompon, c("name", "parent_name", "parent_group", "group", "TP_retDir",
                                                           "TP_name", "retDir", "retDiff", "mzDiff", "formulaDiff",
                                                           "specSimilarity"))

        tabTPs <- merge(tabCompon, tabTPsFeat, by.x = c("group", "TP_name"), by.y = c("group", "susp_name"))
        setnames(tabTPs, "name", "component")
        setnames(tabTPs, "TP_name", "suspect")
        setnames(tabTPs, "parent_name", "parent_suspect")

        parAggr <- function(parentCol) htmlwidgets::JS(sprintf("function(values, rows, groupRows)
{
    return '<i>' + rows[0]['%s'] + '</i>';
}", parentCol))
        
        parAggred <- function(level) htmlwidgets::JS(sprintf("function(cellInfo, state)
{
    return (cellInfo.level === %d) ? cellInfo.value : '';
}", level))

        # NOTE: below values may be in components but then from suspect list
        # UNDONE: convert RT to minutes if needed
        tabTPs[, c("parent_ret", "parent_mz") := groupInfo(objects$fGroups)[parent_group, ]]
        
        for (col in c("parent_ret", "retDiff"))
            tabTPs[, (col) := round(get(col), 2)]
        for (col in c("parent_mz", "mzDiff"))
            tabTPs[, (col) := round(get(col), 5)]

        # add parent intensities & screening info
        rgs <- replicateGroups(objects$fGroups)
        tabTPsPar <- unique(subsetDTColumnsIfPresent(tabTPsFeat, c("group", rgs,
                                                                   paste0("susp_", c("estIDLevel", "neutralMass",
                                                                                     "formula", "d_rt", "d_mz", "sets")))),
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
        hasIDL <- !is.null(tabTPs[["susp_estIDLevel"]])
        if (hasSusp)
        {
            colDefs$susp_neutralMass$name <- "neutralMass"
            if (!is.null(tabTPs[["susp_formula"]]))
                colDefs$susp_formula$name <- "formula"
            colDefs$susp_d_rt$name <- "\U0394 ret"
            colDefs$susp_d_mz$name <- "\U0394 mz"
            if (hasIDL)
            {
                colDefs$susp_estIDLevel$name <- "estIDLevel" # UNDONE: hove-over?
                colDefs$susp_estIDLevel$align <- "right"
            }
        }
        
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
                                                            "specSimilarity"), names(tabTPs))),
            if (hasSusp) reactable::colGroup("screening",
                                             columns = intersect(c("susp_neutralMass", "susp_formula", "susp_d_rt",
                                                                   "susp_d_mz", "susp_estIDLevel"),
                                                                 names(tabTPs))) else NULL,
            reactable::colGroup("intensity", columns = rgs)
        ))
        
        setcolorder(tabTPs, unlist(lapply(colGroups, "[[", "columns")))
        
        makeFeatReactable(tabTPs, "detailsTabTPs", FALSE, groupBy = c("component", "suspect"), columns = colDefs,
                          columnGroups = colGroups, bordered = TRUE)
    }
)

