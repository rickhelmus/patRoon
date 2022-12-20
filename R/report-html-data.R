#' @include main.R
#' @include report-html.R
NULL

getFeatTable <- function(fGroups)
{
    tab <- as.data.table(fGroups, qualities = "score", average = TRUE)
    
    # if (rmdVars$retMin) UNDONE
    tab[, ret := ret / 60]
    
    for (col in names(tab)[(sapply(tab, is.numeric))])
        set(tab, j = col, value = round(tab[[col]], if (col == "mz") 5 else 2))
    
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

genHTMLReportDataFeatures <- function(fGroups, EICs)
{
    ret <- list()
    
    ret$features <- list()
    ret$features$plain <- makeReactable(getFeatTable(), "detailsTabPlain", TRUE)
    ret$features$suspects <- makeReactable(getFeatTable(), "detailsTabSuspects", TRUE)
    ret$features$components <- makeReactable(getFeatTable(), "detailsTabComponents", TRUE)
    return(ret)
    
    tabTPs <- getFeatTable()
    tabCompon <- as.data.table(components)
    tabTPs <- merge(tabCompon[, c("group", setdiff(names(tabCompon), names(tabTPs))), with = FALSE],
                              tabTPs, by = "group")
    setnames(tabTPs, "name", "component")
    setnames(tabTPs, "TP_name", "name")
    
    colDefs <- list()
    for (col in c("group", "name", "formula"))
    {
        colDefs[[col]] <- reactable::colDef(aggregate = parAggr(col), html = TRUE)
        colDefs[[paste0("parent_", col)]] <- reactable::colDef(show = FALSE)
    }
    tabTPs[, c("parent_rt", "parent_mz", "parent_SMILES", "parent_InChI", "parent_InChIKey", "parent_neutralMass",
               "size", "SMILES", "InChI", "InChIKey", "links", "intensity") := NULL]
    
    tabTPsRT <- makeReactable(tabTPs, "detailsTabTPs", FALSE, groupBy = "component", columns = colDefs)
    ret$features$TPs <- tabTPsRT
    
    return(ret)
    
}

generateHTMLReportData <- function(fGroups, MSPeakLists, formulas, compounds, components, TPs, EICs, plots)
{
    return(genHTMLReportDataFeatures(fGroups, EICs))
}


reportHTMLGenerator$methods(
    genFeatTablePlain = function()
    {
        makeFeatReactable(getFeatTable(objects$fGroups), "detailsTabPlain", TRUE)
    },
    genFeatTableSuspects = function()
    {
        makeFeatReactable(getFeatTable(objects$fGroups), "detailsTabSuspects", TRUE)
    },
    genFeatTableComponents = function()
    {
        makeFeatReactable(getFeatTable(objects$fGroups), "detailsTabComponents", TRUE)
    },
    genFeatTableTPs = function()
    {
        tabTPs <- getFeatTable(objects$fGroups)
        tabCompon <- as.data.table(objects$components)
        tabTPs <- merge(tabCompon[, c("group", setdiff(names(tabCompon), names(tabTPs))), with = FALSE],
                        tabTPs, by = "group")
        setnames(tabTPs, "name", "component")
        setnames(tabTPs, "TP_name", "name")
        tabTPs[, c("parent_rt", "parent_mz", "parent_SMILES", "parent_InChI", "parent_InChIKey", "parent_neutralMass",
                   "size", "SMILES", "InChI", "InChIKey", "links", "intensity") := NULL]
        
        parAggr <- function(col) htmlwidgets::JS(sprintf("function(values, rows)
{
    return '<i>' + rows[0].parent_%s + '</i>';
}", col))
        
        colDefs <- list()
        for (col in c("group", "name", "formula"))
        {
            colDefs[[col]] <- reactable::colDef(aggregate = parAggr(col), html = TRUE)
            colDefs[[paste0("parent_", col)]] <- reactable::colDef(show = FALSE)
        }
        
        colDefs$group$cell <- function(value, index)
        {
            htmltools::div(value,
                           htmltools::br(),
                           sparkline::sparkline(EICsTopMost[[value]]$intensity, xvalues = EICsTopMost[[value]]$time,
                                                type = "line"))
        }
        
        makeFeatReactable(tabTPs, "detailsTabTPs", FALSE, groupBy = "component", columns = colDefs)
    }
)

