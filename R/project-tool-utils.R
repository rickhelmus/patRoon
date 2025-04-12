textNote <- function(txt) div(style = "margin: 8px 0 12px; font-size: small", txt)

makeNewProjectHOT <- function(...)
{
    # NOTE: disable column sorting so we don't have to worry about getting correct row index
    # (https://github.com/jrowen/rhandsontable/issues/257)
    # NOTE: preventOverflow should not be set when columnSorting = FALSE
    # https://github.com/handsontable/handsontable/issues/4303
    # NOTE: set selectionMode to range as only row series can currently be queried
    # (https://github.com/jrowen/rhandsontable/issues/313)
    hotOpts <- list(rowHeaderWidth = 40, readOnly = TRUE,
                    sortIndicator = TRUE, selectCallback = TRUE,
                    currentRowClassName = "currentRow", stretchH = "all",
                    selectionMode = "range", outsideClickDeselects = FALSE,
                    contextMenu = FALSE, manualColumnResize = TRUE, ...)
    do.call(rhandsontable::rhandsontable, hotOpts)
}

doObserveSelDir <- function(input, session, textID, buttonID)
{
    observeEvent(input[[buttonID]], {
        d <- rstudioapi::selectDirectory("Select directory", path = input[[textID]])
        if (!is.null(d))
            updateTextInput(session, textID, value = d)
    })
}

moveSelectedTabRows <- function(dir, sel, tab)
{
    # NOTE: assume selection is a block range
    mv <- if (dir == "up")
        seq(min(sel) - 1, max(sel))
    else
        seq(min(sel), max(sel) + 1)
    
    tab <- copy(tab)
    tab[, index := .I]
    tab[mv, index := shift(index, n = if (dir == "up") - 1L else 1L, type = "cyclic")]
    tab <- tab[(index), -"index"]
    return(tab)
}

moveHOTSel <- function(session, dir, sel, id)
{
    newRange <- c(min(sel), max(sel)) + if (dir == "up") -1 else 1
    session$sendCustomMessage("selectHOTRows", list(tab = id, range = newRange))
}

exportShinyOutputVal <- function(output, var, r)
{
    # make reactive available for conditionalPanels: https://stackoverflow.com/a/38899895
    output[[var]] <- r
    outputOptions(output, var, suspendWhenHidden = FALSE)
    return(output)
}

fileSelect <- function(idText, idButton, label, value = "", ...)
{
    fillRow(
        flex = c(1, NA),
        textInput(idText, label, value, width = "100%", ...),
        actionButton(idButton, "", icon("folder-open"), style = "margin: 25px 0 0 10px;")
    )
}

selectSuspList <- function(session, inputName)
{
    sl <- rstudioapi::selectFile("Select suspect list", filter = "csv files (*.csv)")
    if (!is.null(sl))
    {
        csvTab <- tryCatch(fread(sl), error = function(e) FALSE, warning = function(w) FALSE)
        cols <- names(csvTab)
        massCols <- c("mz", "neutralMass", "formula", "SMILES", "InChI")
        
        err <- NULL
        if (is.logical(csvTab))
            err <- "Failed to open/parse selected csv file!"
        else if (nrow(csvTab) == 0)
            err <- "The selected files seems to be empty."
        else if (!"name" %in% cols)
            err <- "The selected file does not have a name column"
        else if (!any(massCols %in% cols))
            err <- paste("The selected file should have at least one of the columns:",
                         paste(massCols, collapse = ", "))
        
        if (!is.null(err))                
            rstudioapi::showDialog("Error", err, "")
        else
            updateTextInput(session, inputName, value = sl)
    }
}
