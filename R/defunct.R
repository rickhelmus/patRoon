#' Defunct functions.
#'
#' These functions do not work anymore and may be updated in the future.
#' future.
#'
#' @name patRoon-defunct
#' @keywords internal
NULL

# nocov start
getCompInfoText <- function(compResults, compIndex, addHTMLURL, normalizeScores, mCompNames)
{
    columns <- names(compResults)
    
    if (normalizeScores)
        compResults <- normalizeCompScores(compResults, mCompNames, FALSE) # UNDONE: select normalization method? Add scoreranges
    
    resultRow <- compResults[compIndex, ]
    
    addValText <- function(curText, fmt, cols)
    {
        cols <- getAllCompCols(cols, columns, mCompNames)
        ret <- ""
        for (cl in cols)
        {
            if (!is.null(resultRow[[cl]]) && !is.na(resultRow[[cl]]) &&
                (!is.character(resultRow[[cl]]) || nzchar(resultRow[[cl]])))
            {
                fm <- sprintf("%s: %s\n", cl, fmt)
                ret <- paste0(ret, sprintf(fm, resultRow[[cl]]))
            }
        }
        
        return(paste0(curText, ret))
    }
    
    ctext <- ""
    
    if (addHTMLURL)
    {
        addIdURL <- function(param, ident, db)
        {
            ident <- as.character(ident)
            
            if (is.na(ident) || !nzchar(ident))
                return("")
            
            # CSI:FingerID might return multiple identifiers
            idlist <- unlist(strsplit(ident, ";"))
            
            if (grepl("pubchem", tolower(db)))
                fmt <- "<a href=\"https://pubchem.ncbi.nlm.nih.gov/compound/%s\">%s</a>"
            else if (tolower(db) == "chemspider")
                fmt <- "<a href=\"http://www.chemspider.com/Search.aspx?q=%s\">%s</a>"
            else
                fmt <- "%s"
            
            return(sprintf("%s: %s\n", param, paste0(sprintf(fmt, idlist, idlist), collapse = "; ")))
        }
        
        if (!is.null(resultRow$identifier)) # compounds were not merged, can use 'regular' column
            ctext <- paste0(ctext, addIdURL("identifier", resultRow$identifier, resultRow$database))
        else
        {
            idcols <- getAllCompCols("identifier", columns, mCompNames)
            dbcols <- getAllCompCols("database", columns, mCompNames)
            
            if (allSame(resultRow[, idcols, with = FALSE])) # no need to show double ids
                ctext <- paste0(ctext, addIdURL("identifier", resultRow[[idcols[1]]], resultRow[[dbcols[1]]]))
            else
            {
                for (i in seq_along(idcols))
                    ctext <- paste0(ctext, addIdURL(idcols[i], resultRow[[idcols[i]]], resultRow[[dbcols[i]]]))
            }
        }
    }
    else
        ctext <- addValText(ctext, "%s", "identifier")
    
    ctext <- addValText(ctext, "%s", c("compoundName", "formula", "SMILES"))
    
    if (length(getAllCompCols("InChIKey", columns, mCompNames)) > 0)
        ctext <- addValText(ctext, "%s", "InChIKey")
    else # only add InChIKey1/2 if full isn't available
        ctext <- addValText(ctext, "%s", c("InChIKey1", "InChIKey2"))
    
    ctext <- addValText(ctext, "%.2f", getCompScoreColNames())
    ctext <- addValText(ctext, "%.2f", c("XlogP", "AlogP"))
    
    # remove trailing newline
    ctext <- gsub("\n$", "", ctext)
    
    return(ctext)
}

getCompViewerUI <- function(pageChoices)
{
    fillPage(
        title = "Compound identification results",
        tags$head(includeScript(system.file("js", "utils-comp.js", package = "patRoon"))),
        
        fillCol(
            flex = c(NA, NA, NA, 1),
            
            fillRow(
                height = 40,
                flex = c(1, 2),
                
                fillRow(
                    width = 95,
                    actionButton("previousGroup", "", icon("arrow-left"), onclick = "selectPrevHot(\"group\");"),
                    actionButton("nextGroup", "", icon("arrow-right"), onclick = "selectNextHot(\"group\");")
                ),
                
                fillRow(
                    flex = c(NA, NA, 1, NA),
                    
                    fillRow(
                        width = 95,
                        actionButton("previousResult", "", icon("arrow-left"), onclick = "selectPrevHot(\"comp\");"),
                        actionButton("nextResult", "", icon("arrow-right"), onclick = "selectNextHot(\"comp\");")
                    ),
                    
                    fillRow(
                        width = 200,
                        flex = c(1, 3),
                        br(),
                        checkboxInput("normalizeScores", "Normalize scores", TRUE)
                    ),
                    
                    fillRow(
                        br()
                    ),
                    
                    fillRow(
                        width = 130,
                        selectInput("page", NULL, pageChoices, width = 120)
                    )
                )
            ),
            
            fillRow(
                flex = c(1, 2),
                height = 405,
                
                fillCol(
                    div(
                        style = "border: 1px solid black; margin: 5px;",
                        rhandsontable::rHandsontableOutput("groupHot")
                    )
                ),
                
                fillCol(
                    div(
                        style = "border: 1px solid black; margin: 5px;",
                        rhandsontable::rHandsontableOutput("metFragHot")
                    )
                )
            ),
            
            hr(),
            
            fillRow(
                flex = c(2, NA, 1),
                
                fillCol(
                    plotly::plotlyOutput("spectrum", height = "100%")
                ),
                
                fillCol(
                    br()
                ),
                
                fillCol(
                    wellPanel(
                        style = "overflow-y: auto; height: 100%; overflow-x: auto;",
                        
                        plotOutput("structure", width = 200, height = 200),
                        tags$head(tags$style(type="text/css", "#structure { width: 60%; display: block; margin-left: auto; margin-right: auto; }")),
                        
                        htmlOutput("compInfo"),
                        tags$head(tags$style(type="text/css", "#compInfo { border: 1px solid black; border-style: dotted; margin: 5px; padding: 5px; overflow-x: auto; white-space: pre; }"))
                    )
                )
            )
        )
    )
}

# @rdname GUI-utils
#' @details The \code{compoundViewer} method is used to view compound
#'   identification results. It will display available candidate information
#'   such as scorings and identifiers, MS/MS spectra with explained peaks and
#'   chemical structures.
#'
#' @param MSPeakLists A \code{\link{MSPeakLists}} object.
#' @param compounds A \code{\link{compounds}} object.
#'
#' @rdname patRoon-defunct
#' @aliases compoundViewer
#' @keywords internal
#' @export
setMethod("compoundViewer", c("featureGroups", "MSPeakLists", "compounds"), function(fGroups, MSPeakLists, compounds)
{
    .Defunct("reportHTML")
    return(NULL)
    
    compTable <- compoundTable(compounds)
    
    fGroups <- fGroups[, groupNames(compounds)]
    
    fTable <- featureTable(fGroups)
    anaInfo <- analysisInfo(fGroups)
    gInfo <- groupInfo(fGroups)
    gTable <- groupTable(fGroups)
    avgGTable <- averageGroups(fGroups)
    ftindex <- groupFeatIndex(fGroups)
    pLists <- peakLists(MSPeakLists)
    
    hotOpts <- list(rowHeaderWidth = 40, readOnly = TRUE, disableVisualSelection = "area",
                    columnSorting = TRUE, sortIndicator = TRUE, selectCallback = TRUE,
                    currentRowClassName = "currentRow", stretchH = "all",
                    preventOverflow = "horizontal", multiSelect = FALSE,
                    outsideClickDeselects = FALSE, manualColumnResize = TRUE,
                    contextMenu = FALSE)
    
    for (g in seq_along(compTable))
        compTable[[g]]$structure <- NA
    
    addResourcePath("struct", tempdir(TRUE))
    
    getPageRange <- function(group, page, resultsPerPage)
    {
        totalr <- nrow(compTable[[group]])
        starti <- ((page - 1) * resultsPerPage) + 1
        endi <- starti + min(resultsPerPage-1, totalr - starti)
        return(c(starti, endi))
    }
    
    generatePageLabels <- function(group, resultsPerPage)
    {
        totalr <- nrow(compTable[[group]])
        pagen <- round(totalr / resultsPerPage + 0.5)
        n <- sapply(seq_len(pagen), function(p)
        {
            pr <- getPageRange(group, p, resultsPerPage)
            return(sprintf("%d - %d", pr[1], pr[2]))
        })
        
        ret <- seq_len(pagen)
        names(ret) <- n
        return(ret)
    }
    
    server <- function(input, output, session)
    {
        rValues <- reactiveValues(currentFGroup = rownames(gInfo)[1],
                                  currentResult = 1,
                                  currentPage = 1,
                                  resultsPerPage = 50)
        
        compData <- reactive({
            pr <- getPageRange(rValues$currentFGroup, rValues$currentPage, rValues$resultsPerPage)
            
            ct <- compTable[[rValues$currentFGroup]]
            if (input$normalizeScores)
                ct <- normalizeCompScores(ct, mergedCompoundNames(compounds), FALSE) # UNDONE: select normalization method? Add scoreranges
            
            ret <- ct[pr[1]:pr[2], ]
            
            # generate structures if necessary
            if (any(is.na(ret$structure)))
            {
                for (i in which(is.na(ret$structure)))
                {
                    f <- tempfile(sprintf("%d-%d-", match(rValues$currentFGroup, colnames(gTable)), i),
                                  fileext = ".png")
                    ret$structure[i] <- paste0("struct/", basename(f)) # server location
                    png(f, 100, 100)
                    par(mai = rep(0, 4))
                    mol <- getMoleculesFromSMILES(ret$SMILES[i])
                    if (isValidMol(mol[[1]])) # this may fail
                        plot(getRCDKStructurePlot(mol[[1]], width = 100, height = 100))
                    dev.off()
                }
            }
            
            allCols <- names(ret)
            getCols <- function(col)
            {
                ret <- getAllCompCols(col, allCols, mergedCompoundNames(compounds))
                return(ret[sapply(ret, function(cl) !all(is.na(cl)))])
            }
            
            keepCols <- c(getCols("identifier"), "structure",
                          getCols(c("compoundName", "formula", "explainedPeaks")))
            
            # optional scoring columns
            scorecols <- getCols(getCompScoreColNames())
            keepCols <- c(keepCols, scorecols)
            
            return(ret[, keepCols, with = FALSE])
        })
        
        fragmentInfo <- reactive({
            compTable[[rValues$currentFGroup]]$fragInfo[[rValues$currentResult]]
        })
        
        observeEvent(input$page, {
            rValues$currentPage <- as.numeric(input$page)
            pr <- getPageRange(rValues$currentFGroup, rValues$currentPage, rValues$resultsPerPage)
            rValues$currentResult <- rValues$currentResult %% rValues$resultsPerPage + pr[1] - 1
        })
        
        observeEvent(input$groupHot_select$select$r, {
            rValues$currentFGroup <- rownames(gInfo)[input$groupHot_select$select$rAll[1]]
            rValues$currentResult <- 1
            updateSelectInput(session, "page", choices = generatePageLabels(rValues$currentFGroup, rValues$resultsPerPage))
        })
        
        observeEvent(input$metFragHot_select$select$r, {
            pr <- getPageRange(rValues$currentFGroup, rValues$currentPage, rValues$resultsPerPage)
            rValues$currentResult <- input$metFragHot_select$select$rAll[1] + pr[1] - 1
        })
        
        output$groupHot <- rhandsontable::renderRHandsontable({
            gData <- data.frame(group = rownames(gInfo), retention = gInfo$rts, mz = gInfo$mzs, stringsAsFactors = FALSE)
            gData[, unique(anaInfo$group)] <- t(avgGTable)
            
            hot <- do.call(rhandsontable::rhandsontable,
                           c(list(gData, colHeaders = c("Feature group", "Retention", "m/z", unique(anaInfo$group)),
                                  width = NULL, height = 400, maxRows = nrow(gData)),
                             hotOpts))
            
            return(hot)
        })
        
        output$metFragHot <- rhandsontable::renderRHandsontable({
            cd <- compData()
            hot <- do.call(rhandsontable::rhandsontable, c(list(cd, width = NULL, height = 400, maxRows = nrow(cd),
                                                                rowHeaders = seq_len(nrow(cd))), hotOpts)) %>%
                rhandsontable::hot_rows(rowHeights = 100) %>%
                rhandsontable::hot_col("structure", width = 110, renderer = "
                        function(instance, td, row, col, prop, value, cellProperties)
                        {
                            var escaped = Handsontable.helper.stringify(value),
                            img = document.createElement('IMG');
                            img.src = value;

                            Handsontable.Dom.addEvent(img, 'mousedown', function (e)
                            {
                                e.preventDefault(); // prevent selection quirk
                            });

                            Handsontable.Dom.empty(td);
                            td.appendChild(img);

                            return td;
                        }")
            return(hot)
        })
        
        output$spectrum <- plotly::renderPlotly({
            spec <- pLists[[rValues$currentFGroup]][["MSMS"]]
            fi <- fragmentInfo()
            
            p <- plotly::plot_ly(type="scatter", mode = "lines", hoverinfo = "text") %>%
                plotly::config(displaylogo = FALSE, scrollZoom = TRUE,
                               modeBarButtonsToRemove = c("hoverClosestCartesian", "hoverCompareCartesian"))
            
            for (i in seq_len(nrow(spec)))
            {
                infoi <- if (is.null(fi) || nrow(fi) == 0) FALSE else match(i, fi$PLIndex, nomatch = FALSE)
                
                htext <- ""
                if (infoi)
                {
                    htext <- sprintf("m/z: %f<br>int.: %.0f<br>formula: %s",
                                     spec$mz[i], spec$intensity[i], fi$formula[infoi])
                    
                    if (!is.null(fi$score) && !is.na(fi$score[infoi]))
                        htext <- paste0(htext, sprintf("<br>score: %f", fi$score[infoi]))
                    if (!is.null(fi$mergedBy))
                        htext <- paste0(htext, sprintf("<br>merged by: %s", fi$mergedBy[[infoi]]))
                }
                
                p <- plotly::add_trace(p, x = rep(spec$mz[i], 2), y = c(0, spec$intensity[i]),
                                       line = list(color = if (infoi) "blue" else "grey",
                                                   width = if (infoi) 3 else 1),
                                       hoverinfo = if (infoi) "text" else "none",
                                       text = htext)
            }
            
            p <- plotly::layout(p, showlegend = FALSE, dragmode = "pan", plot_bgcolor = "#F5FFFA",
                                margin = list(t = 0),
                                xaxis = list(title = "m/z", range = range(spec$mz) * c(0.9, 1.1)),
                                yaxis = list(title = "Intensity", exponentformat = "E", range = c(0, max(spec$intensity) * 1.1)))
            
            # showarrow (like shapes) crash RStudio :(
            if (!is.null(fi) && nrow(fi) > 0)
                p <- plotly::add_annotations(p, x = fi$mz, y = fi$intensity + (max(spec$intensity) * 0.02), xref = "x", yref = "y",
                                             text = fi$formula, showarrow = FALSE) # = TRUE, arrowhead = 7)
            
            # workaround for annoying shiny warnings: https://github.com/ropensci/plotly/issues/985
            p$elementId <- NULL
            
            return(p)
        })
        
        output$structure <- renderPlot({
            mol <- getMoleculesFromSMILES(compTable[[rValues$currentFGroup]]$SMILES[rValues$currentResult])
            if (isValidMol(mol[[1]]))
                plot(getRCDKStructurePlot(mol[[1]], width = 200, height = 200))
        })
        
        output$compInfo <- renderText({
            return(getCompInfoText(compTable[[rValues$currentFGroup]], rValues$currentResult, TRUE, input$normalizeScores,
                                   mergedCompoundNames(compounds)))
        })
    }
    
    runApp(shinyApp(getCompViewerUI(generatePageLabels(names(fGroups)[1], 50)), server))
})

#' @details \code{groupFeaturesScreening} was replaced by
#'   \code{\link{screenSuspects}}. Please refer to this function and the
#'   handbook for the updated interface.
#'
#' @param ... Ignore.
#'
#' @rdname patRoon-defunct
#' @keywords internal
#' @export
groupFeaturesScreening <- function(...) stop("This function was replaced by screenSuspects(). ",
                                             "Please see ?screenSuspects and the handbook for the updated interface.")

#' @details \code{checkChromatograms} was replaced by
#'   \code{\link{checkFeatures}}. Please refer to this function and the
#'   handbook for the updated interface.
#'
#' @rdname patRoon-defunct
#' @keywords internal
#' @export
checkChromatograms <- function(...) stop("This function was replaced by checkFeatures(). ",
                                         "Please see ?checkFeatures and the handbook for the updated interface.")

# nocov end
