#' @include main.R
#' @include compounds.R
#' @include mspeaklists.R
#' @include formulas.R
NULL

# get a vector of all (merged) columns
getAllCompCols <- function(targetCols, allCols, mCompNames)
{
    targetCols <- c(targetCols, sapply(targetCols, function(cl) paste0(cl, "-", mCompNames)))
    return(targetCols[targetCols %in% allCols])
}

mergeFragInfo <- function(fiLeft, fiRight, leftName, rightName, setMergedBy)
{
    fiLeft <- copy(fiLeft)

    if (setMergedBy && is.null(fiLeft[["mergedBy"]]))
        fiLeft[, mergedBy := list(list(leftName))]

    # for overlap: just add label
    fiLeft[PLIndex %in% fiRight$PLIndex,
            mergedBy := lapply(mergedBy, function(cm) list(unique(c(cm, rightName))))]

    # add unique
    fiUnique <- fiRight[!PLIndex %in% fiLeft$PLIndex]
    if (setMergedBy)
        fiUnique[, mergedBy := list(list(rightName))]
    if (nrow(fiUnique) > 0)
        fiLeft <- rbind(fiLeft, fiUnique, fill = TRUE)

    return(fiLeft)
}

getCompScoreColNames <- function()
{
    return(c("score",
             "fragScore",
             "metFusionScore",
             "individualMoNAScore",
             "numberPatents",
             "pubMedReferences",
             "extReferenceCount",
             "dataSourceCount",
             "referenceCount",
             "RSCCount",
             "formulaScore",
             "smartsInclusionScore",
             "smartsExclusionScore",
             "suspectListScore",
             "retentionTimeScore"))
}

normalizeCompScores <- function(compResults, mCompNames)
{
    compResults <- copy(compResults)
    columns <- names(compResults)
    scoreCols <- getAllCompCols(getCompScoreColNames(), columns, mCompNames)
    compResults[, (scoreCols) := lapply(.SD, normalize), .SDcols = scoreCols]
    return(compResults)
}

getCompInfoList <- function(compResults, compIndex, addHTMLURL, mCompNames)
{
    columns <- names(compResults)

    resultRow <- compResults[compIndex, ]

    addValText <- function(curText, fmt, cols)
    {
        cols <- getAllCompCols(cols, columns, mCompNames)
        ret <- character()
        for (cl in cols)
        {
            if (!is.null(resultRow[[cl]]) && !is.na(resultRow[[cl]]) &&
                (!is.character(resultRow[[cl]]) || nzchar(resultRow[[cl]])))
            {
                fm <- sprintf("%s: %s", cl, fmt)
                ret <- c(ret, sprintf(fm, resultRow[[cl]]))
            }
        }

        return(c(curText, ret))
    }

    ctext <- character()

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
                fmt <- "<a target=\"_blank\" href=\"https://pubchem.ncbi.nlm.nih.gov/compound/%s\">%s</a>"
            else if (tolower(db) == "chemspider")
                fmt <- "<a target=\"_blank\" href=\"http://www.chemspider.com/Search.aspx?q=%s\">%s</a>"
            else
                fmt <- "%s"

            return(sprintf("%s: %s", param, paste0(sprintf(fmt, idlist, idlist), collapse = "; ")))
        }

        if (!is.null(resultRow$identifier)) # compounds were not merged, can use 'regular' column
            ctext <- c(ctext, addIdURL("identifier", resultRow$identifier, resultRow$database))
        else
        {
            idcols <- getAllCompCols("identifier", columns, mCompNames)
            dbcols <- getAllCompCols("database", columns, mCompNames)

            if (allSame(resultRow[, idcols, with = FALSE])) # no need to show double ids
                ctext <- c(ctext, addIdURL("identifier", resultRow[[idcols[1]]], resultRow[[dbcols[1]]]))
            else
            {
                for (i in seq_along(idcols))
                    ctext <- c(ctext, addIdURL(idcols[i], resultRow[[idcols[i]]], resultRow[[dbcols[i]]]))
            }
        }
    }
    else
        ctext <- addValText(ctext, "%s", "identifier")

    ctext <- addValText(ctext, "%s", c("trivialName", "formula", "SMILES", "analysis"))

    if (length(getAllCompCols("InChIKey", columns, mCompNames)) > 0)
        ctext <- addValText(ctext, "%s", "InChIKey")
    else # only add InChIKey1/2 if full isn't available
        ctext <- addValText(ctext, "%s", c("InChIKey1", "InChIKey2"))

    ctext <- addValText(ctext, "%.2f", c("XlogP", "AlogP"))

    return(ctext)
}

buildMFLandingURL <- function(mfSettings, peakList, precursorMz)
{
    # Via personal communication Steffen/Emma, see https://github.com/Treutler/MetFamily/blob/22b9f46b2716b805c24c03d260045605c0da8b3e/ClusteringMS2SpectraGUI.R#L2433
    # Code adopted from MetFamily R package: https://github.com/Treutler/MetFamily
    
    if (is.null(mfSettings))
    {
        # no settings given, simply default to PubChem
        mfSettings <- list(MetFragDatabaseType = "PubChem")
    }
    
    mfSettings$IonizedPrecursorMass <- precursorMz
    mfSettings$NeutralPrecursorMass <- "" # make sure user needs to calculate it and remove default
    
    PL <- paste(peakList$mz, peakList$intensity, sep = " ", collapse = "; ")
    mfSettings$PeakList <- PL
    
    if (mfSettings$MetFragDatabaseType == "ExtendedPubChem")
        mfSettings$MetFragDatabaseType <- "PubChem" # user should tick box for now...
    
    # Allowed parameters: list taken from error page when unsupported parameter is given
    mfSettings <- mfSettings[names(mfSettings) %in%
                                 c("FragmentPeakMatchAbsoluteMassDeviation", "FragmentPeakMatchRelativeMassDeviation",
                                   "DatabaseSearchRelativeMassDeviation", "PrecursorCompoundIDs", "IonizedPrecursorMass",
                                   "NeutralPrecursorMass", "NeutralPrecursorMolecularFormula", "PrecursorIonMode",
                                   "PeakList", "MetFragDatabaseType")]
    
    setstr <- paste0(paste0(names(mfSettings), "=", mfSettings), collapse = "&")
    ret <- paste0("https://msbi.ipb-halle.de/MetFragBeta/landing.xhtml?", setstr)
    ret <- sprintf("<a target=\"_blank\" href=\"%s\">MetFragWeb</a>", ret)
    
    return(ret)
}

# nocov start
getCompInfoText <- function(compResults, compIndex, addHTMLURL, normalizeScores, mCompNames)
{
    columns <- names(compResults)
    
    if (normalizeScores)
        compResults <- normalizeCompScores(compResults, mCompNames)
    
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
    
    ctext <- addValText(ctext, "%s", c("trivialName", "formula", "SMILES", "analysis"))
    
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
                        rHandsontableOutput("groupHot")
                    )
                ),

                fillCol(
                    div(
                        style = "border: 1px solid black; margin: 5px;",
                        rHandsontableOutput("metFragHot")
                    )
                )
            ),

            hr(),

            fillRow(
                flex = c(2, NA, 1),

                fillCol(
                    plotlyOutput("spectrum", height = "100%")
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

#' @details The \code{compoundViewer} method is used to view compound
#'   identification results. It will display available candidate information
#'   such as scorings and identifiers, MS/MS spectra with explained peaks and
#'   chemical structures.
#'
#' @param MSPeakLists A \code{\link{MSPeakLists}} object.
#' @param compounds A \code{\link{compounds}} object.
#'
#' @rdname GUI-utils
#' @aliases compoundViewer
#' @export
setMethod("compoundViewer", c("featureGroups", "MSPeakLists", "compounds"), function(fGroups, MSPeakLists, compounds)
{
    compTable <- compoundTable(compounds)

    fGroups <- groupNamesFilter(fGroups, "compoundViewer", names(compoundTable(compounds)), verbose = FALSE)

    fTable <- featureTable(fGroups)
    anaInfo <- analysisInfo(fGroups)
    gInfo <- groupInfo(fGroups)
    gTable <- groups(fGroups)
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

    addResourcePath("struct", tempdir())

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

        getAnalysis <- reactive({ return(compTable[[rValues$currentFGroup]]$analysis[1]) })

        compData <- reactive({
            pr <- getPageRange(rValues$currentFGroup, rValues$currentPage, rValues$resultsPerPage)

            ct <- compTable[[rValues$currentFGroup]]
            if (input$normalizeScores)
                ct <- normalizeCompScores(ct, mergedCompoundNames(compounds))

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
                    mol <- rcdk::parse.smiles(ret$SMILES[i])
                    if (!is.na(mol)) # this may fail
                        rcdkplot(mol[[1]], width = 100, height = 100)
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
                          getCols(c("trivialName", "formula", "explainedPeaks")))

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
            rValues$currentFGroup <- rownames(gInfo)[input$groupHot_select$select$r]
            rValues$currentResult <- 1
            updateSelectInput(session, "page", choices = generatePageLabels(rValues$currentFGroup, rValues$resultsPerPage))
        })

        observeEvent(input$metFragHot_select$select$r, {
            pr <- getPageRange(rValues$currentFGroup, rValues$currentPage, rValues$resultsPerPage)
            rValues$currentResult <- input$metFragHot_select$select$r + pr[1] - 1
        })

        output$groupHot <- renderRHandsontable({
            gData <- data.frame(group = rownames(gInfo), retention = gInfo$rts, mz = gInfo$mzs, stringsAsFactors = FALSE)
            gData[, unique(anaInfo$group)] <- t(avgGTable)

            hot <- do.call(rhandsontable,
                           c(list(gData, colHeaders = c("Feature group", "Retention", "m/z", unique(anaInfo$group)),
                                  width = NULL, height = 400, maxRows = nrow(gData)),
                             hotOpts))

            return(hot)
        })

        output$metFragHot <- renderRHandsontable({
            cd <- compData()
            hot <- do.call(rhandsontable, c(list(cd, width = NULL, height = 400, maxRows = nrow(cd),
                                                 rowHeaders = seq_len(nrow(cd))), hotOpts)) %>%
                hot_rows(rowHeights = 100) %>%
                hot_col("structure", width = 110, renderer = "
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

        output$spectrum <- renderPlotly({
            spec <- pLists[[getAnalysis()]][[rValues$currentFGroup]][["MSMS"]]
            fi <- fragmentInfo()

            p <- plot_ly(type="scatter", mode = "lines", hoverinfo = "text") %>%
                config(displaylogo = FALSE, scrollZoom = TRUE, collaborate = FALSE,
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
                        htext <- paste0(htext, sprintf("<br>merged by: %s", paste0(unlist(fi$mergedBy[[infoi]]), collapse = ",")))
                }

                p <- add_trace(p, x = rep(spec$mz[i], 2), y = c(0, spec$intensity[i]),
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
            mol <- rcdk::parse.smiles(compTable[[rValues$currentFGroup]]$SMILES[rValues$currentResult])
            if (!is.na(mol))
                rcdkplot(mol[[1]], width = 200, height = 200)
        })

        output$compInfo <- renderText({
            return(getCompInfoText(compTable[[rValues$currentFGroup]], rValues$currentResult, TRUE, input$normalizeScores,
                                   mergedCompoundNames(compounds)))
        })
    }

    runApp(shinyApp(getCompViewerUI(generatePageLabels(names(fGroups)[1], 50)), server))
})
# nocov end
