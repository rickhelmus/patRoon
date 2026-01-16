# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' Launch EIC/EIM GUI
#'
#' Launches a Shiny app for generating and plotting EICs and EIMs.
#'
#' @param obj The object to use for analysis information or features/featureGroups.
#' @param suspects For data.frame method, an optional suspect list data.frame.
#' 
#' @templateVar plain TRUE
#' @template adduct-arg
#' 
#' @return A Shiny app object.
#' @export
setMethod("launchEICGUI", "data.frame", function(obj, suspects = NULL, adduct = NULL)
{
    if (!is.null(adduct))
        adduct <- checkAndToAdduct(adduct)
    
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(obj, add = ac)
    assertSuspectList(suspects, needsAdduct = is.null(adduct), skipInvalid = TRUE, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!is.null(suspects))
    {
        # UNDONE: warn that mobility values may be ignored if no adduct is provided
        suspects <- prepareSuspectList(suspects, adduct = adduct, skipInvalid = TRUE, checkDesc = TRUE,
                                       prefCalcChemProps = TRUE, neutralChemProps = FALSE)
        suspects[, mobility_susp := selectFromSuspAdductCol(suspects, "mobility", if (!is.null(adduct)) as.character(adduct))]
        suspects <- expandSuspMobilities(suspects)
    }
    gui <- createEICGUI(obj, analysisInfo, suspects)
    shiny::shinyApp(gui$ui, gui$server)
})

setMethod("launchEICGUI", "features", function(obj)
{
    analysisInfo <- analysisInfo(obj)
    gui <- createEICGUI(obj, analysisInfo)
    shiny::shinyApp(gui$ui, gui$server)
})

setMethod("launchEICGUI", "featureGroups", function(obj)
{
    analysisInfo <- analysisInfo(obj)
    gui <- createEICGUI(obj, analysisInfo)
    shiny::shinyApp(gui$ui, gui$server)
})

setMethod("launchEICGUI", "featureGroupsScreening", function(obj)
{
    analysisInfo <- analysisInfo(obj)
    suspects <- screenInfo(obj)
    gui <- createEICGUI(obj, analysisInfo, suspects)
    shiny::shinyApp(gui$ui, gui$server)
})

createEICGUI <- function(obj, analysisInfo, suspects = NULL)
{
    ui <- shiny::fluidPage(
        shiny::titlePanel("EIC/EIM Viewer"),
        shiny::sidebarLayout(
            shiny::sidebarPanel(width = 3,
                shiny::tabsetPanel(
                    shiny::tabPanel("MZ",
                        shiny::selectInput("analysis", "Analysis", choices = analysisInfo$analysis, selected = analysisInfo$analysis[1]),
                        shiny::radioButtons("mzMode", "MZ Specification", choices = c("Range", "Center +/- Window", "Neutral Mass +/- MZ Window", "Formula +/- MZ Window"), selected = if (is(obj, "features") || is(obj, "featureGroups")) "Center +/- Window" else "Range", inline = TRUE),
                        shiny::conditionalPanel(
                            condition = "input.mzMode == 'Range'",
                            shiny::fluidRow(
                                shiny::column(6, shiny::numericInput("mzMin", "MZ Min", value = 100, min = 0)),
                                shiny::column(6, shiny::numericInput("mzMax", "MZ Max", value = 200, min = 0))
                            )
                        ),
                        shiny::conditionalPanel(
                            condition = "input.mzMode == 'Center +/- Window'",
                            shiny::fluidRow(
                                shiny::column(6, shiny::numericInput("mzCenter", "MZ Center", value = 150, min = 0)),
                                shiny::column(6, shiny::numericInput("mzWindow", "MZ Window", value = 0.005, min = 0))
                            )
                        ),
                        shiny::conditionalPanel(
                            condition = "input.mzMode == 'Neutral Mass +/- MZ Window'",
                            shiny::fluidRow(
                                shiny::column(4, shiny::numericInput("neutralMass", "Neutral Mass", value = 100, min = 0)),
                                shiny::column(4, shiny::selectInput("adductNM", "Adduct", choices = c("[M+H]+", "[M-H]-", "[M+Na]+", "[M+K]+", "[2M+H]+", "[2M-H]-", "[2M+Na]+", "[2M+K]+"))),
                                shiny::column(4, shiny::numericInput("mzWindowNM", "MZ Window", value = 0.005, min = 0))
                            )
                        ),
                        shiny::conditionalPanel(
                            condition = "input.mzMode == 'Formula +/- MZ Window'",
                            shiny::fluidRow(
                                shiny::column(4, shiny::textInput("formula", "Formula", value = "C6H12O6")),
                                shiny::column(4, shiny::selectInput("adductF", "Adduct", choices = c("[M+H]+", "[M-H]-", "[M+Na]+", "[M+K]+", "[2M+H]+", "[2M-H]-", "[2M+Na]+", "[2M+K]+"))),
                                shiny::column(4, shiny::numericInput("mzWindowF", "MZ Window", value = 0.005, min = 0))
                            )
                        )
                    ),
                    shiny::tabPanel("RT",
                        shiny::radioButtons("rtMode", "RT Specification", choices = c("Disabled", "Range", "Value +/- Window"), selected = if (is(obj, "features") || is(obj, "featureGroups")) "Value +/- Window" else "Disabled", inline = TRUE),
                        shiny::conditionalPanel(
                            condition = "input.rtMode == 'Range'",
                            shiny::fluidRow(
                                shiny::column(6, shiny::numericInput("retMin", "RT Min", value = 0, min = 0)),
                                shiny::column(6, shiny::numericInput("retMax", "RT Max", value = 1000, min = 0))
                            )
                        ),
                        shiny::conditionalPanel(
                            condition = "input.rtMode == 'Value +/- Window'",
                            shiny::fluidRow(
                                shiny::column(6, shiny::numericInput("retValue", "RT Value", value = 5, min = 0)),
                                shiny::column(6, shiny::numericInput("retWindow", "RT Window", value = 30, min = 0))
                            )
                        )
                    ),
                    shiny::tabPanel("Mobility",
                        shiny::radioButtons("mobMode", "Mobility Specification", choices = c("Disabled", "Range", "Value +/- Window"), selected = if (is(obj, "features") || is(obj, "featureGroups")) "Value +/- Window" else "Disabled", inline = TRUE),
                        shiny::conditionalPanel(
                            condition = "input.mobMode == 'Range'",
                            shiny::fluidRow(
                                shiny::column(6, shiny::numericInput("mobMin", "Mobility Min", value = 0, min = 0, step = 0.01)),
                                shiny::column(6, shiny::numericInput("mobMax", "Mobility Max", value = 1.3, min = 0, step = 0.01))
                            )
                        ),
                        shiny::conditionalPanel(
                            condition = "input.mobMode == 'Value +/- Window'",
                            shiny::fluidRow(
                                shiny::column(6, shiny::numericInput("mobValue", "Mobility Value", value = 0.5, min = 0, step = 0.01)),
                                shiny::column(6, shiny::numericInput("mobWindow", "Mobility Window", value = 0.1, min = 0, step = 0.01))
                            )
                        )
                    ),
                    shiny::tabPanel("EIM & Options",
                        shiny::checkboxInput("plotEIM", "Also plot EIM", value = FALSE),
                        shiny::conditionalPanel(
                            condition = "input.plotEIM == true",
                            shiny::selectInput("eimSmooth", "EIM Smoothing Method", choices = c("none", "ma", "sg"), selected = "none"),
                            shiny::conditionalPanel(
                                condition = "input.eimSmooth != 'none'",
                                shiny::numericInput("sgLength", "Smoothing Length", value = 5, min = 3, step = 2)
                            ),
                            shiny::conditionalPanel(
                                condition = "input.rtMode == 'Range'",
                                shiny::fluidRow(
                                    shiny::column(6, shiny::numericInput("eimRetMin", "EIM RT Min", value = 0, min = 0)),
                                    shiny::column(6, shiny::numericInput("eimRetMax", "EIM RT Max", value = 10, min = 0))
                                )
                            ),
                            shiny::conditionalPanel(
                                condition = "input.rtMode == 'Value +/- Window'",
                                shiny::fluidRow(
                                    shiny::column(6, shiny::numericInput("eimRetValue", "EIM RT Value", value = 5, min = 0)),
                                    shiny::column(6, shiny::numericInput("eimRetWindow", "EIM RT Window", value = 15, min = 0))
                                )
                            )
                        ),
                        shiny::checkboxInput("plotChromPoints", "Also plot Chromatogram Points", value = FALSE),
                        shiny::checkboxInput("autoGenerate", "Auto-generate plots on input change", value = TRUE),
                        shiny::conditionalPanel(
                            condition = "input.autoGenerate == false",
                            shiny::actionButton("generate", "Generate plot(s)")
                        )
                    ),
                    shiny::tabPanel("Peaks",
                        shiny::radioButtons("peakTarget", "Detect peaks in", choices = c("EIC", "EIM"), selected = "EIC", inline = TRUE),
                        shiny::selectInput("peakAlgorithm", "Algorithm", choices = c("openms", "xcms3", "envipick", "piek"), selected = "piek"),
                        shiny::selectInput("peakIMSType", "IMS Type", choices = c("bruker_ims", "agilent_ims"), selected = "bruker_ims"),
                        shiny::actionButton("editPeakParams", "Advanced peak params"),
                        shiny::actionButton("resetPeakParams", "Restore default peak params"),
                        shiny::verbatimTextOutput("peakParamsSummary"),
                        shiny::conditionalPanel(
                            condition = "input.autoGenerate == false",
                            shiny::actionButton("detectPeaks", "Detect peaks")
                        )
                    )
                )
            ),
            shiny::mainPanel(
                shiny::uiOutput("plotsUI"),
                shiny::uiOutput("tablesUI")
            )
        )
    )
    
    server <- function(input, output, session)
    {
        computeMzRange <- function(input)
        {
            mzMin <- mzMax <- NULL
            if (input$mzMode == "Range")
            {
                mzMin <- input$mzMin
                mzMax <- input$mzMax
            }
            else if (input$mzMode == "Center +/- Window")
            {
                mzMin <- input$mzCenter - input$mzWindow
                mzMax <- input$mzCenter + input$mzWindow
            }
            else if (input$mzMode == "Neutral Mass +/- MZ Window")
            {
                neutralMass <- input$neutralMass
                mzVal <- tryCatch(calculateMasses(neutralMass, as.adduct(input$adductNM), "mz"), error = function(e)
                {
                    shiny::showNotification("Error calculating MZ from neutral mass and adduct", type = "error")
                    NA
                })
                if (!is.na(mzVal))
                {
                    mzMin <- mzVal - input$mzWindowNM
                    mzMax <- mzVal + input$mzWindowNM
                }
            }
            else if (input$mzMode == "Formula +/- MZ Window")
            {
                neutralMass <- tryCatch(getFormulaMass(input$formula), error = function(e)
                {
                    shiny::showNotification("Error parsing formula", type = "error")
                    NA
                })
                if (!is.na(neutralMass))
                {
                    mzVal <- tryCatch(calculateMasses(neutralMass, as.adduct(input$adductF), "mz"), error = function(e)
                    {
                        shiny::showNotification("Error calculating MZ from formula and adduct", type = "error")
                        NA
                    })
                    if (!is.na(mzVal))
                    {
                        mzMin <- mzVal - input$mzWindowF
                        mzMax <- mzVal + input$mzWindowF
                    }
                }
            }
            if (is.null(mzMin) || is.null(mzMax))
            {
                mzMin <- 100
                mzMax <- 200
            }
            list(mzmin = mzMin, mzmax = mzMax)
        }

        plotPeakOverlay <- function(peaks, mat, xName)
        {
            if (!is.null(peaks) && length(peaks) > 0)
            {
                for (r in seq_len(nrow(peaks)))
                {
                    prow <- peaks[r, ]
                    sel <- mat[, xName] >= prow$retmin & mat[, xName] <= prow$retmax
                    EIXFill <- mat[sel, , drop = FALSE]
                    if (nrow(EIXFill) > 0)
                    {
                        col <- adjustcolor("blue", alpha.f = 0.35)
                        polygon(c(EIXFill[, xName], rev(EIXFill[, xName])),
                                c(EIXFill[, "intensity"], rep(0, nrow(EIXFill))),
                                col = col, border = NA)
                    }
                }
            }
        }

        updateInputsFromFeat <- function(feat)
        {
            if (!is.na(feat$mz))
            {
                updateNumericInput(session, "mzCenter", value = round(feat$mz, 5))
                updateNumericInput(session, "mzWindow", value = round((feat$mzmax - feat$mzmin) / 2, 5))
                updateNumericInput(session, "mzMin", value = round(feat$mzmin, 5))
                updateNumericInput(session, "mzMax", value = round(feat$mzmax, 5))
            }
            if (!is.na(feat$ret))
            {
                updateNumericInput(session, "retValue", value = round(feat$ret, 2))
                updateNumericInput(session, "retWindow", value = max(15, round((feat$retmax - feat$retmin) / 2), 2))
                updateNumericInput(session, "retMin", value = round(feat$retmin, 2))
                updateNumericInput(session, "retMax", value = round(feat$retmax, 2))
                updateNumericInput(session, "eimRetValue", value = round(feat$ret, 2))
                updateNumericInput(session, "eimRetWindow", value = round((feat$retmax - feat$retmin) / 2, 2))
                updateNumericInput(session, "eimRetMin", value = round(feat$retmin, 2))
                updateNumericInput(session, "eimRetMax", value = round(feat$retmax, 2))
            }
            if (!is.na(feat$mobility) && !is.na(feat$mobmin))
            {
                updateNumericInput(session, "mobValue", value = round(feat$mobility, 3))
                updateNumericInput(session, "mobWindow", value = round((feat$mobmax - feat$mobmin) / 2, 3))
                updateNumericInput(session, "mobMin", value = round(feat$mobmin, 3))
                updateNumericInput(session, "mobMax", value = round(feat$mobmax, 3))
            }
        }

        updateInputsFromSuspect <- function(susp)
        {
            if (!is.na(susp$mz))
            {
                updateNumericInput(session, "mzCenter", value = round(susp$mz, 5))
                updateNumericInput(session, "mzWindow", value = 0.005)
                updateNumericInput(session, "mzMin", value = round(susp$mz - 0.005, 5))
                updateNumericInput(session, "mzMax", value = round(susp$mz + 0.005, 5))
            }
            if (!is.null(susp[["rt"]]) && !is.na(susp$rt))
            {
                updateNumericInput(session, "retValue", value = round(susp$rt, 2))
                updateNumericInput(session, "retWindow", value = 15)
                updateNumericInput(session, "retMin", value = round(susp$rt - 15, 2))
                updateNumericInput(session, "retMax", value = round(susp$rt + 15, 2))
                updateNumericInput(session, "eimRetValue", value = round(susp$rt, 2))
                updateNumericInput(session, "eimRetWindow", value = 15)
                updateNumericInput(session, "eimRetMin", value = round(susp$rt - 15, 2))
                updateNumericInput(session, "eimRetMax", value = round(susp$rt + 15, 2))
            }
            if (!is.null(susp[["mobility"]]) && !is.na(susp$mobility))
            {
                updateNumericInput(session, "mobValue", value = round(susp$mobility, 3))
                updateNumericInput(session, "mobWindow", value = 0.1)
                updateNumericInput(session, "mobMin", value = round(susp$mobility - 0.1, 3))
                updateNumericInput(session, "mobMax", value = round(susp$mobility + 0.1, 3))
            }
        }

        featureData <- shiny::reactive(
        {
            if (is(obj, "features"))
            {
                ft <- featureTable(obj)[[input$analysis]]
                data.frame(
                    ID = ft$ID,
                    mz = round(ft$mz, 5),
                    rt = round(ft$ret, 2),
                    mobility = if ("mobility" %in% names(ft)) round(ft$mobility, 3) else NA,
                    mzmin = round(ft$mzmin, 5),
                    mzmax = round(ft$mzmax, 5),
                    retmin = round(ft$retmin, 2),
                    retmax = round(ft$retmax, 2),
                    mobmin = if ("mobmin" %in% names(ft)) round(ft$mobmin, 3) else NA,
                    mobmax = if ("mobmax" %in% names(ft)) round(ft$mobmax, 3) else NA
                )
            }
            else if (is(obj, "featureGroups"))
            {
                ft <- featureTable(obj)[[input$analysis]]
                groups <- unique(ft$group)
                data.frame(
                    Group = groups,
                    mz = sapply(groups, function(g) round(ft[group == g][1]$mz, 5)),
                    rt = sapply(groups, function(g) round(ft[group == g][1]$ret, 2)),
                    mobility = sapply(groups, function(g) {
                        mob <- ft[group == g][1]$mobility
                        if (!is.na(mob)) round(mob, 3) else NA
                    }),
                    mzmin = sapply(groups, function(g) round(ft[group == g][1]$mzmin, 5)),
                    mzmax = sapply(groups, function(g) round(ft[group == g][1]$mzmax, 5)),
                    retmin = sapply(groups, function(g) round(ft[group == g][1]$retmin, 2)),
                    retmax = sapply(groups, function(g) round(ft[group == g][1]$retmax, 2)),
                    mobmin = sapply(groups, function(g) {
                        mobmin <- ft[group == g][1]$mobmin
                        if (!is.na(mobmin)) round(mobmin, 3) else NA
                    }),
                    mobmax = sapply(groups, function(g) {
                        mobmax <- ft[group == g][1]$mobmax
                        if (!is.na(mobmax)) round(mobmax, 3) else NA
                    })
                )
            }
            else
                NULL
        })

        output$featureTable <- DT::renderDT({
            req(featureData())
            DT::datatable(featureData(), selection = 'single', options = list(pageLength = 10))
        })

        observeEvent(input$featureTable_rows_selected, {
            req(featureData())
            selected <- input$featureTable_rows_selected
            if (length(selected) == 0) return()
            feat <- featureData()[selected, ]
            if (is(obj, "features"))
            {
                ft <- featureTable(obj)[[input$analysis]]
                feat_full <- ft[ft$ID == feat$ID, ]
            }
            else if (is(obj, "featureGroups"))
            {
                ft <- featureTable(obj)[[input$analysis]]
                feat_full <- ft[group == feat$Group][1]
            }
            if (nrow(feat_full) == 0) return()
            updateInputsFromFeat(feat_full)
        })
        # Sync EIM RT inputs with EIC RT inputs
        observeEvent(c(input$retMin, input$retMax), {
            updateNumericInput(session, "eimRetMin", value = input$retMin)
            updateNumericInput(session, "eimRetMax", value = input$retMax)
        })
        
        observeEvent(c(input$retValue, input$retWindow), {
            updateNumericInput(session, "eimRetValue", value = input$retValue)
            updateNumericInput(session, "eimRetWindow", value = input$retWindow)
        })
        
        suspectData <- shiny::reactive({
            if (!is.null(suspects))
            {
                data.frame(
                    Name = suspects$name,
                    mz = if ("mz" %in% names(suspects)) round(suspects$mz, 5) else NA,
                    rt = if ("rt" %in% names(suspects)) round(suspects$rt, 2) else NA,
                    mobility = if ("mobility" %in% names(suspects)) round(suspects$mobility, 3) else NA,
                    adduct = if ("adduct" %in% names(suspects)) suspects$adduct else NA,
                    formula = if ("formula" %in% names(suspects)) suspects$formula else NA,
                    group = if ("group" %in% names(suspects)) suspects$group else NA
                )
            }
            else
                NULL
        })

        output$suspectTable <- DT::renderDT({
            req(suspectData())
            DT::datatable(suspectData(), selection = 'single', options = list(pageLength = 10))
        })

        eicTableData <- shiny::reactive({
            req(eicData())
            eic <- eicData()[[1]][[1]]  # Assuming single analysis and single range
            if (is.null(eic) || nrow(eic) == 0) return(NULL)
            data.frame(
                Time = eic[, "time"],
                Intensity = eic[, "intensity"]
            )
        })

        output$eicTable <- DT::renderDT({
            req(eicTableData())
            DT::datatable(eicTableData(), options = list(pageLength = 10))
        })

        eimTableData <- shiny::reactive({
            req(input$plotEIM, eimData())
            eim <- eimData()[[1]][[1]]  # Assuming single analysis and single range
            if (is.null(eim) || nrow(eim) == 0) return(NULL)
            data.frame(
                Mobility = eim[, "mobility"],
                Intensity = eim[, "intensity"]
            )
        })

        output$eimTable <- DT::renderDT({
            req(eimTableData())
            DT::datatable(eimTableData(), options = list(pageLength = 10))
        })

        chromPointsTableData <- shiny::reactive({
            req(input$plotChromPoints, chromPointsData())
            chromData <- chromPointsData()
            if (is.null(chromData) || nrow(chromData) == 0) return(NULL)
            data.frame(
                MZ = chromData$mz,
                Intensity = chromData$intensity
            )
        })

        output$chromPointsTable <- DT::renderDT({
            req(chromPointsTableData())
            DT::datatable(chromPointsTableData(), options = list(pageLength = 10))
        })

        observeEvent(input$suspectTable_rows_selected, {
            req(suspectData())
            selected <- input$suspectTable_rows_selected
            if (length(selected) == 0) return()
            susp <- suspects[selected, ]
            if (nrow(susp) == 0) return()
            # For featureGroupsScreening, use the corresponding feature
            if (is(obj, "featureGroupsScreening") && "group" %in% names(suspects) && !is.na(susp$group))
            {
                ft <- featureTable(obj)[[input$analysis]]
                feat <- ft[group == susp$group][1]
                if (nrow(feat) > 0)
                {
                    updateInputsFromFeat(feat)
                }
            }
            else
            {
                updateInputsFromSuspect(susp)
            }
            # Update Adduct and Formula from suspect
            if (!is.null(susp[["adduct"]]) && !is.na(susp$adduct))
            {
                updateSelectInput(session, "adductNM", selected = susp$adduct)
                updateSelectInput(session, "adductF", selected = susp$adduct)
            }
            if (!is.na(susp$formula))
            {
                updateTextInput(session, "formula", value = susp$formula)
            }
        })

        output$plotsUI <- shiny::renderUI({
            if (input$plotEIM && input$plotChromPoints) {
                shiny::fluidRow(
                    shiny::column(4, shiny::plotOutput("eicPlot")),
                    shiny::column(4, shiny::plotOutput("eimPlot")),
                    shiny::column(4, shiny::plotOutput("chromPointsPlot"))
                )
            } else if (input$plotEIM) {
                shiny::fluidRow(
                    shiny::column(6, shiny::plotOutput("eicPlot")),
                    shiny::column(6, shiny::plotOutput("eimPlot"))
                )
            } else if (input$plotChromPoints) {
                shiny::fluidRow(
                    shiny::column(6, shiny::plotOutput("eicPlot")),
                    shiny::column(6, shiny::plotOutput("chromPointsPlot"))
                )
            } else {
                shiny::plotOutput("eicPlot")
            }
        })

        output$tablesUI <- shiny::renderUI({
            shiny::tabsetPanel(
                if (is(obj, "features") || is(obj, "featureGroups"))
                    shiny::tabPanel("Features", DT::dataTableOutput("featureTable")),
                if (!is.null(suspects))
                    shiny::tabPanel("Suspects", DT::dataTableOutput("suspectTable")),
                shiny::tabPanel("EIC", DT::dataTableOutput("eicTable")),
                if (input$plotEIM)
                    shiny::tabPanel("EIM", DT::dataTableOutput("eimTable")),
                if (input$plotChromPoints)
                    shiny::tabPanel("Chromatogram Points", DT::dataTableOutput("chromPointsTable")),
                shiny::tabPanel("Peaks", DT::dataTableOutput("peaksTable"))
            )
        })

        eicData <- shiny::reactive({
            req(input$mzMode, input$rtMode, input$mobMode, input$analysis)
            if (!input$autoGenerate)
                req(input$generate)

            mzRange <- computeMzRange(input)
            mzmin <- mzRange$mzmin
            mzmax <- mzRange$mzmax

            ranges <- data.table::data.table(
                mzmin = mzmin,
                mzmax = mzmax,
                retmin = if (input$rtMode == "Disabled") 0 else if (input$rtMode == "Range") input$retMin else input$retValue - input$retWindow,
                retmax = if (input$rtMode == "Disabled") 0 else if (input$rtMode == "Range") input$retMax else input$retValue + input$retWindow
            )
            if (input$mobMode != "Disabled")
            {
                ranges$mobmin <- if (input$mobMode == "Range") input$mobMin else input$mobValue - input$mobWindow
                ranges$mobmax <- if (input$mobMode == "Range") input$mobMax else input$mobValue + input$mobWindow
            }
            rngList <- list(ranges)
            names(rngList) <- input$analysis
            print(ranges)
            getEICs(analysisInfo[analysis == input$analysis], rngList, gapFactor = 3, output = "pad")
        })
        
        eimData <- shiny::reactive({
            req(input$plotEIM, input$mzMode, input$rtMode, input$mobMode, input$analysis)
            if (!input$autoGenerate)
                req(input$generate)
            
            if (!input$plotEIM)
            {
                return(NULL)
            }

            mzRange <- computeMzRange(input)
            mzmin <- mzRange$mzmin
            mzmax <- mzRange$mzmax

            eimInfo <- data.table::data.table(
                mzmin = mzmin,
                mzmax = mzmax,
                retmin = if (input$rtMode == "Disabled") 0 else if (input$rtMode == "Range") input$eimRetMin else input$eimRetValue - input$eimRetWindow,
                retmax = if (input$rtMode == "Disabled") 0 else if (input$rtMode == "Range") input$eimRetMax else input$eimRetValue + input$eimRetWindow,
                mobmin = if (input$mobMode == "Disabled") 0 else if (input$mobMode == "Range") input$mobMin else input$mobValue - input$mobWindow,
                mobmax = if (input$mobMode == "Disabled") 0 else if (input$mobMode == "Range") input$mobMax else input$mobValue + input$mobWindow
            )
            eimList <- list(eimInfo)
            names(eimList) <- input$analysis
            shiny::validate(shiny::need(input$eimSmooth == 'none' || (input$sgLength %% 2 == 1 && input$sgLength >= 3), "Smoothing length must be odd and >= 3"))
            patRoon:::doGetEIMs(analysisInfo, eimList, minIntensity = 25, smooth = input$eimSmooth, sgOrder = 3, smLength = input$sgLength)
        })

        chromPointsData <- shiny::reactive({
            req(input$plotChromPoints, input$mzMode, input$rtMode, input$analysis)
            if (!input$autoGenerate)
                req(input$generate)

            if (!input$plotChromPoints)
            {
                return(NULL)
            }

            mzRange <- computeMzRange(input)
            mzmin <- mzRange$mzmin
            mzmax <- mzRange$mzmax

            pointsInfo <- data.table::data.table(
                mzmin = mzmin,
                mzmax = mzmax,
                retmin = if (input$rtMode == "Disabled") 0 else if (input$rtMode == "Range") input$retMin else input$retValue - input$retWindow,
                retmax = if (input$rtMode == "Disabled") 0 else if (input$rtMode == "Range") input$retMax else input$retValue + input$retWindow,
                mobmin = if (input$mobMode == "Disabled") 0 else if (input$mobMode == "Range") input$mobMin else input$mobValue - input$mobWindow,
                mobmax = if (input$mobMode == "Disabled") 0 else if (input$mobMode == "Range") input$mobMax else input$mobValue + input$mobWindow
            )
            
            if (pointsInfo$mobmin == 0 && pointsInfo$mobmax == 0)
                pointsInfo[, c("mobmin", "mobmax") := NULL]
            printf("points: %s\n", paste(capture.output(print(pointsInfo)), collapse = "; "))
            
            pointsList <- list(pointsInfo)
            names(pointsList) <- input$analysis

            chromPoints <- doGetChromPoints(analysisInfo[analysis == input$analysis], pointsList)

            # Process the data: aggregate by MZ, sum intensities
            if (length(chromPoints) > 0 && !is.null(chromPoints[[1]][[1]]) && nrow(chromPoints[[1]][[1]]) > 0)
            {
                points <- as.data.table(chromPoints[[1]][[1]])
                aggData <- points[, .(intensity = sum(intensity)), by = mz][order(mz)]
                aggData <- unique(aggData, by = "mz")
                return(aggData)
            }
            return(NULL)
        })

        ## Peak detection reactive storage (persist per target / algorithm / IMS type)
        peaksRV <- shiny::reactiveVal(NULL)
        peaksDT <- shiny::reactiveVal(data.table::data.table())
        peakParams <- shiny::reactiveVal(NULL)

        # store: nested list stored as reactiveVal; keys: target -> algorithm -> imstype (or "default")
        peakParamsStore <- shiny::reactiveVal(list())

        getStoreParams <- function(target, algorithm, imstype) {
          store <- peakParamsStore()
          if (!is.list(store)) return(NULL)
          if (!is.null(store[[target]]) && !is.null(store[[target]][[algorithm]])) {
            # imstype-aware: prefer exact imstype, then fallback to algorithm-level
            algNode <- store[[target]][[algorithm]]
            if (!is.null(algNode[[imstype]])) return(algNode[[imstype]])
            if (!is.null(algNode[["default"]])) return(algNode[["default"]])
          }
          NULL
        }

        setStoreParams <- function(target, algorithm, imstype, params) {
          store <- peakParamsStore()
          if (!is.list(store)) store <- list()
          if (is.null(store[[target]])) store[[target]] <- list()
          if (is.null(store[[target]][[algorithm]])) store[[target]][[algorithm]] <- list()
          # for EIM we store under imstype, for EIC we store under "default" (or still imstype)
          key <- if (is.null(imstype) || imstype == "") "default" else as.character(imstype)
          store[[target]][[algorithm]][[key]] <- params
          peakParamsStore(store)
        }

        # initialize default peak params when algorithm/target/type changes
        shiny::observeEvent(list(input$peakAlgorithm, input$peakTarget, input$peakIMSType), {
          target <- as.character(input$peakTarget)
          algorithm <- as.character(input$peakAlgorithm)
          imstype <- ifelse(is.null(input$peakIMSType), "default", as.character(input$peakIMSType))
          # try store first
          sp <- tryCatch(getStoreParams(target, algorithm, imstype), error = function(e) NULL)
          if (!is.null(sp)) {
            peakParams(sp)
            return()
          }
          # fallback to default provider
          typeArg <- if (isTRUE(target == "EIC")) "chrom" else imstype
          p <- tryCatch(getDefPeakParams(type = typeArg, algorithm = algorithm), error = function(e) list(algorithm = algorithm))
          peakParams(p)
          # also store as default for this key so user's subsequent edits persist even if they haven't saved
          setStoreParams(target, algorithm, imstype, p)
        }, ignoreNULL = FALSE)
        
        output$peakParamsSummary <- shiny::renderPrint({
            req(peakParams())
            print(peakParams())
        })
        
        # modal for editing peak params (dynamic per-parameter inputs)
        shiny::observeEvent(input$editPeakParams, {
            p <- peakParams()
            if (is.null(p))
            {
                shiny::showNotification("Peak parameters not initialized", type = "error")
                return()
            }
            # build dynamic inputs based on parameter types
            inputs <- lapply(names(p), function(nm) {
                id <- paste0("pp_", nm)
                val <- p[[nm]]
                # logical scalar
                if (is.logical(val) && length(val) == 1)
                    return(shiny::checkboxInput(id, nm, value = val))
                # numeric scalar
                if (is.numeric(val) && length(val) == 1)
                    return(shiny::numericInput(id, nm, value = val))
                # numeric vector length 2 -> two numeric inputs
                if (is.numeric(val) && length(val) == 2)
                    return(shiny::fluidRow(
                        shiny::column(6, shiny::numericInput(paste0(id, "_1"), paste0(nm, " (min)"), value = val[1])),
                        shiny::column(6, shiny::numericInput(paste0(id, "_2"), paste0(nm, " (max)"), value = val[2]))
                    ))
                # character scalar
                if (is.character(val) && length(val) == 1)
                    return(shiny::textInput(id, nm, value = val))
                # fallback: JSON textarea for complex/list values
                shiny::tagList(
                    paste(id, nm),
                    shiny::textAreaInput(id, NULL, value = tryCatch(jsonlite::toJSON(val, auto_unbox = TRUE, pretty = TRUE), error = function(e) ""),
                                         rows = 4, width = "100%")
                )
            })
            shiny::showModal(shiny::modalDialog(
                title = "Edit peak detection parameters",
                shiny::tagList(inputs),
                footer = shiny::tagList(
                    shiny::modalButton("Cancel"),
                    shiny::actionButton("savePeakParams", "Save")
                ),
                size = "l",
                easyClose = TRUE
            ))
        })
        
        shiny::observeEvent(input$savePeakParams, {
            oldp <- peakParams()
            if (is.null(oldp))
                return()
            newp <- list()
            for (nm in names(oldp))
            {
                id <- paste0("pp_", nm)
                # vector inputs
                val1 <- NULL; val2 <- NULL
                if (!is.null(input[[paste0(id, "_1")]]) || !is.null(input[[paste0(id, "_2")]]))
                {
                    val1 <- input[[paste0(id, "_1")]]
                    val2 <- input[[paste0(id, "_2")]]
                    newp[[nm]] <- c(val1, val2)
                    next
                }
                # simple scalar inputs
                if (!is.null(input[[id]]))
                {
                    v <- input[[id]]
                    # if textarea was used for complex values, attempt JSON parse
                    if (is.character(v) && grepl("^\\s*\\{", v) || grepl("^\\s*\\[", v))
                    {
                        parsed <- tryCatch(jsonlite::fromJSON(v, simplifyVector = FALSE), error = function(e) v)
                        newp[[nm]] <- parsed
                    }
                    else
                        newp[[nm]] <- v
                    next
                }
                # fallback to old value
                newp[[nm]] <- oldp[[nm]]
            }
            # ensure algorithm key present
            if (!is.null(oldp$algorithm) && is.null(newp$algorithm))
                newp$algorithm <- oldp$algorithm
            # persist to store
            target <- as.character(input$peakTarget)
            algorithm <- as.character(ifelse(is.null(input$peakAlgorithm), "unknown", input$peakAlgorithm))
            imstype <- ifelse(is.null(input$peakIMSType), "default", as.character(input$peakIMSType))
            setStoreParams(target, algorithm, imstype, newp)
            peakParams(newp)
            shiny::removeModal()
        })

        # restore defaults observer
        shiny::observeEvent(input$resetPeakParams, {
            target <- as.character(input$peakTarget)
            algorithm <- as.character(input$peakAlgorithm)
            imstype <- ifelse(is.null(input$peakIMSType), "default", as.character(input$peakIMSType))
            # get defaults
            typeArg <- if (isTRUE(target == "EIC")) "chrom" else imstype
            defaults <- tryCatch(getDefPeakParams(type = typeArg, algorithm = algorithm), error = function(e) list(algorithm = algorithm))
            # overwrite stored settings
            setStoreParams(target, algorithm, imstype, defaults)
            # update current
            peakParams(defaults)
            shiny::showNotification("Peak parameters restored to defaults", type = "message")
        })

        # re-run peak detection whenever peak parameters are updated (only when auto-generate is enabled)
        shiny::observeEvent(list(peakParams(), input$autoGenerate), {
          req(peakParams())
          if (isTRUE(input$autoGenerate)) {
            runDetect()
          }
        }, ignoreNULL = TRUE)

        # Run detection function
        runDetect <- function()
        {
            p <- peakParams()
            if (is.null(p))
                return(NULL)
            
            # prepare EIC-like list for findPeaks
            EICs <- NULL
            fillEICs <- NULL
            if (isTRUE(input$peakTarget == "EIC"))
            {
                ed <- eicData()
                if (is.null(ed) || length(ed) == 0)
                {
                    shiny::showNotification("No EIC data available for peak detection", type = "warning")
                    return(NULL)
                }
                # assume single analysis
                EICs <- ed[[1]]
                fillEICs <- TRUE
            }
            else # EIM
            {
                ed <- eimData()
                if (is.null(ed) || length(ed) == 0)
                {
                    shiny::showNotification("No EIM data available for peak detection", type = "warning")
                    return(NULL)
                }
                eims <- ed[[1]]
                # convert mobility -> time for peak detection
                EICs <- lapply(eims, function(mat) {
                    if (is.null(mat) || nrow(mat) == 0) return(mat)
                    m <- as.matrix(mat)
                    cn <- colnames(m)
                    if ("mobility" %in% cn)
                        cn[cn == "mobility"] <- "time"
                    colnames(m) <- cn
                    return(m)
                })
                fillEICs <- FALSE
            }
            
            # ensure named list
            if (is.null(names(EICs)))
                names(EICs) <- as.character(seq_along(EICs))
            
            logPath <- tempfile(fileext = ".log")
            peaks <- tryCatch({
                findPeaks(EICs, fillEICs, p, logPath)
            }, error = function(e) {
                shiny::showNotification(sprintf("Peak finding failed: %s", e$message), type = "error")
                return(NULL)
            })
            if (is.null(peaks))
                return(NULL)
            
            if (length(peaks) == 0 || nrow(peaks[[1]]) == 0)
            {
                peaksRV(data.table::data.table())
                peaksDT(data.table::data.table())
            }
            else
            {
                peaks <- peaks[[1]]
                peaksRV(peaks)
                peaksdt <- copy(peaks)
                peaksdt[, analysis := input$analysis]
                peaksdt[, algorithm := input$peakAlgorithm]
                peaksDT(peaksdt)
            }
            invisible(NULL)
        }
        
        # detect on demand or when autoGenerate is TRUE and inputs/data change
        shiny::observeEvent(input$detectPeaks, {
            runDetect()
        })
        shiny::observeEvent(list(eicData(), input$peakAlgorithm, input$peakTarget, input$peakIMSType), {
            if (input$autoGenerate)
                runDetect()
        }, ignoreNULL = FALSE)
        shiny::observeEvent(eicData(), {
            if (input$autoGenerate && input$peakTarget == "EIM")
                runDetect()
        }, ignoreNULL = FALSE)
        
        output$peaksTable <- DT::renderDT({
            req(peaksDT())
            DT::datatable(peaksDT(), selection = 'single', options = list(pageLength = 10))
        })
        
        output$eicPlot <- shiny::renderPlot(
        {
            req(eicData())
            eic <- eicData()[[1]][[1]]  # Assuming single analysis and single range
            if (is.null(eic) || nrow(eic) == 0)
            {
                plot.new()
                text(0.5, 0.5, "No EIC data", cex = 2)
                return()
            }
            plot(eic[, "time"], eic[, "intensity"], type = "l", xlab = "Retention Time", ylab = "Intensity",
                 main = "Extracted Ion Chromatogram")
            
            peaks <- peaksRV()
            plotPeakOverlay(peaks, eic, "time")
        })
        
        output$eimPlot <- shiny::renderPlot(
        {
            req(eimData())
            eim <- eimData()[[1]][[1]]  # Assuming single analysis and single range
            if (is.null(eim) || nrow(eim) == 0)
            {
                plot.new()
                text(0.5, 0.5, "No EIM data", cex = 2)
                return()
            }
            plot(eim[, "mobility"], eim[, "intensity"], type = "l", xlab = "Mobility", ylab = "Intensity",
                 main = "Extracted Ion Mobilogram")

            peaks <- peaksRV()
            plotPeakOverlay(peaks, eim, "mobility")
        })

        output$chromPointsPlot <- shiny::renderPlot(
        {
            req(chromPointsData())
            chromData <- chromPointsData()
            if (is.null(chromData) || nrow(chromData) == 0)
            {
                plot.new()
                text(0.5, 0.5, "No Chromatogram Points data", cex = 2)
                return()
            }
            plot(chromData$mz, chromData$intensity, type = "l", xlab = "m/z", ylab = "Intensity",
                 main = "Chromatogram Points")
        })
    }
    
    return(list(ui = ui, server = server))
}
