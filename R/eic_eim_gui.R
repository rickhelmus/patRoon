#' Launch EIC/EIM GUI
#'
#' Launches a Shiny app for generating and plotting EICs and EIMs.
#'
#' @param analysisInfo The analysis information object.
#' @return A Shiny app object.
#' @export
launchEICGUI <- function(analysisInfo)
{
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo)
    
    ui <- shiny::fluidPage(
        shiny::titlePanel("EIC/EIM Viewer"),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::selectInput("analysis", "Analysis", choices = analysisInfo$analysis, selected = analysisInfo$analysis[1]),
                shiny::radioButtons("mzMode", "MZ Specification", choices = c("Range", "Center +/- Window", "Neutral Mass +/- MZ Window", "Formula +/- MZ Window"), selected = "Range", inline = TRUE),
                shiny::conditionalPanel(
                    condition = "input.mzMode == 'Range'",
                    shiny::fluidRow(
                        shiny::column(6, shiny::numericInput("mzmin", "MZ Min", value = 100, min = 0)),
                        shiny::column(6, shiny::numericInput("mzmax", "MZ Max", value = 200, min = 0))
                    )
                ),
                shiny::conditionalPanel(
                    condition = "input.mzMode == 'Center +/- Window'",
                    shiny::fluidRow(
                        shiny::column(6, shiny::numericInput("mzcenter", "MZ Center", value = 150, min = 0)),
                        shiny::column(6, shiny::numericInput("mzwindow", "MZ Window", value = 0.005, min = 0))
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
                ),
                shiny::checkboxInput("useRet", "Use Retention Time Range", value = FALSE),
                shiny::conditionalPanel(
                    condition = "input.useRet == true",
                    shiny::fluidRow(
                        shiny::column(6, shiny::numericInput("retmin", "RT Min", value = 0, min = 0)),
                        shiny::column(6, shiny::numericInput("retmax", "RT Max", value = 0, min = 0))
                    )
                ),
                shiny::checkboxInput("useMob", "Use Mobility Range", value = FALSE),
                shiny::conditionalPanel(
                    condition = "input.useMob == true",
                    shiny::fluidRow(
                        shiny::column(6, shiny::numericInput("mobmin", "Mobility Min", value = 0, min = 0, step = 0.01)),
                        shiny::column(6, shiny::numericInput("mobmax", "Mobility Max", value = 1.3, min = 0, step = 0.01))
                    )
                ),
                shiny::checkboxInput("plotEIM", "Also plot EIM", value = FALSE),
                shiny::conditionalPanel(
                    condition = "input.plotEIM == true",
                    shiny::checkboxInput("enableEIMSmoothing", "Enable EIM Smoothing", value = FALSE),
                    shiny::conditionalPanel(
                        condition = "input.enableEIMSmoothing == true",
                        shiny::numericInput("sgLength", "Smoothing Length", value = 5, min = 3, step = 2)
                    )
                ),
                shiny::actionButton("generate", "Generate plot(s)")
            ),
            shiny::mainPanel(
                shiny::plotOutput("eicPlot"),
                shiny::conditionalPanel(
                    condition = "input.plotEIM == true",
                    shiny::plotOutput("eimPlot")
                )
            )
        )
    )
    
    server <- function(input, output, session)
    {
        eicData <- shiny::eventReactive(input$generate,
        {
            mzmin <- mzmax <- NULL
            if (input$mzMode == "Range")
            {
                mzmin <- input$mzmin
                mzmax <- input$mzmax
            }
            else if (input$mzMode == "Center +/- Window")
            {
                mzmin <- input$mzcenter - input$mzwindow
                mzmax <- input$mzcenter + input$mzwindow
            }
            else if (input$mzMode == "Neutral Mass +/- MZ Window")
            {
                neutralMass <- input$neutralMass
                mz <- tryCatch(calculateMasses(neutralMass, as.adduct(input$adductNM), "mz"), error = function(e)
                {
                    shiny::showNotification("Error calculating MZ from neutral mass and adduct", type = "error")
                    NA
                })
                if (!is.na(mz))
                {
                    mzmin <- mz - input$mzWindowNM
                    mzmax <- mz + input$mzWindowNM
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
                    mz <- tryCatch(calculateMasses(neutralMass, as.adduct(input$adductF), "mz"), error = function(e)
                    {
                        shiny::showNotification("Error calculating MZ from formula and adduct", type = "error")
                        NA
                    })
                    if (!is.na(mz))
                    {
                        mzmin <- mz - input$mzWindowF
                        mzmax <- mz + input$mzWindowF
                    }
                }
            }
            if (is.null(mzmin) || is.null(mzmax))
            {
                mzmin <- 100
                mzmax <- 200
            }
            ranges <- data.table::data.table(
                mzmin = mzmin,
                mzmax = mzmax,
                retmin = if (input$useRet) input$retmin else 0,
                retmax = if (input$useRet) input$retmax else 0
            )
            if (input$useMob)
            {
                ranges$mobmin <- input$mobmin
                ranges$mobmax <- input$mobmax
            }
            rngList <- list(ranges)
            names(rngList) <- input$analysis
            getEICs(analysisInfo[analysis == input$analysis], rngList, gapFactor = 3, output = "pad")
        })
        
        eimData <- shiny::eventReactive(input$generate,
        {
            if (!input$plotEIM)
            {
                return(NULL)
            }
            mzmin <- mzmax <- NULL
            if (input$mzMode == "Range")
            {
                mzmin <- input$mzmin
                mzmax <- input$mzmax
            }
            else if (input$mzMode == "Center +/- Window")
            {
                mzmin <- input$mzcenter - input$mzwindow
                mzmax <- input$mzcenter + input$mzwindow
            }
            else if (input$mzMode == "Neutral Mass +/- MZ Window")
            {
                neutralMass <- input$neutralMass
                mz <- tryCatch(calculateMasses(neutralMass, as.adduct(input$adductNM), "mz"), error = function(e)
                {
                    shiny::showNotification("Error calculating MZ from neutral mass and adduct", type = "error")
                    NA
                })
                if (!is.na(mz))
                {
                    mzmin <- mz - input$mzWindowNM
                    mzmax <- mz + input$mzWindowNM
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
                    mz <- tryCatch(calculateMasses(neutralMass, as.adduct(input$adductF), "mz"), error = function(e)
                    {
                        shiny::showNotification("Error calculating MZ from formula and adduct", type = "error")
                        NA
                    })
                    if (!is.na(mz))
                    {
                        mzmin <- mz - input$mzWindowF
                        mzmax <- mz + input$mzWindowF
                    }
                }
            }
            if (is.null(mzmin) || is.null(mzmax))
            {
                mzmin <- 100
                mzmax <- 200
            }
            eimInfo <- data.table::data.table(
                mzmin = mzmin,
                mzmax = mzmax,
                retmin = if (input$useRet) input$retmin else 0,
                retmax = if (input$useRet) input$retmax else 0,
                mobmin = if (input$useMob) input$mobmin else 0,
                mobmax = if (input$useMob) input$mobmax else 0
            )
            eimList <- list(eimInfo)
            names(eimList) <- input$analysis
            sgLength <- if (input$enableEIMSmoothing) input$sgLength else 0
            shiny::validate(shiny::need(!input$enableEIMSmoothing || (sgLength %% 2 == 1 && sgLength >= 3), "Smoothing length must be odd and >= 3"))
            patRoon:::doGetEIMs(analysisInfo, eimList, minIntensity = 25, sgOrder = 3, sgLength = sgLength)
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
        })
    }
    
    shiny::shinyApp(ui, server)
}
