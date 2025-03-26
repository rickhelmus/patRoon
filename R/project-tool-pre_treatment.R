splitConvTypeFormat <- function(typeFormat) strsplit(typeFormat, "_")

newProjectPreTreatUI <- function(id)
{
    ns <- NS(id)
    htmltools::tagList(
        miniUI::miniContentPanel(
            fillCol(
                flex = c(NA, 1),
                fillRow(
                    height = 30,
                    strong("MS data conversion steps")
                ),
                fillRow(
                    rhandsontable::rHandsontableOutput(ns("MSConversionHOT"))
                )
            )
        ),
        miniUI::miniButtonBlock(
            actionButton(ns("addMSConversion"), label = "Add", icon = icon("plus"),
                         title = "Add MS conversion step"),
            shinyjs::disabled(actionButton(ns("removeMSConversionRows"), label = "Remove", icon = icon("trash"),
                                           title = "Remove selected row(s)")),
            shinyjs::disabled(actionButton(ns("MSConversionRowsUp"), label = "Move up", icon = icon("arrow-up"),
                                           title = "Move selected row(s) up")),
            shinyjs::disabled(actionButton(ns("MSConversionRowsDown"), label = "Move down", icon = icon("arrow-down"),
                                           title = "Move selected row(s) down")),
            actionButton(ns("configConvPaths"), label = "Output paths", icon = icon("folder-open")),
            actionButton(ns("configBrukerCalib"), label = "Bruker m/z calibration", icon = icon("cog"),
                         title = "Configure Bruker DataAnalysis m/z calibration")
        )
    )
}

newProjectPreTreatServer <- function(id, ionization, settings)
{
    ns <- NS(id)
    
    groupConvTypesFormats <- function(algo, inOut)
    {
        fTypes <- intersect(getMSFileTypes(), if (inOut == "input") validConvFromTypes(algo) else validConvToTypes(algo))
        return(sapply(fTypes, function(ft)
        {
            formats <- if (inOut == "input") getMSInConversionFormats(algo, ft) else getMSOutConversionFormats(algo, ft)
            return(setNames(paste0(ft, "_", formats), paste0(formats, " (", ft, ")")))
        }, simplify = FALSE))
    }
    
    convTypeFormatToLabel <- function(typeFormat)
    {
        tf <- splitConvTypeFormat(typeFormat)
        return(paste0(tf[[1]][1], " (", tf[[1]][2], ")"))
    }
    
    moduleServer(id, function(input, output, session)
    {
        rValues <- reactiveValues(
            triggerMSConvHOTUpdate = 0,
            steps = data.table(algorithm = character(0), from = character(0), to = character(0)),
            output = list(centroid = file.path("converted", "centroid"), profile = file.path("converted", "profile"),
                          ims = file.path("converted", "ims")),
            brukerCalib = list(enabled = FALSE, method = "", methodPos = "", methodNeg = "")
        )
        
        triggerMSConvHOTUpdate <- function() rValues$triggerMSConvHOTUpdate <- rValues$triggerMSConvHOTUpdate + 1
        
        observeEvent(settings(), {
            rValues$steps <- as.data.table(settings()$steps)
            rValues$output <- settings()$output
            rValues$brukerCalib <- settings()$brukerCalib
            triggerMSConvHOTUpdate()
        })
        
        makeMSConversionHOT <- function()
        {
            rValues$triggerMSConvHOTUpdate
            tab <- isolate(copy(rValues$steps))
            tab[, from := sapply(from, convTypeFormatToLabel)]
            tab[, to := sapply(to, convTypeFormatToLabel)]
            makeNewProjectHOT(tab, height = 250, maxRows = nrow(tab), columnSorting = FALSE)
        }
        moveSelectedConversions <- function(dir)
        {
            # NOTE: assume selection is a block range
            sel <- input$MSConversionHOT_select$select$rAll
            rValues$steps <- moveSelectedTabRows(dir, sel, rValues$steps)
            triggerMSConvHOTUpdate()
            moveHOTSel(session, dir, sel, ns("MSConversionHOT"))
        }
        selectDAMethod <- function(inputName)
        {
            dm <- rstudioapi::selectDirectory("Select DataAnalysis method")
            if (!is.null(dm))
            {
                if (!file.exists(file.path(dm, "DataAnalysis.method")))
                {
                    rstudioapi::showDialog("Invalid DataAnalysis method", "Please select a valid DataAnalysis method!", "")
                    dm <- NULL
                }
            }
            
            if (!is.null(dm))
                updateTextInput(session, inputName, value = dm)
        }
        
        observeEvent(input$addMSConversion, {
            showModal(modalDialog(
                title = "Add MS data conversion step",
                fillCol(
                    flex = 1,
                    height = 150,
                    width = 600,
                    fillCol(
                        selectInput(ns("convAlgo"), "Algorithm",
                                    c("ProteoWizard" = "pwiz", "OpenMS" = "openms", "Bruker" = "bruker",
                                      "IM collapse" = "im_collapse", "TIMSCONVERT" = "timsconvert"),
                                    multiple = FALSE, width = "100%")
                    ),
                    fillRow(
                        flex = c(1, NA, 1),
                        selectInput(ns("convFrom"), "From", groupConvTypesFormats("pwiz", "input"), multiple = FALSE,
                                    width = "100%"),
                        fillCol(width = "20px"),
                        selectInput(ns("convTo"), "To", groupConvTypesFormats("pwiz", "output"), multiple = FALSE,
                                    width = "100%")
                    )
                ),
                easyClose = TRUE,
                footer = tagList(
                    modalButton("Cancel"),
                    actionButton(ns("addMSConversionOK"), "OK")
                )
            ))
        })
        doObserveSelDir(input, session, "MSConversionPath", "MSConversionPathButton")
        observeEvent(input$configConvPaths, {
            showModal(modalDialog(
                title = "Set MS conversion output paths",
                fillCol(
                    flex = c(1, 1, 1, NA),
                    height = 275,
                    width = 400,
                    fileSelect(ns("MSConvOutputCentroid"), ns("MSConvOutputCentroidButton"), "Centroid",
                               rValues$output$centroid),
                    fileSelect(ns("MSConvOutputProfile"), ns("MSConvOutputProfileButton"), "Profile",
                               rValues$output$profile),
                    fileSelect(ns("MSConvOutputIMS"), ns("MSConvOutputIMSButton"), "IMS", rValues$output$ims),
                    fillCol(
                        height = 25,
                        textNote("Relative paths are relative to the project directory.")
                    )
                ),
                easyClose = TRUE,
                footer = tagList(
                    modalButton("Cancel"),
                    actionButton(ns("MSConvOutputOK"), "OK")
                )
            ))
        })
        doObserveSelDir(input, session, "MSConvOutputRaw", "MSConvOutputRawButton")
        doObserveSelDir(input, session, "MSConvOutputCentroid", "MSConvOutputCentroidButton")
        doObserveSelDir(input, session, "MSConvOutputProfile", "MSConvOutputProfileButton")
        doObserveSelDir(input, session, "MSConvOutputIMS", "MSConvOutputIMSButton")
        
        observeEvent(input$configBrukerCalib, {
            showModal(modalDialog(
                title = "Configure Bruker DataAnalysis m/z re-calibration",
                fillCol(
                    flex = NA,
                    width = 600,
                    fillCol(
                        height = 50,
                        checkboxInput(ns("doBrukerCalib"), "Perform m/z re-calibration",
                                      value = rValues$brukerCalib$enabled)
                    ),
                    conditionalPanel(
                        condition = "input.doBrukerCalib",
                        ns = ns,
                        fillCol(
                            flex = NA,
                            height = 90,
                            fillCol(
                                flex = NA,
                                height = 50,
                                conditionalPanel(
                                    condition = "output.ionization != \"both\"",
                                    ns = ns,
                                    fileSelect(ns("DAMethod"), ns("DAMethodButton"), "DataAnalysis method",
                                               rValues$brukerCalib$method),
                                ),
                                conditionalPanel(
                                    condition = "output.ionization == \"both\"",
                                    ns = ns,
                                    fillRow(
                                        fillCol(
                                            width = "95%",
                                            fileSelect(ns("DAMethodPos"), ns("DAMethodButtonPos"),
                                                       "DataAnalysis method (positive)",
                                                       rValues$brukerCalib$methodPos)
                                        ),
                                        fillCol(
                                            width = "95%",
                                            fileSelect(ns("DAMethodNeg"), ns("DAMethodButtonNeg"),
                                                       "DataAnalysis method (negative)",
                                                       rValues$brukerCalib$methodNeg)
                                        )
                                    )
                                )
                            ),
                            textNote("Leaving this blank will not set any method.")
                        )
                    ),
                    fillCol(
                        height = 40,
                        textNote(HTML("Only supported with bruker data and if DataAnalysis is installed.<br>This step is always executed before other conversion steps."))
                    )
                ),
                easyClose = TRUE,
                footer = tagList(
                    modalButton("Cancel"),
                    actionButton(ns("brukerCalibOK"), "OK")
                )
            ))
        })
        
        observeEvent(input$convAlgo, {
            f <- groupConvTypesFormats(input$convAlgo, "input")
            updateSelectInput(session, "convFrom", choices = f, selected = f[[1]][1])
            f <- groupConvTypesFormats(input$convAlgo, "output")
            updateSelectInput(session, "convTo", choices = f, selected = f[[1]][1])
        })
        observeEvent(input$DAMethodButton, selectDAMethod("DAMethod"))
        observeEvent(input$DAMethodButtonPos, selectDAMethod("DAMethodPos"))
        observeEvent(input$DAMethodButtonNeg, selectDAMethod("DAMethodNeg"))
        observeEvent(input$addMSConversionOK, {
            removeModal()
            rValues$steps <- rbind(rValues$steps, data.table(algorithm = input$convAlgo, from = input$convFrom,
                                                             to = input$convTo))
            triggerMSConvHOTUpdate()
        })
        observeEvent(input$MSConversionHOT, {
            # HACK: maxRows: make sure we don't have empty table as hot_to_r errors otherwise
            if (is.null(input$MSConversionHOT) || input$MSConversionHOT$params$maxRows == 0)
            {
                shinyjs::disable("removeMSConversionRows")
                shinyjs::disable("MSConversionRowsUp")
                shinyjs::disable("MSConversionRowsDown")
            }
        })
        observeEvent(input$MSConversionHOT_select$select$r, {
            sel <- input$MSConversionHOT_select$select$rAll
            e <- nrow(rValues$steps) > 0 && length(sel) > 0
            shinyjs::toggleState("removeMSConversionRows", condition = e)
            shinyjs::toggleState("MSConversionRowsUp", condition = e && min(sel) > 1)
            shinyjs::toggleState("MSConversionRowsDown", condition = e && max(sel) < nrow(rValues$steps))
        })
        observeEvent(input$removeMSConversionRows, {
            sel <- input$MSConversionHOT_select$select$rAll
            rValues$steps <- rValues$steps[-sel]
            triggerMSConvHOTUpdate()
        })
        observeEvent(input$MSConversionRowsUp, { moveSelectedConversions("up") })
        observeEvent(input$MSConversionRowsDown, { moveSelectedConversions("down") })
        observeEvent(input$MSConvOutputOK, {
            removeModal()
            rValues$output$centroid <- input$MSConvOutputCentroid
            rValues$output$profile <- input$MSConvOutputProfile
            rValues$output$ims <- input$MSConvOutputIMS
        })
        observeEvent(input$brukerCalibOK, {
            removeModal()
            rValues$brukerCalib <- list(enabled = input$doBrukerCalib, method = input$DAMethod,
                                        methodPos = input$DAMethodPos, methodNeg = input$DAMethodNeg)
        })
        
        output$MSConversionHOT <- rhandsontable::renderRHandsontable(makeMSConversionHOT())
        
        output <- exportShinyOutputVal(output, "ionization", ionization)
        
        list(
            valid = \() TRUE,
            settings = reactive(list(
                steps = rValues$steps,
                output = rValues$output,
                brukerCalib = rValues$brukerCalib
            ))
        )
    })
}
