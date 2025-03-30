newProjectFeaturesUI <- function(id)
{
    ns <- NS(id)
    
    htmltools::tagList(
        miniUI::miniContentPanel(
            fillCol(
                flex = NA,
                fillCol(
                    height = 90,
                    flex = NA,
                    fillRow(
                        selectInput(ns("featFinder"), "Feature finder", c("OpenMS", "XCMS", "enviPick", "SIRIUS", "KPIC2",
                                                                          "Bruker DataAnalysis" = "Bruker"),
                                    multiple = FALSE, width = "95%"),
                        fillCol(
                            flex = c(1, NA),
                            height = 90,
                            selectInput(ns("featGrouper"), "Feature grouper", c("OpenMS", "XCMS", "KPIC2", "SIRIUS"),
                                        multiple = FALSE, width = "100%"),
                            conditionalPanel(
                                condition = "input.featGrouper == \"SIRIUS\"",
                                ns = ns,
                                textNote(HTML("This will always find <b>and</b> group features with SIRIUS."))
                            )
                        )
                    )
                ),
                fillCol(
                    height = 125,
                    flex = NA,
                    conditionalPanel(
                        condition = "output.ionization != \"both\"",
                        ns = ns,
                        fileSelect(ns("suspectList"), ns("suspectListButton"), "Suspect list",
                                   placeholder = "Leave empty for no suspect screening")
                    ),
                    conditionalPanel(
                        condition = "output.ionization == \"both\"",
                        ns = ns,
                        fillRow(
                            height = 60,
                            fillCol(
                                width = "95%",
                                fileSelect(ns("suspectListPos"), ns("suspectListButtonPos"), "Suspect list (positive)",
                                           placeholder = "Leave empty for no suspect screening")
                            ),
                            fillCol(
                                width = "95%",
                                fileSelect(ns("suspectListNeg"), ns("suspectListButtonNeg"), "Suspect list (negative)",
                                           placeholder = "Leave empty if same as positive")
                            )
                        )
                    ),
                    fillRow(
                        height = 50,
                        checkboxInput(ns("exSuspList"), "Example suspect list(s)")
                    )
                )
            )
        ),
        miniUI::miniButtonBlock(
            actionButton(ns("advancedOpen"), "Advanced", icon = icon("cog"))
        )
    )
}

newProjectFeaturesServer <- function(id, ionization, settings)
{
    ns <- NS(id)
    
    rangeNumeric <- function(id, label, value, minVal = 0, maxVal = 0, ...)
    {
        fillRow(
            numericInput(paste0(id, "-min"), paste("Min.", label), value = value[1], ..., width = "95%"),
            numericInput(paste0(id, "-max"), paste("Max.", label), value = value[2], ..., width = "100%")
        )
    }
    
    moduleServer(id, function(input, output, session)
    {
        rValues <- reactiveValues(advanced = defaultFeaturesSettings()$advanced)
        
        observeEvent(settings(), {
            updateSelectInput(session, "featFinder", selected = settings()$featFinder)
            updateSelectInput(session, "featGrouper", selected = settings()$featGrouper)
            updateTextInput(session, "suspectList", value = settings()$suspectList)
            updateTextInput(session, "suspectListPos", value = settings()$suspectListPos)
            updateTextInput(session, "suspectListNeg", value = settings()$suspectListNeg)
            updateCheckboxInput(session, "exSuspList", value = settings()$exSuspList)
            rValues$advanced <- settings()$advanced
        })

        observeEvent(input$suspectListButton, selectSuspList(session, "suspectList"))
        observeEvent(input$suspectListButtonPos, selectSuspList(session, "suspectListPos"))
        observeEvent(input$suspectListButtonNeg, selectSuspList(session, "suspectListNeg"))
        observeEvent(input$exSuspList, {
            for (id in c("suspectList", "suspectListButton", "suspectListPos", "suspectListButtonPos",
                         "suspectListNeg", "suspectListButtonNeg"))
                shinyjs::toggleState(id, !input$exSuspList)
        })
        
        observeEvent(input$ISTDListButton, selectSuspList(session, "ISTDList"))
        observeEvent(input$ISTDListButtonPos, selectSuspList(session, "ISTDListPos"))
        observeEvent(input$ISTDListButtonNeg, selectSuspList(session, "ISTDListNeg"))
        
        observeEvent(input$featGrouper, {
            if (ionization() == "both" && input$featGrouper == "SIRIUS")
            {
                rstudioapi::showDialog("Not supported", "Grouping features with SIRIUS is currently not supported with sets", "")
                updateSelectInput(inputId = "featGrouper", selected = "OpenMS")
            }
            else
                shinyjs::toggleState("featFinder", input$featGrouper != "SIRIUS")
        })
        
        observeEvent(input$advancedOpen, {
            showModal(modalDialog(
                title = "Advanced feature settings",
                fillCol(
                    flex = NA,
                    width = 500,
                    height = 400,
                    style = "overflow-y: auto;",
                    fillCol(
                        flex = NA,
                        height = 70,
                        strong("Post-Filtering of feature groups"),
                        textNote("Set below values to zero to disable a particular filter.")
                    ),
                    fillCol(
                        height = 325,
                        fillRow(
                            numericInput(ns("preIntThr"), "Pre-Intensity threshold", rValues$advanced$preIntThr, 0,
                                         step = 100, width = "95%"),
                            numericInput(ns("intThr"), "Intensity threshold", rValues$advanced$intThr, 0, step = 1000,
                                         width = "100%")
                        ),
                        fillRow(
                            numericInput(ns("repAbundance"), "Min. replicate abundance (relative)",
                                         rValues$advanced$repAbundance, 0, 1.0, 0.1, width = "95%"),
                            numericInput(ns("maxRepRSD"), "Max. replicate intensity RSD", rValues$advanced$maxRepRSD, 0,
                                         step = 0.1, width = "100%")
                        ),
                        fillRow(
                            numericInput(ns("blankThr"), "Min. blank threshold", rValues$advanced$blankThr, 0, step = 1,
                                         width = "95%"),
                            checkboxInput(ns("removeBlanks"), "Discard blanks after filtering",
                                          rValues$advanced$removeBlanks)
                        ),
                        rangeNumeric(ns("retention"), "retention time (s)", step = 10,
                                     value = rValues$advanced$retention),
                        rangeNumeric(ns("mz"), "m/z", step = 10, value = rValues$advanced$mz)
                    ),
                    hr(),
                    fillCol(
                        height = 125,
                        flex = NA,
                        
                        selectInput(ns("featNorm"), "Feature normalization", c("None" = "none",
                                                                               "Internal standard" = "istd",
                                                                               "Internal standard concentration" = "conc",
                                                                               "TIC" = "tic"),
                                    rValues$advanced$featNorm),
                        checkboxInput(ns("groupNorm"), "Group normalization", rValues$advanced$groupNorm),
                        
                        conditionalPanel(
                            condition = "input.featNorm == \"istd\" && output.ionization != \"both\"",
                            ns = ns,
                            fileSelect(ns("ISTDList"), ns("ISTDListButton"), "Internal standard list",
                                       rValues$advanced$ISTDList, placeholder = "Leave empty for example list")
                        ),
                        conditionalPanel(
                            condition = "input.featNorm == \"istd\" && output.ionization == \"both\"",
                            ns = ns,
                            fillRow(
                                height = 60,
                                fillCol(
                                    width = "95%",
                                    fileSelect(ns("ISTDListPos"), ns("ISTDListButtonPos"),
                                               "Internal standard list (positive)",
                                               rValues$advanced$ISTDListPos,
                                               placeholder = "Leave empty for example list")
                                ),
                                fillCol(
                                    width = "95%",
                                    fileSelect(ns("ISTDListNeg"), ns("ISTDListButtonNeg"),
                                               "Internal standard list (negative)",
                                               rValues$advanced$ISTDListNeg,
                                               placeholder = "Leave empty if same as positive")
                                )
                            )
                        ),
                        conditionalPanel(
                            condition = "input.featNorm == \"istd\" || input.featNorm == \"conc\"",
                            ns = ns,
                            textNote(HTML("Please make sure that the <i>norm_conc</i> column of the analysis information is set."))
                        )
                    )
                ),
                easyClose = TRUE,
                footer = tagList(
                    modalButton("Cancel"),
                    actionButton(ns("advancedOK"), "OK")
                )
            ))
        })
        
        observeEvent(input$advancedOK, {
            removeModal()
            rValues$advanced <- list(
                preIntThr = input$preIntThr,
                intThr = input$intThr,
                repAbundance = input$repAbundance,
                maxRepRSD = input$maxRepRSD,
                blankThr = input$blankThr,
                removeBlanks = input$removeBlanks,
                retention = c(input[["retention-min"]], input[["retention-max"]]),
                mz = c(input[["mz-min"]], input[["mz-max"]]),
                featNorm = input$featNorm,
                groupNorm = input$groupNorm,
                ISTDList = input$ISTDList,
                ISTDListPos = input$ISTDListPos,
                ISTDListNeg = input$ISTDListNeg
            )
        })
        
        output <- exportShinyOutputVal(output, "ionization", ionization)
        
        list(
            valid = \() TRUE,
            settings = reactive(list(
                featFinder = input$featFinder,
                featGrouper = input$featGrouper,
                suspectList = input$suspectList,
                suspectListPos = input$suspectListPos,
                suspectListNeg = input$suspectListNeg,
                exSuspList = input$exSuspList,
                advanced = rValues$advanced
            ))
        )
    })
}

defaultFeaturesSettings <- function()
{
    return(list(
        featFinder = "OpenMS",
        featGrouper = "OpenMS",
        suspectList = "",
        suspectListPos = "",
        suspectListNeg = "",
        exSuspList = FALSE,
        advanced = list(
            preIntThr = 1E2,
            intThr = 1E4,
            repAbundance = 1,
            maxRepRSD = 0.75,
            blankThr = 5,
            removeBlanks = TRUE,
            retention = c(0, 0),
            mz = c(0, 0),
            featNorm = "none",
            groupNorm = FALSE,
            ISTDList = "",
            ISTDListPos = "",
            ISTDListNeg = ""
        )
    ))
}

upgradeFeaturesSettings <- function(settings)
{
    # NOTE: this updates from first file version
    ret <- modifyList(defaultFeaturesSettings(),
                      settings[c("featFinder", "featGrouper", "suspectList", "suspectListPos", "suspectListNeg",
                                 "exSuspList")])
    ret$retention <- c(settings[["retention-min"]], settings[["retention-max"]])
    ret$mz <- c(settings[["mz-min"]], settings[["mz-max"]])
    ret$advanced <- settings[c("preIntThr", "intThr", "repAbundance", "maxRepRSD", "blankThr", "removeBlanks",
                               "retention", "mz", "featNorm", "groupNorm", "ISTDList", "ISTDListPos", "ISTDListNeg")]
    return(ret)
}
