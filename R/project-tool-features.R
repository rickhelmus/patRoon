newProjectFeaturesUI <- function(id)
{
    ns <- NS(id)
    
    miniUI::miniContentPanel(
        fillCol(
            flex = NA,
            fillRow(
                height = 90,
                selectInput(ns("featAlgo"), "Feature finding algorithm",
                            c("OpenMS", "XCMS", "enviPick", "SIRIUS", "KPIC2", EICs = "EIC",
                              "Bruker DataAnalysis" = "Bruker"), multiple = FALSE, width = "95%"),
                conditionalPanel(
                    condition = "input.featAlgo == \"EIC\"",
                    ns = ns,
                    style = "margin-top: 25px;", # align with select label: https://stackoverflow.com/a/70543848
                    actionButton(ns("featEICParams"), "", icon = icon("cog"), title = "EIC parameters")
                )
            ),
            fillRow(
                height = 90,
                fillCol(
                    flex = c(1, NA),
                    selectInput(ns("fGroupsAlgo"), "Feature grouping algorithm", c("OpenMS", "XCMS", "KPIC2", "SIRIUS"),
                                multiple = FALSE, width = "95%"),
                    conditionalPanel(
                        condition = "input.fGroupsAlgo == \"SIRIUS\"",
                        ns = ns,
                        textNote(HTML("This will always find <b>and</b> group features with SIRIUS."))
                    )
                ),
                fillCol(
                    style = "margin-top: 25px;",
                    actionButton(ns("fGroupsAdv"), "", icon = icon("cog"), title = "Advanced settings")
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
                                       placeholder = "Leave empty for same as positive")
                        )
                    )
                ),
                fillRow(
                    height = 50,
                    checkboxInput(ns("exSuspList"), "Example suspect list(s)")
                )
            )
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
    
    getEICIMSMethods <- function(methodMZ)
    {
        switch(methodMZ,
               bins = "bins",
               suspects = c("bins", "suspects"),
               ms2 = c("bins", "ms2"))
    }
    
    moduleServer(id, function(input, output, session)
    {
        rValues <- reactiveValues(
            featEICParams = defaultFeaturesSettings()$featEICParams,
            fGroupsAdv = defaultFeaturesSettings()$fGroupsAdv
        )
        
        observeEvent(settings(), {
            updateSelectInput(session, "featAlgo", selected = settings()$featAlgo)
            rValues$featEICParams <- settings()$featEICParams
            updateSelectInput(session, "fGroupsAlgo", selected = settings()$fGroupsAlgo)
            rValues$fGroupsAdv <- settings()$fGroupsAdv
            updateTextInput(session, "suspectList", value = settings()$suspectList)
            updateTextInput(session, "suspectListPos", value = settings()$suspectListPos)
            updateTextInput(session, "suspectListNeg", value = settings()$suspectListNeg)
            updateCheckboxInput(session, "exSuspList", value = settings()$exSuspList)
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
        
        observeEvent(input$fGroupsAlgo, {
            if (ionization() == "both" && input$fGroupsAlgo == "SIRIUS")
            {
                rstudioapi::showDialog("Not supported", "Grouping features with SIRIUS is currently not supported with sets", "")
                updateSelectInput(inputId = "fGroupsAlgo", selected = "OpenMS")
            }
            else
                shinyjs::toggleState("featAlgo", input$fGroupsAlgo != "SIRIUS")
        })

        observeEvent(input$featEICParams, {
            # UNDONE: update methodIMS options based on methodMZ selection
            showModal(modalDialog(
                title = "EIC feature parameters",
                fillCol(
                    flex = NA,
                    width = 600,
                    fillRow(
                        height = 75,
                        selectInput(ns("methodMZ"), "EIC method (m/z)", c("bins", "suspects", "ms2"),
                                    selected = rValues$featEICParams$methodMZ, width = "95%"),
                        # UNDONE: only in post-IMS workflows
                        selectInput(ns("methodIMS"), "EIC method (IMS)",
                                    getEICIMSMethods(rValues$featEICParams$methodMZ),
                                    selected = rValues$featEICParams$methodIMS, width = "95%"),
                    ),
                    conditionalPanel(
                        condition = "input.methodMZ == \"suspects\" && output.ionization != \"both\"",
                        ns = ns,
                        fillRow(
                            height = 75,
                            width = "97.5%", # 97.5% is from the 2x 95% widths for the method selectInputs
                            fileSelect(ns("featEICSuspectList"), ns("featEICSuspectListButton"), "Suspect list",
                                       rValues$featEICParams$suspects,
                                       placeholder = "Empty for same as Features tab")
                        )
                    ),
                    conditionalPanel(
                        condition = "input.methodMZ == \"suspects\" && output.ionization == \"both\"",
                        ns = ns,
                        fillRow(
                            height = 75,
                            fillCol(
                                width = "95%",
                                fileSelect(ns("featEICSuspectListPos"), ns("featEICSuspectListButtonPos"),
                                           "Suspect list (positive)", rValues$featEICParams$suspectsPos,
                                           placeholder = "Empty for same as Features tab")
                            ),
                            fillCol(
                                width = "95%",
                                fileSelect(ns("featEICSuspectListNeg"), ns("featEICSuspectListButtonNeg"),
                                           "Suspect list (negative)", rValues$featEICParams$suspectsNeg,
                                           placeholder = "Empty for same as positive")
                            )
                        )
                    ),
                    fillRow(
                        height = 75,
                        # NOTE: 97.5% is from the 2x 95% widths for the method selectInputs
                        selectInput(ns("peaksAlgo"), "Peak detection algorithm",
                                    c(Dietrich = "dietrich", "OpenMS" = "openms", "XCMS3" = "xcms3", "enviPick" = "envipick"),
                                    selected = rValues$featEICParams$peaksAlgo, width = "97.5%")
                    )
                ),
                easyClose = TRUE,
                footer = tagList(
                    modalButton("Cancel"),
                    actionButton(ns("featEICParamsOK"), "OK")
                )
            ))
        })
        
        observeEvent(input$methodMZ, {
            ch <- getEICIMSMethods(input$methodMZ)
            sel <- input$methodIMS
            if (!sel %in% ch)
                sel <- ch[1]
            updateSelectInput(session, "methodIMS", choices = ch, selected = sel)
        })
        
        observeEvent(input$featEICParamsOK, {
            removeModal()
            rValues$featEICParams <- list(
                methodMZ = input$methodMZ,
                methodIMS = input$methodIMS,
                suspects = input$featEICSuspectList,
                suspectsPos = input$featEICSuspectListPos,
                suspectsNeg = input$featEICSuspectListNeg,
                peaksAlgo = input$peaksAlgo
            )
        })
        
        observeEvent(input$fGroupsAdv, {
            showModal(modalDialog(
                title = "Advanced feature group settings",
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
                            numericInput(ns("preIntThr"), "Pre-Intensity threshold", rValues$fGroupsAdv$preIntThr, 0,
                                         step = 100, width = "95%"),
                            numericInput(ns("intThr"), "Intensity threshold", rValues$fGroupsAdv$intThr, 0, step = 1000,
                                         width = "100%")
                        ),
                        fillRow(
                            numericInput(ns("repAbundance"), "Min. replicate abundance (relative)",
                                         rValues$fGroupsAdv$repAbundance, 0, 1.0, 0.1, width = "95%"),
                            numericInput(ns("maxRepRSD"), "Max. replicate intensity RSD", rValues$fGroupsAdv$maxRepRSD, 0,
                                         step = 0.1, width = "100%")
                        ),
                        fillRow(
                            numericInput(ns("blankThr"), "Min. blank threshold", rValues$fGroupsAdv$blankThr, 0, step = 1,
                                         width = "95%"),
                            checkboxInput(ns("removeBlanks"), "Discard blanks after filtering",
                                          rValues$fGroupsAdv$removeBlanks)
                        ),
                        rangeNumeric(ns("retention"), "retention time (s)", step = 10,
                                     value = rValues$fGroupsAdv$retention),
                        rangeNumeric(ns("mz"), "m/z", step = 10, value = rValues$fGroupsAdv$mz)
                    ),
                    hr(),
                    fillCol(
                        height = 125,
                        flex = NA,
                        
                        selectInput(ns("featNorm"), "Feature normalization", c("None" = "none",
                                                                               "Internal standard" = "istd",
                                                                               "Internal standard concentration" = "conc",
                                                                               "TIC" = "tic"),
                                    rValues$fGroupsAdv$featNorm),
                        checkboxInput(ns("groupNorm"), "Group normalization", rValues$fGroupsAdv$groupNorm),
                        
                        conditionalPanel(
                            condition = "input.featNorm == \"istd\" && output.ionization != \"both\"",
                            ns = ns,
                            fileSelect(ns("ISTDList"), ns("ISTDListButton"), "Internal standard list",
                                       rValues$fGroupsAdv$ISTDList, placeholder = "Leave empty for example list")
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
                                               rValues$fGroupsAdv$ISTDListPos,
                                               placeholder = "Leave empty for example list")
                                ),
                                fillCol(
                                    width = "95%",
                                    fileSelect(ns("ISTDListNeg"), ns("ISTDListButtonNeg"),
                                               "Internal standard list (negative)",
                                               rValues$fGroupsAdv$ISTDListNeg,
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
                    actionButton(ns("fGroupsAdvOK"), "OK")
                )
            ))
        })
        
        observeEvent(input$fGroupsAdvOK, {
            removeModal()
            rValues$fGroupsAdv <- list(
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
                featAlgo = input$featAlgo,
                fGroupsAlgo = input$fGroupsAlgo,
                suspectList = input$suspectList,
                suspectListPos = input$suspectListPos,
                suspectListNeg = input$suspectListNeg,
                exSuspList = input$exSuspList,
                fGroupsAdv = rValues$fGroupsAdv
            ))
        )
    })
}

defaultFeaturesSettings <- function()
{
    return(list(
        featAlgo = "OpenMS",
        fGroupsAlgo = "OpenMS",
        suspectList = "",
        suspectListPos = "",
        suspectListNeg = "",
        exSuspList = FALSE,
        featEICParams = list(
            methodMZ = "bins",
            methodIMS = "bins",
            suspects = "",
            suspectsPos = "",
            suspectsNeg = "",
            peaksAlgo = "dietrich"
        ),
        fGroupsAdv = list(
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
                      settings[c("suspectList", "suspectListPos", "suspectListNeg", "exSuspList")])
    ret$featAlgo <- settings$featFinder
    ret$fGroupsAlgo <- settings$featGrouper
    ret$retention <- c(settings[["retention-min"]], settings[["retention-max"]])
    ret$mz <- c(settings[["mz-min"]], settings[["mz-max"]])
    ret$fGroupsAdv <- settings[c("preIntThr", "intThr", "repAbundance", "maxRepRSD", "blankThr", "removeBlanks",
                                 "retention", "mz", "featNorm", "groupNorm", "ISTDList", "ISTDListPos", "ISTDListNeg")]
    return(ret)
}
