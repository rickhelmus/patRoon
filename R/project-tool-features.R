getFeatAlgoSelections <- function()
{
    c("OpenMS", "XCMS", "enviPick", "SIRIUS", "KPIC2", "piek", "Bruker DataAnalysis" = "Bruker")
}

newProjectFeaturesUI <- function(id)
{
    ns <- NS(id)
 
    getPeakAlgoSelections <- function() c("piek", openms = "OpenMS", xcms3 = "XCMS", envipick = "enviPick")
    
    tagList(
        miniUI::miniContentPanel(
            fillCol(
                flex = NA,
                fillRow(
                    flex = c(1, NA, 1),
                    height = 130,
                    fillCol(
                        fillRow(
                            flex = c(1, NA),
                            selectInput(ns("featAlgo"), "Feature finding algorithm", getFeatAlgoSelections(),
                                        width = "100%"),
                            conditionalPanel(
                                condition = "input.featAlgo == \"piek\"",
                                ns = ns,
                                # align with select label: https://stackoverflow.com/a/70543848
                                style = "margin-top: 25px; margin-left: 10px;",
                                actionButton(ns("piekParams"), "", icon = icon("cog"), title = "piek parameters")
                            )
                        ),
                        fillCol(
                            flex = c(1, NA),
                            selectInput(ns("fGroupsAlgo"), "Feature grouping algorithm",
                                        c("OpenMS", "XCMS", "KPIC2", "SIRIUS"), multiple = FALSE, width = "100%"),
                            conditionalPanel(
                                condition = "input.fGroupsAlgo == \"SIRIUS\"",
                                ns = ns,
                                textNote(HTML("This will always find <b>and</b> group features with SIRIUS."))
                            )
                        )
                    ),
                    div(style = "width: 25px;"),
                    conditionalPanel(
                        condition = "output.IMSMode == \"post\"",
                        ns = ns,
                        fillCol(
                            height = 130,
                            selectInput(ns("IMSPeaksMob"), "Mobility peak detection", getPeakAlgoSelections(),
                                        width = "90%"),
                            selectInput(ns("IMSPeaksChrom"), "Mobility feature peak re-integration",
                                        c(None = "", getPeakAlgoSelections()), width = "90%")
                        )
                    )
                ),
                hr(),
                fillCol(
                    flex = NA,
                    height = 130,
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
                    ),
                    conditionalPanel(
                        condition = "output.IMSMode != \"none\"",
                        ns = ns,
                        fillRow(
                            height = 75,
                            radioButtons(ns("IMSSuspCCSPred"), "CCS prediction", getCCSPredSelections(), inline = TRUE)
                        )
                    )
                )
            )
        ),
        miniUI::miniButtonBlock(
            actionButton(ns("fGroupsAdv"), "Advanced", icon = icon("cog"))
        )
    )
}

newProjectFeaturesServer <- function(id, ionization, IMSMode, settings)
{
    ns <- NS(id)
    
    rangeNumeric <- function(id, label, value, minVal = 0, maxVal = 0, ...)
    {
        fillRow(
            numericInput(paste0(id, "-min"), paste("Min.", label), value = value[1], ..., width = "95%"),
            numericInput(paste0(id, "-max"), paste("Max.", label), value = value[2], ..., width = "100%")
        )
    }
    
    getPiekEICIMSMethods <- function(methodMZ)
    {
        switch(methodMZ,
               bins = "bins",
               suspects = c("bins", "suspects"),
               ms2 = c("bins", "ms2"))
    }
    
    moduleServer(id, function(input, output, session)
    {
        rValues <- reactiveValues(
            piekParams = defaultFeaturesSettings()$piekParams,
            fGroupsAdv = defaultFeaturesSettings()$fGroupsAdv
        )
        
        observeEvent(settings(), {
            updateSelectInput(session, "featAlgo", selected = settings()$featAlgo)
            rValues$piekParams <- settings()$piekParams
            updateSelectInput(session, "fGroupsAlgo", selected = settings()$fGroupsAlgo)
            updateSelectInput(session, "IMSPeaksMob", selected = settings()$IMSPeaksMob)
            updateSelectInput(session, "IMSPeaksChrom", selected = settings()$IMSPeaksChrom)
            rValues$fGroupsAdv <- settings()$fGroupsAdv
            updateTextInput(session, "suspectList", value = settings()$suspects$single)
            updateTextInput(session, "suspectListPos", value = settings()$suspects$sets$pos)
            updateTextInput(session, "suspectListNeg", value = settings()$suspects$sets$neg)
            updateCheckboxInput(session, "exSuspList", value = settings()$exSuspList)
            updateRadioButtons(session, "IMSSuspCCSPred", selected = settings()$IMSSuspCCSPred)
        })
        
        observeEvent(IMSMode(), {
            sels <- getFeatAlgoSelections()
            curSel <- input$featAlgo
            if (IMSMode() == "direct")
            {
                sels <- sels[sels == "piek"]
            }
            if (!curSel %in% sels)
                curSel <- sels[1]
            updateSelectInput(session, "featAlgo", choices = sels, selected = curSel)
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

        observeEvent(input$piekParams, {
            showModal(modalDialog(
                title = "piek feature parameters",
                fillCol(
                    flex = NA,
                    width = 600,
                    fillRow(
                        height = 75,
                        selectInput(ns("methodMZ"), "EIC method (m/z)", c("bins", "suspects", "ms2"),
                                    selected = rValues$piekParams$methodMZ, width = "95%"),
                        
                    ),
                    conditionalPanel(
                        condition = "output.IMSMode == \"direct\"",
                        ns = ns,
                        fillRow(
                            height = 75,
                            selectInput(ns("methodIMS"), "EIC method (IMS)",
                                        getPiekEICIMSMethods(rValues$piekParams$methodMZ),
                                        selected = rValues$piekParams$methodIMS, width = "95%")
                        )
                    ),
                    conditionalPanel(
                        condition = "input.methodMZ == \"suspects\" && output.ionization != \"both\"",
                        ns = ns,
                        fillRow(
                            height = 75,
                            width = "95%",
                            fileSelect(ns("piekSuspectList"), ns("piekSuspectListButton"), "Suspect list",
                                       rValues$piekParams$suspects$single,
                                       placeholder = "Empty for same as Features tab")
                        )
                    ),
                    conditionalPanel(
                        condition = "input.methodMZ == \"suspects\" && output.ionization == \"both\"",
                        ns = ns,
                        fillRow(
                            height = 75,
                            fillCol(
                                width = "90%",
                                fileSelect(ns("piekSuspectListPos"), ns("piekSuspectListButtonPos"),
                                           "Suspect list (positive)", rValues$piekParams$suspects$sets$pos,
                                           placeholder = "Empty for same as Features tab")
                            ),
                            fillCol(
                                width = "90%",
                                fileSelect(ns("piekSuspectListNeg"), ns("piekSuspectListButtonNeg"),
                                           "Suspect list (negative)", rValues$piekParams$suspects$sets$neg,
                                           placeholder = "Empty for same as positive")
                            )
                        )
                    ),
                    fillRow(
                        height = 75,
                        selectInput(ns("peaksAlgo"), "Peak detection algorithm",
                                    c("piek", "OpenMS" = "openms", "XCMS3" = "xcms3", "enviPick" = "envipick"),
                                    selected = rValues$piekParams$peaksAlgo, width = "95%")
                    )
                ),
                easyClose = TRUE,
                footer = tagList(
                    modalButton("Cancel"),
                    actionButton(ns("piekParamsOK"), "OK")
                )
            ))
        })
        
        observeEvent(input$methodMZ, {
            ch <- getPiekEICIMSMethods(input$methodMZ)
            sel <- input$methodIMS
            if (!sel %in% ch)
                sel <- ch[1]
            updateSelectInput(session, "methodIMS", choices = ch, selected = sel)
        })
        
        observeEvent(input$piekParamsOK, {
            removeModal()
            rValues$piekParams <- list(
                methodMZ = input$methodMZ,
                methodIMS = input$methodIMS,
                suspects = list(single = input$piekSuspectList,
                                sets = list(pos = input$piekSuspectListPos, neg = input$piekSuspectListNeg)),
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
                                       rValues$fGroupsAdv$ISTDLists$single, placeholder = "Leave empty for example list")
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
                                               rValues$fGroupsAdv$ISTDLists$sets$pos,
                                               placeholder = "Leave empty for example list")
                                ),
                                fillCol(
                                    width = "95%",
                                    fileSelect(ns("ISTDListNeg"), ns("ISTDListButtonNeg"),
                                               "Internal standard list (negative)",
                                               rValues$fGroupsAdv$ISTDLists$sets$neg,
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
                ISTDLists = list(single = input$ISTDList,
                                 sets = list(pos = input$ISTDListPos, neg = input$ISTDListNeg))
            )
        })
        
        output <- exportShinyOutputVal(output, "ionization", ionization)
        output <- exportShinyOutputVal(output, "IMSMode", IMSMode)
        
        list(
            valid = reactive({
                suspEmpty <- (ionization() != "both" && !nzchar(input$suspectList)) ||
                    (ionization() == "both" && !nzchar(input$suspectListPos))
                suspPiekEmpty <- (ionization() != "both" && !nzchar(rValues$piekParams$suspects$single)) ||
                    (ionization() == "both" && !nzchar(rValues$piekParams$suspects$pos))
                if (input$featAlgo == "piek" && rValues$piekParams$methodMZ == "suspects" && suspEmpty &&
                    suspPiekEmpty)
                {
                    list(title = "No suspect list",
                         msg = "Please select suspect data for the piek feature finding algorithm!")
                }
                else
                    TRUE
            }),
            settings = reactive(list(
                featAlgo = input$featAlgo,
                piekParams = rValues$piekParams,
                fGroupsAlgo = input$fGroupsAlgo,
                IMSPeaksMob = input$IMSPeaksMob,
                IMSPeaksChrom = input$IMSPeaksChrom,
                suspects = list(single = input$suspectList,
                                sets = list(pos = input$suspectListPos, neg = input$suspectListNeg)),
                exSuspList = input$exSuspList,
                IMSSuspCCSPred = input$IMSSuspCCSPred,
                fGroupsAdv = rValues$fGroupsAdv
            ))
        )
    })
}

defaultFeaturesSettings <- function()
{
    return(list(
        featAlgo = "OpenMS",
        piekParams = list(
            methodMZ = "bins",
            methodIMS = "bins",
            suspects = list(single = "", sets = list(pos = "", neg = "")),
            peaksAlgo = "piek"
        ),
        fGroupsAlgo = "OpenMS",
        IMSPeaksMob = "piek",
        IMSPeaksChrom = "piek",
        suspects = list(single = "", sets = list(pos = "", neg = "")),
        exSuspList = FALSE,
        IMSSuspCCSPred = "none",
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
            ISTDLists = list(single = "", sets = list(pos = "", neg = ""))
        )
    ))
}

upgradeFeaturesSettings <- function(settings)
{
    # NOTE: this updates from first file version
    ret <- defaultFeaturesSettings()
    ret$suspects <- list(single = settings$suspectList,
                         sets = list(pos = settings$suspectListPos, neg = settings$suspectListNeg))
    ret$exSuspList <- settings$exSuspList
    ret$featAlgo <- settings$featFinder
    ret$fGroupsAlgo <- settings$featGrouper
    ret$retention <- c(settings[["retention-min"]], settings[["retention-max"]])
    ret$mz <- c(settings[["mz-min"]], settings[["mz-max"]])
    ret$fGroupsAdv <- settings[c("preIntThr", "intThr", "repAbundance", "maxRepRSD", "blankThr", "removeBlanks",
                                 "featNorm", "groupNorm")]
    ret$fGroupsAdv$retention <- c(settings[["retention-min"]], settings[["retention-max"]])
    ret$fGroupsAdv$mz <- c(settings[["mz-min"]], settings[["mz-max"]])
    ret$fGroupsAdv$ISTDLists <- list(single = settings$ISTDList,
                                     sets = list(pos = settings$ISTDListPos, neg = settings$ISTDListNeg))
    return(ret)
}
