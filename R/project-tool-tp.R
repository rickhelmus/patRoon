# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

getTPGenInputs <- function(isLib)
{
    ret <- c("Suspect list" = "suspects",
             "Screening results" = "screening")
    if (isLib)
        ret <- c(ret, "All library parents" = "all")
    return(ret)
}

newProjectTPUI <- function(id)
{
    ns <- NS(id)
    
    miniUI::miniContentPanel(
        fillCol(
            flex = NA,
            fillCol(
                height = 100,
                fillRow(
                    height = 75,
                    selectInput(ns("TPsAlgo"), "TP algorithm",
                                c("None" = "", "BioTransformer" = "biotransformer", "CTS" = "cts",
                                  "Library" = "library", "Logic" = "logic",
                                  "From formula annotations" = "ann_form",
                                  "From compound annotations" = "ann_comp"), width = "95%")
                ),
                conditionalPanel(
                    condition = "['biotransformer', 'cts', 'library'].includes(input.TPsAlgo)",
                    ns = ns,
                    fillRow(
                        height = 25,
                        checkboxInput(ns("TPDoMFDB"), "Generate TP MetFrag database", TRUE)
                    )
                )
            ),
            conditionalPanel(
                condition = "input.TPsAlgo != '' && input.TPsAlgo != 'logic'",
                ns = ns,
                fillRow(
                    selectInput(ns("TPGenInput"), "Parent input", getTPGenInputs(FALSE), width = "95%"),
                    conditionalPanel(
                        condition = "input.TPGenInput == \"suspects\"",
                        ns = ns,
                        fillRow(
                            width = "95%",
                            fileSelect(ns("TPSuspectList"), ns("TPSuspButton"), "Parent suspect list",
                                       placeholder = "Please specify parent suspect list")
                        )
                    )
                )
            )
        )
    )
}

newProjectTPServer <- function(id, hasSuspects, formulasAlgo, compoundsAlgo, settings)
{
    moduleServer(id, function(input, output, session)
    {
        observeEvent(settings(), {
            updateSelectInput(session, "TPsAlgo", selected = settings()$TPsAlgo)
            updateSelectInput(session, "TPGenInput", selected = settings()$TPGenInput)
            updateTextInput(session, "TPSuspectList", value = settings()$TPSuspectList)
            updateCheckboxInput(session, "TPDoMFDB", value = settings()$TPDoMFDB)
        })
        
        observeEvent(input$TPsAlgo, {
            updateSelectInput(inputId = "TPGenInput", choices = getTPGenInputs(input$TPsAlgo == "library"))
        })
        
        observeEvent(input$TPGenInput, {
            if (input$TPGenInput == "screening" && !hasSuspects())
            {
                rstudioapi::showDialog("Enable suspect screening",
                                       "This requires a workflow with suspect screening. Please configure in the Features tab.", "")
                updateSelectInput(inputId = "TPGenInput", selected = "suspects")
            }
        })
        
        observeEvent(input$TPSuspButton, selectSuspList(session, "TPSuspectList"))
        
        list(
            valid = reactive({
                if ((nzchar(input$TPsAlgo) && input$TPsAlgo != "logic") && input$TPGenInput == "suspects" &&
                    !nzchar(input$TPSuspectList))
                {
                    list(title = "No suspect list", msg = "Please select a suspect list!")
                }
                else if ((input$TPsAlgo == "ann_form" && !nzchar(formulasAlgo())) ||
                         (input$TPsAlgo == "ann_comp" && !nzchar(compoundsAlgo())))
                {
                    fc <- if (input$TPsAlgo == "ann_form") "formula" else "compound"
                    list(
                        title = sprintf("No %s annotation algorithm set", fc),
                        msg = sprintf("Please set a %s annotation algorithm in the Annotation tab when using the %s TP algorithm",
                                      fc, input$TPsAlgo)
                    )
                }
                else
                    TRUE
            }),
            settings = reactive(list(
                TPsAlgo = input$TPsAlgo,
                TPGenInput = input$TPGenInput,
                TPSuspectList = input$TPSuspectList,
                TPDoMFDB = input$TPDoMFDB
            ))
        )
    })
}

defaultTPSettings <- function()
{
    return(list(TPsAlgo = "", TPGenInput = "suspects", TPSuspectList = "", TPDoMFDB = TRUE))
}

upgradeTPsSettings <- function(settings)
{
    # NOTE: this updates from first file version
    ret <- modifyList(defaultTPSettings(), settings[c("TPGenInput", "TPSuspectList", "TPDoMFDB")])
    if (settings$doTPs)
        ret$TPsAlgo <- tolower(settings$TPGen)
    return(ret)
}
