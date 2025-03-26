newProjectGeneralUI <- function(id, destPath)
{
    ns <- NS(id)
    tagList(
        miniUI::miniContentPanel(
            fillCol(
                flex = NA,
                fileSelect(ns("destination"), ns("destinationButton"), "Project destination",
                           if (is.null(destPath)) "~/" else destPath),
                br(),
                fillRow(
                    height = 100,
                    radioButtons(ns("outputScriptTo"), "Insert code into", c("New file" = "newFile",
                                                                             "Current file" = "curFile")),
                    conditionalPanel(
                        condition = "input.outputScriptTo == \"newFile\"",
                        ns = ns,
                        textInput(ns("scriptFile"), "Script file", "process.R", width = "80%"),
                        checkboxInput(ns("createRStudioProj"), "Create (and open) RStudio project", value = TRUE)
                    )
                ),
                br(),
                fillRow(
                    height = 100,
                    radioButtons(ns("ionization"), "Ionization", c("positive", "negative", "both (sets)" = "both"))
                )
            )
        ),
        miniUI::miniButtonBlock(
            actionButton(ns("loadParams"), "Load parameters"),
            actionButton(ns("saveParams"), "Save parameters")
        )
    )
}

newProjectGeneralServer <- function(id, settings)
{
    moduleServer(id, function(input, output, session)
    {
        doObserveSelDir(input, session, "destination", "destinationButton")
        
        observeEvent(settings(), {
            # UNDONE: update destination?
            updateRadioButtons(session, "outputScriptTo", selected = settings()$outputScriptTo)
            updateTextInput(session, "scriptFile", value = settings()$scriptFile)
            updateCheckboxInput(session, "createRStudioProj", value = settings()$createRStudioProj)
            updateRadioButtons(session, "ionization", selected = settings()$ionization)
        })
        
        list(
            valid = reactive({
                if (!nzchar(input$destination))
                    list(title = "Invalid destination", msg = "Please select a destination path!")
                else if (input$outputScriptTo != "curFile" && !nzchar(input$scriptFile))
                    list(title = "No script file", msg = "Please select a destination script file!")
                else
                    TRUE
            }),
            settings = reactive(list(
                destination = input$destination,
                outputScriptTo = input$outputScriptTo,
                scriptFile = input$scriptFile,
                createRStudioProj = input$createRStudioProj,
                ionization = input$ionization
            )),
            loadParams = reactive(input$loadParams),
            saveParams = reactive(input$saveParams)
        )
    })
}
