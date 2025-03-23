newProjectGeneralUI <- function(id, destPath)
{
    ns <- NS(id)
    tagList(
        miniUI::miniContentPanel(
            fillCol(
                flex = NA,
                fileSelect(ns("destinationPath"), ns("destinationPathButton"), "Project destination",
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
        doObserveSelDir(input, session, "destinationPath", "destinationPathButton")
        
        observeEvent(settings(), {
            updateRadioButtons(session, "ionization", selected = settings()$ionization)
        })
        
        list(
            destination = reactive(input$destinationPath),
            outputScriptTo = reactive(input$outputScriptTo),
            scriptFile = reactive(input$scriptFile),
            createRStudioProj = reactive(input$createRStudioProj),
            ionization = reactive(input$ionization),
            loadParams = reactive(input$loadParams),
            saveParams = reactive(input$saveParams)
        )
    })
}
