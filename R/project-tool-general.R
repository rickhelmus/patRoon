newProjectGeneralUI <- function(id)
{
    ns <- NS(id)
    tagList(
        miniUI::miniContentPanel(
            fillCol(
                flex = c(1, NA, 1, NA, 1),
                fillRow(
                    fillCol(
                        height = 105,
                        width = "90%",
                        fileSelect(ns("destination"), ns("destinationButton"), "Project destination"),
                        checkboxInput(ns("createRStudioProj"), "Create (and open) RStudio project")
                    ),
                    fillCol(
                        height = 100,
                        textInput(ns("scriptFile"), label = "Script file", width = "90%"),
                        textNote("Make empty to add code to currently opened file.")
                    )
                ),
                hr(),
                fillRow(
                    flex = c(1, 2),
                    width = "90%",
                    radioButtons(ns("IMS"), "Mobility asignment", c("none", "from features" = "direct",
                                                                    "post features" = "post")),
                    conditionalPanel(
                        condition = "input.IMS != \"none\"",
                        ns = ns,
                        fillRow(
                            height = 130,
                            fillCol(
                                radioButtons(ns("IMSLimits"), "IMS configuration",
                                             c("Bruker" = "bruker", "Agilent" = "agilent")),
                                textNote("Configures the limits.yml file.")
                            ),
                            fillCol(
                                selectInput(ns("CCSMethod"), "CCS calculation method",
                                            choices = c("none",
                                                        "Bruker (TIMS)" = "bruker",
                                                        "Mason-Schamp 1/k (TIMS)" = "mason-schamp_1/k",
                                                        "Agilent (DTIMS)" = "agilent",
                                                        "Mason-Schamp" = "mason-schamp_k"),
                                            selected = "none"),
                                conditionalPanel(
                                    condition = "input.CCSMethod == 'agilent'",
                                    ns = ns,
                                    fileSelect(ns("CCSCalibrant"), ns("CCSCalibrantButton"), NULL,
                                               placeholder = "Calibrant .d file")
                                )
                            )
                        )
                    )
                ),
                hr(),
                radioButtons(ns("ionization"), "Ionization", c("positive", "negative", "both (sets)" = "both"))
            )
        ),
        miniUI::miniButtonBlock(
            actionButton(ns("loadParams"), "Load parameters", icon("upload")),
            actionButton(ns("saveParams"), "Save parameters", icon("download")),
        )
    )
}

newProjectGeneralServer <- function(id, settings)
{
    moduleServer(id, function(input, output, session)
    {
        doObserveSelDir(input, session, "destination", "destinationButton")
        doObserveSelDir(input, session, "CCSCalibrant", "CCSCalibrantButton")
        
        observeEvent(settings(), {
            updateTextInput(session, "destination", value = settings()$destination)
            updateRadioButtons(session, "outputScriptTo", selected = settings()$outputScriptTo)
            updateTextInput(session, "scriptFile", value = settings()$scriptFile)
            updateCheckboxInput(session, "createRStudioProj", value = settings()$createRStudioProj)
            updateRadioButtons(session, "IMS", selected = settings()$IMS)
            updateRadioButtons(session, "IMSLimits", selected = settings()$IMSLimits)
            updateSelectInput(session, "CCSMethod", selected = settings()$CCSMethod)
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
                IMS = input$IMS,
                IMSLimits = input$IMSLimits,
                CCSMethod = input$CCSMethod,
                ionization = input$ionization
            )),
            CCSCalibrant = reactive(input$CCSCalibrant),
            loadParams = reactive(input$loadParams),
            saveParams = reactive(input$saveParams)
        )
    })
}

defaultGeneralSettings <- function(destPath)
{
    return(list(
        destination = if (is.null(destPath)) "~/" else destPath,
        outputScriptTo = "newFile",
        scriptFile = "process.R",
        createRStudioProj = TRUE,
        IMS = "none",
        IMSLimits = "bruker",
        CCSMethod = "none",
        ionization = "positive"
    ))
}

upgradeGeneralSettings <- function(settings, destPath)
{
    # NOTE: this updates from first file version
    return(modifyList(defaultGeneralSettings(destPath),
                      settings[c("outputScriptTo", "scriptFile", "createRStudioProj", "ionization")]))
}
