getIDConfUIOpts <- function()
{
    return(c(
        "Suspects" = "suspects",
        "Formula candidates" = "formulas",
        "Compound candidates" = "compounds"
    ))
}

newProjectAnnotationUI <- function(id)
{
    ns <- NS(id)
    
    miniUI::miniContentPanel(
        fillCol(
            flex = NA,
            fillCol(
                height = 55,
                selectInput(ns("componAlgo"), "Component generation",
                            c("None" = "", "RAMClustR", "CAMERA", "OpenMS", "CliqueMS", "nontarget"),
                            multiple = FALSE, width = "97.5%")
            ),
            conditionalPanel(
                condition = "input.componAlgo != \"nontarget\" && input.componAlgo != \"\"",
                ns = ns,
                fillRow(
                    height = 25,
                    checkboxInput(ns("selectIons"), "Select feature adduct ions")
                )
            ),
            br(),
            fillCol(
                height = 55,
                selectInput(ns("formulasAlgo"), "Formula generation",
                            c("None" = "", "GenForm", "SIRIUS", "Bruker DataAnalysis" = "Bruker"),
                            multiple = FALSE, width = "97.5%")
            ),
            conditionalPanel(
                condition = "input.formulasAlgo == \"Bruker\"",
                ns = ns,
                fillCol(
                    height = 25,
                    textNote("DataAnalysis only works with features from DataAnalysis")    
                )
            ),
            br(),
            fillCol(
                height = 60,
                selectInput(ns("compoundsAlgo"), "Compound identification",
                            c("None" = "", "SIRIUS+CSI:FingerID" = "SIRIUS", "MetFrag", "Library"),
                            multiple = FALSE, width = "97.5%")
            ),
            conditionalPanel(
                condition = "input.compoundsAlgo == \"Library\"",
                ns = ns,
                fillRow(
                    height = 60,
                    fillCol(
                        width = "95%",
                        selectInput(ns("MSLibraryFormat"), "Library format",
                                    c("MSP" = "msp", "MoNA JSON" = "json"), multiple = FALSE, width = "100%")
                    ),
                    fillCol(
                        width = "95%",
                        fileSelect(ns("MSLibraryPath"), ns("MSLibraryPathButton"), "MS library path")
                    )
                )
            ),
            hr(),
            fillRow(
                flex = c(2, 1),
                fillCol(
                    flex = NA,
                    fillCol(
                        height = 90,
                        checkboxGroupInput(ns("estIDConf"), "Identification confidence estimation (if data is present)",
                                           getIDConfUIOpts(), getIDConfUIOpts(), width = "100%"),
                    ),
                    fillCol(
                        height = 25,
                        textNote("ID confidence estimation is currently only optimized for GenForm/MetFrag")
                    )
                ),
                conditionalPanel(
                    condition = "output.IMSMode != 'none'",
                    ns = ns,
                    fillRow(
                        height = 75,
                        radioButtons(ns("compCCSPred"), "CCS prediction", getCCSPredSelections())
                    )
                )
            )
        )
    )
}

newProjectAnnotationServer <- function(id, hasSuspects, IMSMode, settings)
{
    moduleServer(id, function(input, output, session)
    {
        observeEvent(settings(), {
            updateSelectInput(session, "componAlgo", selected = settings()$componAlgo)
            updateCheckboxInput(session, "selectIons", value = settings()$selectIons)
            updateSelectInput(session, "formulasAlgo", selected = settings()$formulasAlgo)
            updateSelectInput(session, "compoundsAlgo", selected = settings()$compoundsAlgo)
            updateSelectInput(session, "MSLibraryFormat", selected = settings()$MSLibraryFormat)
            updateTextInput(session, "MSLibraryPath", value = settings()$MSLibraryPath)
            updateCheckboxGroupInput(session, "estIDConf", selected = settings()$estIDConf)
            updateRadioButtons(session, "compCCSPred", selected = settings()$compCCSPred)
        })
        
        observeEvent(input$MSLibraryPathButton, {
            MSFile <- rstudioapi::selectFile(path = "~/")
            if (!is.null(MSFile))
                updateTextInput(session, "MSLibraryPath", value = MSFile)
        })

        output <- exportShinyOutputVal(output, "hasSuspects", hasSuspects)
        output <- exportShinyOutputVal(output, "IMSMode", IMSMode)
        
        list(
            valid = reactive({
                if (input$compoundsAlgo == "Library" && !nzchar(input$MSLibraryPath))
                    list(title = "No library path", msg = "Please select a library path!")
                else
                    TRUE
            }),
            settings = reactive(list(
                componAlgo = input$componAlgo,
                selectIons = input$selectIons,
                formulasAlgo = input$formulasAlgo,
                compoundsAlgo = input$compoundsAlgo,
                MSLibraryFormat = input$MSLibraryFormat,
                MSLibraryPath = input$MSLibraryPath,
                estIDConf = input$estIDConf,
                compCCSPred = input$compCCSPred
            ))
        )
    })
}

defaultAnnotationSettings <- function()
{
    return(list(
        componAlgo = "",
        selectIons = TRUE,
        formulasAlgo = "",
        compoundsAlgo = "",
        MSLibraryFormat = "msp",
        MSLibraryPath = "",
        estIDConf = unname(getIDConfUIOpts()),
        compCCSPred = "none"
    ))
}

upgradeAnnotationSettings <- function(settings)
{
    # NOTE: this updates from first file version
    # NOTE: some settings are ignored (eg MSPL algos, DIA)
    ret <- modifyList(defaultAnnotationSettings(), settings[c("selectIons", "MSLibraryFormat", "MSLibraryPath")])
    ret$componAlgo <- settings$components
    ret$formulasAlgo <- ret$formulaGen
    ret$compoundsAlgo <- ret$compIdent
    if (!settings$annotateSus)
        ret$estIDConf <- setdiff(ret$estIDConf, "suspects")
    return(ret)
}
