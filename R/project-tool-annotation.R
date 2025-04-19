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
                height = 100,
                selectInput(ns("componAlgo"), "Component generation",
                            c("None" = "", "RAMClustR", "CAMERA", "OpenMS", "CliqueMS", "nontarget"),
                            multiple = FALSE, width = "97.5%"),
                conditionalPanel(
                    condition = "input.componAlgo != \"nontarget\" && input.componAlgo != \"\"",
                    ns = ns,
                    checkboxInput(ns("selectIons"), "Select feature adduct ions")
                )
            ),
            fillRow(
                height = 120,
                fillCol(
                    flex = NA,
                    fillCol(
                        height = 50,
                        selectInput(ns("formulasAlgo"), "Formula generation",
                                    c("None" = "", "GenForm", "SIRIUS", "Bruker DataAnalysis" = "Bruker"),
                                    multiple = FALSE, width = "95%")
                    ),
                    fillCol(
                        height = 25,
                        textNote("DataAnalysis only works with features from DataAnalysis")    
                    )
                ),
                selectInput(ns("compoundsAlgo"), "Compound identification",
                            c("None" = "", "SIRIUS+CSI:FingerID" = "SIRIUS", "MetFrag", "Library"),
                            multiple = FALSE, width = "95%")
            ),
            conditionalPanel(
                condition = "input.compoundsAlgo == \"Library\"",
                ns = ns,
                fillRow(
                    height = 90,
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
            )
        )
    )
}

newProjectAnnotationServer <- function(id, hasSuspects, settings)
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
        })
        
        observeEvent(input$MSLibraryPathButton, {
            MSFile <- rstudioapi::selectFile(path = "~/")
            if (!is.null(MSFile))
                updateTextInput(session, "MSLibraryPath", value = MSFile)
        })

        output <- exportShinyOutputVal(output, "hasSuspects", hasSuspects)
        
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
                estIDConf = input$estIDConf
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
        estIDConf = unname(getIDConfUIOpts())
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
