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
                            multiple = FALSE, width = "100%"),
                conditionalPanel(
                    condition = "input.componAlgo != \"nontarget\" && input.componAlgo != \"\"",
                    ns = ns,
                    checkboxInput(ns("selectIons"), "Select feature adduct ions")
                )
            ),
            fillRow(
                height = 90,
                fillCol(
                    flex = c(1, NA),
                    selectInput(ns("formulasAlgo"), "Formula generation",
                                c("None" = "", "GenForm", "SIRIUS", "Bruker DataAnalysis" = "Bruker"),
                                multiple = FALSE, width = "95%"),
                    textNote("DataAnalysis only works with features from DataAnalysis")
                ),
                selectInput(ns("compoundsAlgo"), "Compound identification",
                            c("None" = "", "SIRIUS+CSI:FingerID" = "SIRIUS", "MetFrag", "Library"),
                            multiple = FALSE, width = "100%")
            ),
            conditionalPanel(
                condition = "input.formulasAlgo != \"\" || input.compoundsAlgo != \"\"",
                ns = ns,
                fillRow(
                    height = 110,
                    fillCol(
                        selectInput(ns("peakListGen"), "Peak list generator",
                                    c("mzR", "Bruker DataAnalysis" = "Bruker"),
                                    multiple = FALSE, width = "95%"),
                        checkboxInput(ns("DIA"), "Data Independent Acquisition (DIA) MS/MS")
                    ),
                    conditionalPanel(
                        condition = "input.peakListGen == \"mzR\" && !input.DIA",
                        ns = ns,
                        numericInput(ns("precursorMzWindow"), "MS/MS precursor m/z search window", 4, width = "100%"),
                    )
                )
            ),
            conditionalPanel(
                condition = "input.compoundsAlgo == \"Library\"",
                ns = ns,
                fillRow(
                    height = 90,
                    fillCol(
                        width = "95%",
                        selectInput(ns("MSLibraryFormat"), "Library format",
                                    c("MSP" = "msp", "MoNA JSON" = "json"), multiple = FALSE)
                    ),
                    fillCol(
                        width = "95%",
                        fileSelect(ns("MSLibraryPath"), ns("MSLibraryPathButton"), "MS library path")
                    )
                )
            ),
            conditionalPanel(
                condition = "output.hasSuspects && (input.formulasAlgo != \"\" || input.compoundsAlgo != \"\")",
                ns = ns,
                fillRow(
                    height = 90,
                    fillCol(
                        strong("Suspect annotation"),
                        checkboxInput(ns("annotateSus"), "Annotate suspects", TRUE, width = "100%"),
                        conditionalPanel(
                            condition = "input.annotateSus",
                            ns = ns,
                            checkboxInput(ns("genIDLevelFile"), "Generate template file with configurable identification levels",
                                          width = "100%")
                        ),
                        textNote("Suspect annotation is currently only optimized for GenForm/MetFrag")
                    )
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
            updateSelectInput(session, "peakListGen", selected = settings()$peakListGen)
            updateCheckboxInput(session, "DIA", value = settings()$DIA)
            updateNumericInput(session, "precursorMzWindow", value = settings()$precursorMzWindow)
            updateSelectInput(session, "MSLibraryFormat", selected = settings()$MSLibraryFormat)
            updateTextInput(session, "MSLibraryPath", value = settings()$MSLibraryPath)
            updateCheckboxInput(session, "annotateSus", value = settings()$annotateSus)
            updateCheckboxInput(session, "genIDLevelFile", value = settings()$genIDLevelFile)
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
                peakListGen = input$peakListGen,
                DIA = input$DIA,
                precursorMzWindow = input$precursorMzWindow,
                MSLibraryFormat = input$MSLibraryFormat,
                MSLibraryPath = input$MSLibraryPath,
                annotateSus = input$annotateSus,
                genIDLevelFile = input$genIDLevelFile
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
        peakListGen = "mzR",
        DIA = FALSE,
        precursorMzWindow = 4,
        MSLibraryFormat = "msp",
        MSLibraryPath = "",
        annotateSus = TRUE,
        genIDLevelFile = TRUE
    ))
}

upgradeAnnotationSettings <- function(settings)
{
    # NOTE: this updates from first file version
    ret <- modifyList(defaultAnnotationSettings(),
                      settings[c("selectIons", "peakListGen", "DIA", "precursorMzWindow", "MSLibraryFormat",
                                 "MSLibraryPath", "annotateSus", "genIDLevelFile")])
    ret$componAlgo <- settings$components
    ret$formulasAlgo <- ret$formulaGen
    ret$compoundsAlgo <- ret$compIdent
    return(ret)
}
