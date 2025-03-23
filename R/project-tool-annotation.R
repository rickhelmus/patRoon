newProjectAnnotationUI <- function(id)
{
    ns <- NS(id)
    
    miniUI::miniContentPanel(
        fillCol(
            flex = NA,
            
            fillCol(
                height = 100,
                selectInput(ns("components"), "Component generation",
                            c("None" = "", "RAMClustR", "CAMERA", "OpenMS", "CliqueMS", "nontarget"),
                            multiple = FALSE, width = "100%"),
                conditionalPanel(
                    condition = "input.components != \"nontarget\" && input.components != \"\"",
                    ns = ns,
                    checkboxInput(ns("selectIons"), "Select feature adduct ions", value = TRUE)
                )
            ),
            fillRow(
                height = 90,
                fillCol(
                    flex = c(1, NA),
                    selectInput(ns("formulaGen"), "Formula generation",
                                c("None" = "", "GenForm", "SIRIUS", "Bruker DataAnalysis" = "Bruker"),
                                multiple = FALSE, width = "95%"),
                    textNote("DataAnalysis only works with features from DataAnalysis")
                ),
                selectInput(ns("compIdent"), "Compound identification",
                            c("None" = "", "SIRIUS+CSI:FingerID" = "SIRIUS", "MetFrag", "Library"),
                            multiple = FALSE, width = "100%")
            ),
            conditionalPanel(
                condition = "input.formulaGen != \"\" || input.compIdent != \"\"",
                ns = ns,
                fillRow(
                    height = 110,
                    fillCol(
                        selectInput(ns("peakListGen"), "Peak list generator",
                                    c("mzR", "Bruker DataAnalysis" = "Bruker"),
                                    "mzR", multiple = FALSE, width = "95%"),
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
                condition = "input.compIdent == \"Library\"",
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
                condition = "output.hasSuspects && (input.formulaGen != \"\" || input.compIdent != \"\")",
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
                                          TRUE, width = "100%")
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
        observeEvent(input$MSLibraryPathButton, {
            MSFile <- rstudioapi::selectFile(path = "~/")
            if (!is.null(MSFile))
                updateTextInput(session, "MSLibraryPath", value = MSFile)
        })
        
        output <- exportShinyOutputVal(output, "hasSuspects", hasSuspects)
    })
}
