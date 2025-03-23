# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
NULL

getNewProjectUI <- function(destPath)
{
    miniUI::miniPage(
        shinyjs::useShinyjs(),
        
        # based on https://stackoverflow.com/a/63882648
        tags$style(
            type = 'text/css',
            '.modal-dialog { width: fit-content !important; height: fit-content !important; }'
        ),
        
        tags$script(htmlwidgets::JS('
                Shiny.addCustomMessageHandler("selectHOTRows", function(data)
                {
                    // get rhot instance: https://github.com/jrowen/rhandsontable/issues/97
                    var ht = HTMLWidgets.getInstance(window[data.tab]).hot;
                    ht.selectRows(data.range[0] - 1, data.range[1] - 1);
                });
            ')
        ),
        
        miniUI::gadgetTitleBar("Create project tool", right = miniUI::miniTitleBarButton("create", "Create", TRUE)),

        miniUI::miniTabstripPanel(
            miniUI::miniTabPanel(
                "General", icon = icon("save"),
                newProjectGeneralUI("general", destPath)
            ),
            miniUI::miniTabPanel(
                "Analyses", icon = icon("folder-open"),
                newProjectAnalysesUI("analyses")
            ),
            miniUI::miniTabPanel(
                "Data pre-treatment", icon = icon("upload"),
                newProjectPreTreatUI("pretreat")
            ),
            miniUI::miniTabPanel(
                "Features", icon = icon("chart-area"),
                newProjectFeaturesUI("features")
            ),
            miniUI::miniTabPanel(
                "Annotation", icon = icon("chart-bar"),
                newProjectAnnotationUI("annotation")
            ),
            miniUI::miniTabPanel(
                "TP screening", icon = icon("react"),
                newProjectTPUI("tp")
            ),
            miniUI::miniTabPanel(
                "Reporting", icon = icon("file-medical-alt"),
                newProjectReportUI("report")
            )
        )
    )
}

getNewProjectWidgetTypes <- function()
{
    list(
        outputScriptTo = "radio",
        scriptFile = "text",
        createRStudioProj = "check",
        ionization = "radio",
        analysisTableFile = "text",
        analysisTableFilePos = "text",
        analysisTableFileNeg = "text",
        convAlgo = "select",
        convFrom = "select",
        convTo = "select",
        peakPicking = "check",
        peakPickingVendor = "check",
        DAMethod = "text", # UNDONE: remove
        DAMethodPos = "text", # UNDONE: remove
        DAMethodNeg = "text", # UNDONE: remove
        doBrukerCalib = "check", # UNDONE: remove
        featFinder = "select",
        featGrouper = "select",
        suspectList = "text",
        suspectListPos = "text",
        suspectListNeg = "text",
        exSuspList = "check",
        preIntThr = "numeric",
        intThr = "numeric",
        repAbundance = "numeric",
        maxRepRSD = "numeric",
        blankThr = "numeric",
        "retention-min" = "numeric",
        "retention-max" = "numeric",
        "mz-min" = "numeric",
        "mz-max" = "numeric",
        removeBlanks = "check",
        featNorm = "select",
        groupNorm = "check",
        ISTDList = "text",
        ISTDListPos = "text",
        ISTDListNeg = "text",
        components = "select",
        selectIons = "check",
        formulaGen = "select",
        compIdent = "select",
        peakListGen = "select",
        DIA = "check",
        precursorMzWindow = "numeric",
        MSLibraryFormat = "select",
        MSLibraryPath = "text",
        annotateSus = "check",
        genIDLevelFile = "check",
        doTPs = "check",
        TPGen = "select",
        TPGenInput = "select",
        TPSuspectList = "text",
        TPDoMFDB = "check",
        reportGen = "checkGroup",
        reportLegacy = "checkGroup"
    )
}

loadNewProjectParams <- function(file, input, session)
{
    wtypes <- getNewProjectWidgetTypes()
    values <- readYAML(file)
    for (param in names(values))
    {
        if (wtypes[[param]] == "radio")
            updateRadioButtons(session, param, selected = values[[param]])
        else if (wtypes[[param]] == "text")
            updateTextInput(session, param, value = values[[param]])
        else if (wtypes[[param]] == "check")
            updateCheckboxInput(session, param, value = values[[param]])
        else if (wtypes[[param]] == "checkGroup")
            updateCheckboxGroupInput(session, param, selected = values[[param]])
        else if (wtypes[[param]] == "select")
            updateSelectInput(session, param, selected = values[[param]])
        else if (wtypes[[param]] == "numeric")
            updateNumericInput(session, param, value = values[[param]])
    }
}

saveNewProjectParams <- function(file, input)
{
    values <- getNewProjectWidgetTypes()
    values <- isolate(reactiveValuesToList(input))[names(values)]
    writeYAML(values, file)
}

#' Easily create new \pkg{patRoon} projects
#' 
#' The \code{newProject} function is used to quickly generate a processing R script. This tool allows the user to
#' quickly select the targeted analyses, workflow steps and configuring some of their common parameters. This function
#' requires to be run within a \href{https://www.rstudio.com/}{RStudio} session. The resulting script is either added to
#' the current open file or to a new file. The \link[=analysis-information]{analysis information} will be written to a
#' \file{.csv} file so that it can easily be modified afterwards.
#'
#' @param destPath Set destination path value to this value (useful for debugging). Set to \code{NULL} for a default
#'   value.
#'
#' @export
newProject <- function(destPath = NULL)
{
    rstudioapi::verifyAvailable()

    # UNDONE: warning/message about empty groups
    
    server <- function(input, output, session)
    {
        rValues <- reactiveValues(generalSettings = list())

        general <- newProjectGeneralServer("general", reactive(rValues$generalSettings))
        analyses <- newProjectAnalysesServer("analyses", general$ionization, reactive(rValues$analysesSettings))
        MSConversion <- newProjectPreTreatServer("pretreat", general$ionization, reactive(rValues$preTreatSettings))
        features <- newProjectFeaturesServer("features", general$ionization, reactive(rValues$featureSettings))
        hasSusp <- reactive({
            features$exSuspList() || (general$ionization() != "both" && nzchar(features$suspectList())) ||
                (general$ionization() == "both" && nzchar(features$suspectListPos()))
        })
        annotations <- newProjectAnnotationServer("annotation", hasSusp, reactive(rValues$annotationSettings))
        TPs <- newProjectTPServer("tp", hasSusp, reactive(rValues$TPSettings))
        report <- newProjectReportServer("report", reactive(rValues$reportSettings))

        verifyAnalysesOK <- function()
        {
            verifyAny <- function(pol)
            {
                if (nrow(anaInfoTabs[[pol]]) == 0 && input$ionization %in% c(pol, "both"))
                {
                    rstudioapi::showDialog("No analyses", paste0("Please select some analyses for ", pol, " mode!"), "")
                    return(FALSE)
                }
            }
            verifyAny("positive"); verifyAny("negative")

            if (input$generateAnaInfo == "table")
            {
                n <- paste0("analysisTableFile", input$analysisTableFileType)
                checkAnas <- if (input$ionization != "both")
                    input[[n]]
                else
                    input[paste0(n, c("Pos", "Neg"))]
                for (f in checkAnas)
                {
                    p <- file.path(input$destinationPath, f)
                    if (file.exists(p))
                    {
                        ov <- rstudioapi::showQuestion("Analysis table file already exists",
                                                       sprintf("Analysis table file already exists: '%s'.\nOverwrite?", p),
                                                       "Yes", "No")
                        if (!ov)
                            return(FALSE)
                    }
                }
            }
            
            return(TRUE)
        }
        
        observeEvent(input$create, {
            if (!nzchar(input$destinationPath))
                rstudioapi::showDialog("Invalid destination", "Please select a destination path!", "")
            else if (input$outputScriptTo != "curFile" && !nzchar(input$scriptFile))
                rstudioapi::showDialog("No script file", "Please select a destination script file!", "")
            else if (input$generateAnaInfo == "table" && !verifyAnalysesOK())
            {}
            else if (input$outputScriptTo != "curFile" && file.exists(file.path(input$destinationPath, input$scriptFile)) &&
                     !rstudioapi::showQuestion("Script file already exists",
                                               sprintf("Script file already exists: '%s'.\nOverwrite?",
                                                       file.path(input$destinationPath, input$scriptFile)),
                                               "Yes", "No"))
            {}
            else if (input$compIdent == "Library" && !nzchar(input$MSLibraryPath))
                rstudioapi::showDialog("No MS library", "Please select an MS library!", "")
            else if (input$doTPs && input$TPGen != "Logic" && input$TPGenInput == "suspects" &&
                     !nzchar(input$TPSuspectList))
                rstudioapi::showDialog("No parent suspect list", "Please select a parent suspect list!", "")
            else if ("legacy" %in% input$reportGen && length(input$reportLegacy) == 0)
                rstudioapi::showDialog("No legacy format", "Please select at least one legacy reporting format!", "")
            else
            {
                doCreateProject(input, anaInfoTabs, MSConvSettings$steps)
                stopApp(TRUE)
            }
        })

        observeEvent(general$loadParams(), {
            sl <- rstudioapi::selectFile("Select parameter file", filter = "yml files (*.yml)")
            if (!is.null(sl))
                loadNewProjectParams(sl, input, session)
        })
        
        observeEvent(input$saveParams, {
            sl <- rstudioapi::selectFile("Select parameter file", filter = "yml files (*.yml)", existing = FALSE)
            if (!is.null(sl))
                saveNewProjectParams(sl, input)
        })
    }

    runGadget(getNewProjectUI(destPath), server, viewer = dialogViewer("Create new project", width = 800, height = 600))
    # runGadget(getNewProjectUI(), server, viewer = paneViewer())
}
