# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
NULL

getNewProjectUI <- function()
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
                newProjectGeneralUI("general")
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
    
    checkExistingScript <- function(settings)
    {
        p <- file.path(settings$destination, settings$scriptFile)
        return(settings$outputScriptTo == "curFile" || !file.exists(p) ||
                   rstudioapi::showQuestion("Script file already exists",
                                            sprintf("Script file already exists: '%s'.\nOverwrite?", p), "Yes", "No"))
    }
    
    checkExistingAnaInfo <- function(settingsGen, settingsAna)
    {
        if (settingsAna$generateAnaInfo == "table")
        {
            n <- paste0("analysisTableFile", settingsAna$analysisTableFileType)
            checkAnas <- if (settingsGen$ionization != "both")
                input[[n]]
            else
                input[paste0(n, c("Pos", "Neg"))]
            for (f in checkAnas)
            {
                p <- file.path(settingsGen$destination, f)
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
    
    server <- function(input, output, session)
    {
        loadedSettings <- reactiveValues(general = defaultGeneralSettings(destPath),
                                         analyses = defaultAnalysesSettings(),
                                         preTreatment = defaultPreTreatSettings(),
                                         features = defaultFeaturesSettings(),
                                         annotation = defaultAnnotationSettings(),
                                         TP = defaultTPSettings(),
                                         report = defaultReportSettings())
        data <- list()
        
        data$general <- newProjectGeneralServer("general", reactive(loadedSettings$general))
        ionization <- reactive(data$general$settings()$ionization)
        data$analyses <- newProjectAnalysesServer("analyses", ionization, reactive(loadedSettings$analyses))
        data$preTreatment <- newProjectPreTreatServer("pretreat", ionization, reactive(loadedSettings$preTreatment))
        data$features <- newProjectFeaturesServer("features", ionization, reactive(loadedSettings$features))
        hasSusp <- reactive({
            data$features$settings()$exSuspList || (ionization() != "both" && nzchar(data$features$settings()$suspects$single) ||
                (ionization() == "both" && nzchar(data$features$settings()$suspects$positive)))
        })
        data$annotations <- newProjectAnnotationServer("annotation", hasSusp, reactive(loadedSettings$annotation))
        data$TPs <- newProjectTPServer("tp", hasSusp, reactive(data$annotations$settings()$formulasAlgo),
                                       reactive(data$annotations$settings()$compoundsAlgo), reactive(loadedSettings$TP))
        data$report <- newProjectReportServer("report", reactive(loadedSettings$report))
        
        getSettings <- function() sapply(data, \(d) d$settings(), simplify = FALSE)

        validate <- function()
        {
            for (d in data)
            {
                v <- d$valid()
                if (!isTRUE(v))
                {
                    rstudioapi::showDialog(v$title, v$msg, "")
                    return(FALSE)
                }
            }
            return(TRUE)
        }
        
        observeEvent(input$create, {
            if (!validate())
            {}
            else if (!checkExistingScript(data$general$settings()))
            {}
            else if (!checkExistingAnaInfo(data$general$settings(), data$analyses$settings()))
            {}
            else
            {
                doCreateProject(data$analyses$anaInfoTabs(), getSettings())
                stopApp(TRUE)
            }
        })

        observeEvent(data$general$loadParams(), {
            sl <- rstudioapi::selectFile("Select parameter file", filter = "yml files (*.yml)")
            if (!is.null(sl))
            {
                ymlSettings <- readYAML(sl)
                if (is.null(ymlSettings[["version"]]))
                    ymlSettings$version <- 1L # first file version was unversioned
                if (ymlSettings$version < 2L) 
                {
                    ymlSettings <- list(general = upgradeGeneralSettings(ymlSettings, destPath),
                                        analyses = upgradeAnalysesSettings(ymlSettings),
                                        preTreatment = upgradePreTreatSettings(ymlSettings),
                                        feature = upgradeFeatureSettings(ymlSettings),
                                        annotation = upgradeAnnotationSettings(ymlSettings),
                                        TP = upgradeTPSettings(ymlSettings),
                                        report = upgradeReportSettings(ymlSettings))
                }
                for (s in names(ymlSettings))
                    loadedSettings[[s]] <- ymlSettings[[s]]
            }
        })
        
        observeEvent(data$general$saveParams(), {
            sl <- rstudioapi::selectFile("Select parameter file", filter = "yml files (*.yml)", existing = FALSE)
            if (!is.null(sl))
            {
                out <- getSettings()
                out$version <- 2L
                writeYAML(out, sl)
            }
        })
    }

    runGadget(getNewProjectUI(), server, viewer = dialogViewer("Create new project", width = 800, height = 600))
    # runGadget(getNewProjectUI(), server, viewer = paneViewer())
}
