# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
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
                "Pre-treatment", icon = icon("file-export"),
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

# NOTE: the read/write functions below are also used by tests
readProjectSettings <- function(file, destPath)
{
    settings <- readYAML(file)
    
    if (is.null(settings[["version"]]))
        settings$version <- 1L # first file version was unversioned
    if (settings$version < 2L) 
    {
        settings <- list(general = upgradeGeneralSettings(settings, destPath),
                         analyses = upgradeAnalysesSettings(settings),
                         preTreatment = upgradePreTreatSettings(settings),
                         features = upgradeFeaturesSettings(settings),
                         annotations = upgradeAnnotationsSettings(settings),
                         TPs = upgradeTPsSettings(settings),
                         report = upgradeReportSettings(settings))
    }
    
    # HACK: make sure format/content matches regular settings. If this gets more complex, separate this into
    # functions.
    
    settings <- settings[names(settings) != "version"] # added in when saving YML
    
    # steps is a DT, which isn't transferred well to/from YML
    settings$preTreatment$steps <- data.table(algorithm = as.character(settings$preTreatment$steps$algorithm),
                                              from = as.character(settings$preTreatment$steps$from),
                                              to = as.character(settings$preTreatment$steps$to))
    
    # handle empty vector to list converion
    settings$annotations$estIDConf <- as.character(settings$annotations$estIDConf)
    
    return(settings)
}

writeProjectSettings <- function(settings, file)
{
    settings$version <- 2L
    writeYAML(settings, file)
}

newProjectServer <- function(destPath)
{
    checkExistingScript <- function(settings)
    {
        if (!nzchar(settings$scriptFile))
            return(TRUE)
        p <- file.path(settings$destination, settings$scriptFile)
        return(!file.exists(p) || rstudioapi::showQuestion("Script file already exists",
                                                           sprintf("Script file already exists: '%s'.\nOverwrite?", p),
                                                           "Yes", "No"))
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
                                         annotations = defaultAnnotationSettings(),
                                         TPs = defaultTPSettings(),
                                         report = defaultReportSettings())
        data <- list()
        
        data$general <- newProjectGeneralServer("general", reactive(loadedSettings$general))
        ionization <- reactive(data$general$settings()$ionization)
        IMSMode <- reactive(data$general$settings()$IMS$mode)
        data$analyses <- newProjectAnalysesServer("analyses", ionization, reactive(loadedSettings$analyses))
        data$preTreatment <- newProjectPreTreatServer("pretreat", ionization, reactive(loadedSettings$preTreatment))
        data$features <- newProjectFeaturesServer("features", ionization, IMSMode, reactive(loadedSettings$features))
        hasSusp <- reactive({
            data$features$settings()$exSuspList || (ionization() != "both" && nzchar(data$features$settings()$suspects$single) ||
                                                        (ionization() == "both" && nzchar(data$features$settings()$suspects$positive)))
        })
        data$annotations <- newProjectAnnotationServer("annotation", hasSusp, IMSMode, reactive(loadedSettings$annotations))
        data$TPs <- newProjectTPServer("tp", hasSusp, reactive(data$annotations$settings()$formulasAlgo),
                                       reactive(data$annotations$settings()$compoundsAlgo), reactive(loadedSettings$TPs))
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
                doCreateProject(data$general$CCSCalibrant(), data$analyses$anaInfoTabs(), getSettings())
                stopApp(TRUE)
            }
        })
        
        observeEvent(data$general$loadParams(), {
            sl <- rstudioapi::selectFile("Select parameter file", filter = "yml files (*.yml)")
            if (!is.null(sl))
            {
                ymlSettings <- readProjectSettings(sl, destPath)
                for (s in names(ymlSettings))
                    loadedSettings[[s]] <- ymlSettings[[s]]
            }
        })
        
        observeEvent(data$general$saveParams(), {
            sl <- rstudioapi::selectFile("Select parameter file", filter = "yml files (*.yml)", existing = FALSE)
            if (!is.null(sl))
            {
                out <- getSettings()
                writeProjectSettings(out, sl)
            }
        })
    }
    
    return(server)
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
    
    runGadget(getNewProjectUI(), newProjectServer(destPath),
              viewer = dialogViewer("Create new project", width = 800, height = 600))
    # runGadget(getNewProjectUI(), server, viewer = paneViewer())
}
