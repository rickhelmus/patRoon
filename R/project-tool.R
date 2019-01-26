#' @include main.R
NULL

getScriptCode <- function(destination, generateAnaInfo, analysisTableFile, analyses, DAMethod, dataPretreatments,
                          featFinder, featGrouper, intThreshold, replThreshold,
                          blankThreshold, filterRepetitions, peakListGenerator, precursorMzWindow, formulaGenerator,
                          compIdentifier, componentGenerator, polarity, reportFormats)
{
    optionalCodeBlock <- function(e) if (e) "<<startCodeBlock>>" else "<<skipCodeBlock>>"
    endCodeBlock <- function() "<<endCodeBlock>>"
    optionalLine <- function(e) if (!e) "<<skipThisLine>>" else ""

    # Can't set lists in template() call (bug?)
    dataPretreatmentOpts = list(DAMethod = DAMethod, steps = dataPretreatments)
    featFinderOpts = list(algo = featFinder)
    featGrouperOpts = list(algo = featGrouper)
    filterFGroupsOpts = list(intThr = intThreshold, replThr = replThreshold,
                             blankThr = blankThreshold, filterRepetitions = filterRepetitions)
    peakListOpts = list(algo = peakListGenerator)
    formulaOpts = list(algo = formulaGenerator)
    identOpts = list(algo = compIdentifier)
    componentOpts = list(algo = componentGenerator)
    sldest <- gsub("\\", "/", destination, fixed = TRUE) # BUG: tmpl() seems to eat backslashes!? Convert to forward slashes for now
    template <- tmpl(readAllFile(system.file("templates", "main_script.R", package = "patRoon")),
                     destination = sldest, generateAnaInfo = generateAnaInfo, analysisTableFile = analysisTableFile,
                     analyses = analyses,
                     doMSPeakFind = formulaGenerator != "" || compIdentifier != "",
                     precursorMzWindow = precursorMzWindow,
                     polarity = polarity, reportFormats = reportFormats)

    ret <- template

    # remove blocks marked as '<<skipCodeBlock>>' till <<endCodeBlock>> + leading whitespace
    ret <- gsub("[\t ]*<<skipCodeBlock>>[\\s\\S]*?<<endCodeBlock>>[.\r]*\n?", "", ret, perl = TRUE)

    # remove remaining startCodeBlock and endCodeBlock markers
    ret <- gsub("[\t ]*(<<startCodeBlock>>|<<endCodeBlock>>)[.\r]*\n?", "", ret, perl = TRUE)

    # remove disabled optional lines
    ret <- gsub("([\n]*).*<<skipThisLine>>[.\r]*\n?", "\\1", ret, perl = TRUE)

    return(ret)
}

doCreateProject <- function(destination, scriptFile, createRStudioProj, generateAnaInfo, analyses, analysisTableFile,
                            DAMethod, dataPretreatments, featFinder, featGrouper, intThreshold,
                            replThreshold, blankThreshold, filterRepetitions, peakListGenerator, precursorMzWindow,
                            formulaGenerator, compIdentifier, componentGenerator, polarity, reportFormats)
{
    mkdirp(destination)

    analyses <- copy(analyses)
    analyses[, group := ifelse(!nzchar(group), analysis, group)]

    # Make analysis table
    if (generateAnaInfo == "table")
        write.csv(analyses[, c("path", "analysis", "group", "ref")], file.path(destination, analysisTableFile),
                  row.names = FALSE)

    code <- getScriptCode(destination, generateAnaInfo, analysisTableFile, analyses, DAMethod, dataPretreatments,
                          featFinder, featGrouper, intThreshold, replThreshold, blankThreshold,
                          filterRepetitions, peakListGenerator, precursorMzWindow, formulaGenerator,
                          compIdentifier, componentGenerator, polarity, reportFormats)
    if (is.null(scriptFile))
    {
        # insert at end of current document
        insertText(Inf, code, getSourceEditorContext()$id)
    }
    else
    {
        sp <- file.path(destination, scriptFile)
        writeChar(code, file.path(destination, scriptFile))

        if (createRStudioProj)
        {
            rstudioapi::initializeProject(destination)
            rstudioapi::openProject(destination)
        }
        else
            navigateToFile(sp)
    }
}

getNewProjectUI <- function()
{
    textNote <- function(txt) div(style = "margin: 8px 0 12px; font-size: small", txt)
    fileSelect <- function(idText, idButton, label, value = "")
    {
        fillRow(
            flex = c(1, NA),
            textInput(idText, label, value, width = "100%"),
            actionButton(idButton, "", icon("folder-open"), style = "margin: 25px 0 0 15px")
        )
    }

    miniPage(
        gadgetTitleBar("Create project tool", right = miniTitleBarButton("create", "Create", TRUE)),

        miniTabstripPanel(
            miniTabPanel("Destination", icon = icon("save"),
                         miniContentPanel(
                             fillCol(
                                 fileSelect("destinationPath", "projectDestButton", "Project destination", "~/"),
                                 fillRow(
                                     radioButtons("outputScriptTo", "Insert code into", c("New file" = "newFile",
                                                                                          "Current file" = "curFile")),
                                     conditionalPanel(
                                         condition = "input.outputScriptTo == \"newFile\"",
                                         textInput("scriptFile", "Script file", "process.R"),
                                         checkboxInput("createRStudioProj", "Create RStudio project", value = TRUE)
                                     )
                                 )
                             )
                         )
            ),

            miniTabPanel("Analyses", icon = icon("folder-open"),
                         miniContentPanel(
                             fillCol(
                                 flex = c(NA, 1),
                                 fillRow(
                                     height = 125,
                                     radioButtons("generateAnaInfo", "Generate analysis information",
                                                  c("None" = "none", "From new csv file" = "table",
                                                    "Load in script" = "script")),
                                     conditionalPanel(
                                         condition = "input.generateAnaInfo == \"table\"",
                                         textInput("analysisTableFile", "Analysis table output file", "analyses.csv")
                                     )
                                 ),
                                 conditionalPanel(
                                     condition = "input.generateAnaInfo != \"none\"",
                                     rHandsontableOutput("analysesHot")
                                 )
                             )
                         ),

                         conditionalPanel(
                             condition = "input.generateAnaInfo != \"none\"",
                             miniButtonBlock(
                                 actionButton("addAnalysesDir", "Add analyses from directory ..."),
                                 actionButton("addAnalysesCSV", "Add analyses from csv file ..."),
                                 actionButton("removeAnalyses", "Remove")
                             )
                         )
            ),
            miniTabPanel("Bruker data pretreatment", icon = icon("upload"),
                         miniContentPanel(
                             fillCol(
                                 flex = NA,
                                 fillCol(
                                     flex = c(1, NA),
                                     height = 90,
                                     fileSelect("DAMethod", "DAMethodButton", "Set DataAnalysis method to"),
                                     textNote("Leaving this blank will not set any method")
                                 ),
                                 fillCol(
                                     flex = c(1, NA),
                                     height = 90,
                                     selectInput("pretreat", "Pretreatment steps", c("Recalibration" = "recalibrate",
                                                                                     "Export mzML" = "expMzML",
                                                                                     "Export mzXML" = "expMzXML"),
                                                 multiple = TRUE, width = "100%"),
                                     textNote("enviPick needs mzXML, XCMS and OpenMS need mzML")
                                 )
                             )
                         )
            ),
            miniTabPanel("Workflow", icon = icon("code-fork"),
                         miniContentPanel(
                             fillCol(
                                 flex = NA,

                                 fillRow(
                                     height = 75,
                                     selectInput("featFinder", "Feature finder", c("OpenMS", "XCMS", "enviPick",
                                                                                   "Bruker DataAnalysis" = "Bruker"),
                                                 "OpenMS", FALSE, width = "95%"),
                                     selectInput("featGrouper", "Feature grouper", c("OpenMS", "XCMS"),
                                                 "OpenMS", FALSE, width = "100%")
                                 ),
                                 fillRow(
                                     height = 90,
                                     fillCol(
                                         flex = c(1, NA),
                                         selectInput("formulaGen", "Formula generation",
                                                     c("None" = "", "GenForm", "SIRIUS", "Bruker DataAnalysis" = "Bruker"),
                                                     multiple = FALSE, width = "95%"),
                                         textNote("DataAnalysis only works with features from DataAnalysis")
                                     ),
                                     selectInput("compIdent", "Compound identification",
                                                 c("None" = "", "SIRIUS+CSI:FingerID" = "SIRIUS", "MetFrag"),
                                                 multiple = FALSE, width = "100%")
                                 ),
                                 selectInput("components", "Component generation",
                                             c("None" = "", "RAMClustR", "CAMERA", "nontarget"),
                                             multiple = FALSE, width = "100%"),
                                 conditionalPanel(
                                     condition = "input.formulaGen != \"\" || input.compIdent != \"\"",
                                     fillRow(
                                         height = 90,
                                         selectInput("peakListGen", "Peak list generator",
                                                     c("mzR", "Bruker DataAnalysis" = "Bruker"),
                                                     "mzR", multiple = FALSE, width = "95%"),
                                         conditionalPanel(
                                             condition = "input.peakListGen == \"mzR\"",
                                             numericInput("precursorMzWindow", "MS/MS precursor m/z search window", 8, width = "100%"),
                                             textNote("The precursor m/z search window applied when finding MS/MS spectra. Set to zero for DIA experiments.")
                                         )
                                     )
                                 ),
                                 conditionalPanel(
                                     condition = "input.formulaGen != \"\" || input.compIdent == \"\" || input.components != \"\"",
                                     selectInput("polarity", "Polarity", c("positive", "negative"), "positive",
                                                 multiple = FALSE, width = "100%")
                                 ),
                                 checkboxGroupInput("report", "Report generation",
                                                    c("CSV (text tables)" = "CSV", "PDF (basic plots)" = "PDF",
                                                      "HTML (easy browsable plots, bit slower to PDF)" = "MD"),
                                                    c("CSV", "MD"))
                             )
                         )
            ),
            miniTabPanel("Miscellaneous", icon = icon("sliders"),
                         miniContentPanel(
                             fillCol(
                                 numericInput("intThr", "Intensity threshold", 1E4, 0, step = 1000, width = "100%"),
                                 numericInput("replThr", "Min. replicate abundance", 0.75, 0, 1.0, 0.1, width = "100%"),
                                 numericInput("blankThr", "Min. blank threshold", 5, 0, step = 1, width = "100%"),
                                 numericInput("filterRepetitions", "Repetitions", 2, 1, step = 1, width = "100%")
                             )
                         )
            )
        )
    )
}

#' @details The \code{newProject} function is used to quickly generate a
#'   processing R script. This tool allows the user to quickly select the
#'   targeted analyses, workflow steps and configuring some of their common
#'   parameters. This function requires to be run within a
#'   \href{https://www.rstudio.com/}{RStudio} session. The resulting script is
#'   either added to the current open file or to a new file. The
#'   \link[=analysis-information]{analysis information} will be written to a \file{.csv}
#'   file so that it can easily be modified afterwards.
#'
#' @rdname GUI-utils
#' @export
newProject <- function()
{
    # UNDONE: warning/message about empty groups

    hotOpts <- list(rowHeaderWidth = 40, readOnly = TRUE,
                    columnSorting = TRUE, sortIndicator = TRUE, selectCallback = TRUE,
                    currentRowClassName = "currentRow", stretchH = "all",
                    preventOverflow = "horizontal", multiSelect = TRUE,
                    outsideClickDeselects = FALSE,
                    contextMenu = FALSE, manualColumnResize = TRUE)

    server <- function(input, output, session)
    {
        rValues <- reactiveValues(analyses = data.table(analysis = character(0), type = character(0),
                                                        group = character(0), ref = character(0), path = character(0)))

        observeEvent(input$create, {
            if (input$destinationPath == "")
                showDialog("Invalid destination", "Please select a destination path!", "")
            else if (input$outputScriptTo != "curFile" && input$scriptFile == "")
                showDialog("No script file", "Please select a destination script file!", "")
            else if (input$generateAnaInfo != "none" && nrow(rValues$analyses) == 0)
                showDialog("No analyses selected", "Please select some analyses!", "")
            else if (input$generateAnaInfo == "table" && file.exists(file.path(input$destinationPath, input$analysisTableFile)) &&
                     !showQuestion("Analysis table file already exists",
                                   sprintf("Analysis table file already exists: '%s'.\nOverwrite?",
                                           file.path(input$destinationPath, input$analysisTableFile)),
                                   "Yes", "No"))
            {}
            else if (input$outputScriptTo != "curFile" && file.exists(file.path(input$destinationPath, input$scriptFile)) &&
                     !showQuestion("Script file already exists:",
                                   sprintf("Script file already exists: '%s'.\nOverwrite?",
                                           file.path(input$destinationPath, input$scriptFile)),
                                   "Yes", "No"))
            {}
            else
            {
                doCreateProject(input$destinationPath,
                                if (input$outputScriptTo == "curFile") NULL else input$scriptFile,
                                input$createRStudioProj, input$generateAnaInfo, rValues$analyses,
                                input$analysisTableFile, input$DAMethod, input$pretreat, input$featFinder,
                                input$featGrouper, input$intThr, input$replThr, input$blankThr, input$filterRepetitions,
                                input$peakListGen, input$precursorMzWindow, input$formulaGen, input$compIdent, input$components,
                                input$polarity, input$report)
                stopApp(TRUE)
            }
        })

        observeEvent(input$projectDestButton, {
            dest <- selectDirectory("Select destination directory", path = input$destinationPath)
            if (!is.null(dest))
                updateTextInput(session, "destinationPath", value = dest)
        })

        observeEvent(input$analysesHot, {
            # HACK: input$analysesHot$params$maxRows: make sure we don't have empty table as hot_to_r errors otherwise
            if (input$analysesHot$params$maxRows > 0)
            {
                df <- hot_to_r(input$analysesHot)
                rValues$analyses[, c("group", "ref") := .(df$group, df$ref)]
            }
        })

        observeEvent(input$addAnalysesDir, {
            anaDir <- selectDirectory(path = "~/")
            if (!is.null(anaDir))
            {
                files <- list.files(anaDir, "\\.(mzML|mzXML|d)$", full.names = TRUE)

                dt <- data.table(path = dirname(files), analysis = simplifyAnalysisNames(files),
                                  group = "", ref = "", ext = tools::file_ext(files))
                dt[ext == "d", ext := "bruker"]
                dt[, type := paste0(.SD$ext, collapse = ", "), by = .(path, analysis)]
                dt[, ext := NULL]
                dt <- unique(dt, by = c("analysis", "path"))
                setcolorder(dt, c("analysis", "type", "group", "ref", "path"))

                rValues$analyses <- rbind(rValues$analyses, dt)
            }
        })

        observeEvent(input$addAnalysesCSV, {
            csvFile <- selectFile(path = "~/", filter = "csv files (*.csv)")
            if (!is.null(csvFile))
            {
                csvTab <- tryCatch(fread(csvFile, select = c("path", "analysis", "group", "ref")),
                                   error = function(e) FALSE)
                if (is.logical(csvTab))
                    showDialog("Error", "Failed to open/parse selected csv file!", "")
                else
                {
                    exts <- mapply(csvTab$analysis, csvTab$path, FUN = function(ana, path)
                    {
                        fp <- file.path(path, ana)
                        ret <- c("mzML", "mzXML", "d")
                        ret <- ret[file.exists(paste0(fp, ".", ret))]
                        ret[ret == "d"] <- "bruker"
                        return(paste0(ret, collapse = ", "))
                    })

                    csvTab[, type := exts]
                    rValues$analyses <- rbind(rValues$analyses, csvTab)
                }
            }

        })

        observeEvent(input$removeAnalyses, {
            rValues$analyses <- rValues$analyses[-seq.int(input$analysesHot_select$select$r,
                                                          input$analysesHot_select$select$r2)]
        })

        observeEvent(input$DAMethodButton, {
            dm <- selectDirectory("Select DataAnalysis method")
            if (!is.null(dm))
            {
                if (!file.exists(file.path(dm, "DataAnalysis.method")))
                    showDialog("Invalid DataAnalysis method", "Please select a valid DataAnalysis method!", "")
                else
                    updateTextInput(session, "DAMethod", value = dm)
            }
        })

        output$analysesHot <- renderRHandsontable({
            hot <- do.call(rhandsontable,
                           c(list(rValues$analyses, height = 350, maxRows = nrow(rValues$analyses)),
                             hotOpts)) %>%
                hot_col(c("group", "ref"), readOnly = FALSE, type = "text")

            return(hot)
        })
    }

    runGadget(getNewProjectUI(), server, viewer = dialogViewer("Create new project", width = 800, height = 600))
    # runGadget(getNewProjectUI(), server, viewer = paneViewer())
}
