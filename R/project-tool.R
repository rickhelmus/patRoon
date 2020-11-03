#' @include main.R
NULL

getScriptCode <- function(input, analyses)
{
    optionalCodeBlock <- function(e) if (e) "<<startCodeBlock>>" else "<<skipCodeBlock>>"
    endCodeBlock <- function() "<<endCodeBlock>>"
    optionalLine <- function(e) if (!e) "<<skipThisLine>>" else ""
    header <- function(title)
    {
        hd <- paste("#", strrep("-", 25))
        paste(hd, paste("#", title), hd, sep = "\n")
    }

    if (input$peakPicking)
        centroid <- if (input$peakPickingVendor && input$convAlgo == "pwiz") "\"vendor\"" else TRUE
    else
        centroid <- FALSE

    # Can't set lists in template() call (bug?)
    preTreatOpts = list(convAlgo = input$convAlgo,
                        convFrom = input$convFrom,
                        convTo = input$convTo, centroid = centroid,
                        DAMethod = input$DAMethod, doDACalib = input$doDACalib,
                        do = nzchar(input$convAlgo) || nzchar(input$DAMethod) || input$doDACalib)
    featFinderOpts = list(algo = input$featFinder)
    featGrouperOpts = list(algo = input$featGrouper)

    filterFGroupsOpts = list(preIntThr = input$preIntThr, intThr = input$intThr,
                             repAbundance = input$repAbundance, maxRepRSD = input$maxRepRSD,
                             blankThr = input$blankThr, removeBlanks = input$removeBlanks)
    filterFGroupsOpts <- lapply(filterFGroupsOpts, function(x) if (is.numeric(x) && x == 0) "NULL" else x)

    retRange <- c(input[["retention-min"]], input[["retention-max"]])
    if (all(retRange == 0))
        retRange <- NULL
    else if (retRange[2] == 0)
        retRange[2] <- Inf
    mzRange <- c(input[["mz-min"]], input[["mz-max"]])
    if (all(mzRange == 0))
        mzRange <- NULL
    else if (mzRange[2] == 0)
        mzRange[2] <- Inf
    filterFGroupsOpts$retRange <- retRange; filterFGroupsOpts$mzRange <- mzRange

    peakListOpts = list(algo = input$peakListGen)
    formulaOpts = list(algo = input$formulaGen)
    identOpts = list(algo = input$compIdent)
    componentOpts = list(algo = input$components)
    sldest <- gsub("\\", "/", input$destinationPath, fixed = TRUE) # BUG: tmpl() seems to eat backslashes!? Convert to forward slashes for now

    template <- templates::tmpl(readAllFile(system.file("templates", "main_script.R", package = "patRoon")),
                                destination = sldest, generateAnaInfo = input$generateAnaInfo, analysisTableFile = input$analysisTableFile,
                                analyses = analyses, suspectList = input$suspectList, suspectAdduct = input$suspectAdduct,
                                doMSPeakFind = (input$formulaGen != "" && input$formulaGen != "Bruker") || input$compIdent != "",
                                precursorMzWindow = input$precursorMzWindow, polarity = input$polarity,
                                annotateSus = input$annotateSus, genIDLevelFile = input$genIDLevelFile,
                                reportFormats = input$report)

    ret <- template

    # remove blocks marked as '<<skipCodeBlock>>' till <<endCodeBlock>> + leading whitespace
    ret <- gsub("[\t ]*<<skipCodeBlock>>[\\s\\S]*?<<endCodeBlock>>[.\r]*\n?", "", ret, perl = TRUE)

    # remove remaining startCodeBlock and endCodeBlock markers
    ret <- gsub("[\t ]*(<<startCodeBlock>>|<<endCodeBlock>>)[.\r]*\n?", "", ret, perl = TRUE)

    # remove disabled optional lines
    ret <- gsub("([\n]*).*<<skipThisLine>>[.\r]*\n?", "\\1", ret, perl = TRUE)

    # carriage returns seem to mess up cat when writing code
    ret <- gsub("[\r]", "", ret)

    return(ret)
}

doCreateProject <- function(input, analyses)
{
    mkdirp(input$destinationPath)

    analyses <- copy(analyses)
    analyses[, group := ifelse(!nzchar(group), analysis, group)]

    # Make analysis table
    if (input$generateAnaInfo == "table")
        write.csv(analyses[, c("path", "analysis", "group", "blank")],
                  file.path(input$destinationPath, input$analysisTableFile), row.names = FALSE)

    if (input$genIDLevelFile)
        write.csv(defaultIDLevelRules(),
                  file.path(input$destinationPath, "idlevelrules.csv"), row.names = FALSE)
    
    code <- getScriptCode(input, analyses)
    if (input$outputScriptTo == "curFile")
    {
        # insert at end of current document
        rstudioapi::insertText(Inf, code, getSourceEditorContext()$id)
    }
    else
    {
        sp <- file.path(input$destinationPath, input$scriptFile)
        cat(code, file = sp, sep = "")

        if (input$createRStudioProj)
        {
            rstudioapi::initializeProject(input$destinationPath)
            rstudioapi::openProject(input$destinationPath)
        }
        else
            rstudioapi::navigateToFile(sp)
    }
}

getNewProjectUI <- function(destPath)
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
    rangeNumeric <- function(id, label, minVal = 0, maxVal = 0, ...)
    {
        fillRow(
            numericInput(paste0(id, "-min"), paste("Min.", label), value = minVal, ..., width = "95%"),
            numericInput(paste0(id, "-max"), paste("Max.", label), value = maxVal, ..., width = "100%")
        )
    }

    miniUI::miniPage(
        miniUI::gadgetTitleBar("Create project tool", right = miniUI::miniTitleBarButton("create", "Create", TRUE)),

        miniUI::miniTabstripPanel(
            miniUI::miniTabPanel(
                "Destination", icon = icon("save"),
                miniUI::miniContentPanel(
                    fillCol(
                        fileSelect("destinationPath", "projectDestButton", "Project destination",
                                   if (is.null(destPath)) "~/" else destPath),
                        fillRow(
                            radioButtons("outputScriptTo", "Insert code into", c("New file" = "newFile",
                                                                                 "Current file" = "curFile")),
                            conditionalPanel(
                                condition = "input.outputScriptTo == \"newFile\"",
                                textInput("scriptFile", "Script file", "process.R"),
                                checkboxInput("createRStudioProj", "Create (and open) RStudio project", value = TRUE)
                            )
                        )
                    )
                )
            ),
            
            miniUI::miniTabPanel(
                "Analyses", icon = icon("folder-open"),
                miniUI::miniContentPanel(
                    fillCol(
                        flex = c(NA, NA, 1),
                        fillRow(
                            height = 120,
                            radioButtons("generateAnaInfo", "Generate analysis information",
                                         c("None" = "none", "From new csv file" = "table",
                                           "Load in script" = "script", "Example data" = "example")),
                            conditionalPanel(
                                condition = "input.generateAnaInfo == \"table\"",
                                textInput("analysisTableFile", "Analysis table output file", "analyses.csv")
                            )
                        ),
                        conditionalPanel(
                            condition = "input.generateAnaInfo == \"table\" || input.generateAnaInfo == \"script\"",
                            fillRow(
                                height = 30,
                                textNote("Make sure to consider data conversion if data files are not yet in mzXML/mzML format.")
                            )
                        ),
                        conditionalPanel(
                            condition = "input.generateAnaInfo == \"example\"",
                            fillRow(
                                height = 30,
                                textNote("Make sure that the patRoonData package is installed.")
                            )
                        ),
                        conditionalPanel(
                            condition = "input.generateAnaInfo == \"table\" || input.generateAnaInfo == \"script\"",
                            rhandsontable::rHandsontableOutput("analysesHot")
                        )
                    )
                ),
                
                conditionalPanel(
                    condition = "input.generateAnaInfo == \"table\" || input.generateAnaInfo == \"script\"",
                    miniUI::miniButtonBlock(
                        actionButton("addAnalysesDir", "Add analyses from directory ..."),
                        actionButton("addAnalysesCSV", "Add analyses from csv file ..."),
                        actionButton("removeAnalyses", "Remove")
                    )
                )
            ),
            miniUI::miniTabPanel(
                "Data pre-treatment", icon = icon("upload"),
                miniUI::miniContentPanel(
                    fillCol(
                        flex = NA,
                        height = 220,
                        fillRow(
                            height = 75,
                            selectInput("convAlgo", "Data Conversion Algorithm", c("None" = "", "ProteoWizard" = "pwiz",
                                                                                   "Bruker DataAnalysis" = "bruker",
                                                                                   "OpenMS" = "openms"),
                                        width = "100%")
                        ),
                        conditionalPanel(
                            "input.convAlgo != \"\"",
                            fillRow(
                                height = 90,
                                selectInput("convFrom", "Input format", MSFileFormats(), multiple = FALSE,
                                            width = "95%"),
                                fillCol(
                                    flex = c(1, NA),
                                    selectInput("convTo", "Output format", "mzML", multiple = FALSE,
                                                selected = "mzML", width = "100%"),
                                    textNote("enviPick/XCMS support mzXML, XCMS/OpenMS support mzML")
                                )
                            )
                        ),
                        conditionalPanel(
                            "input.convAlgo == \"pwiz\" || input.convAlgo == \"bruker\"",
                            fillRow(
                                height = 60,
                                checkboxInput("peakPicking", "Perform peak picking (line spectra)", value = TRUE),
                                conditionalPanel(
                                    condition = "input.convAlgo == \"pwiz\"",
                                    checkboxInput("peakPickingVendor", "Use vendor algorithm for peak picking", value = TRUE)
                                )
                            )
                        )
                    ),
                    hr(),
                    fillCol(
                        flex = NA,
                        height = 70,
                        strong("Bruker DataAnalysis options"),
                        textNote("Only supported with bruker data and if DataAnalysis is installed.")
                    ),
                    fillCol(
                        flex = NA,
                        height = 120,
                        fillCol(
                            flex = c(1, NA),
                            height = 90,
                            fileSelect("DAMethod", "DAMethodButton", "Set DataAnalysis method"),
                            textNote("Leaving this blank will not set any method")
                        ),
                        checkboxInput("doDACalib", "Perform m/z re-calibration")
                    )
                )
            ),
            miniUI::miniTabPanel(
                "Features", icon = icon("chart-area"),
                miniUI::miniContentPanel(
                    fillRow(
                        height = 75,
                        selectInput("featFinder", "Feature finder", c("OpenMS", "XCMS", "enviPick",
                                                                      "Bruker DataAnalysis" = "Bruker"),
                                    "OpenMS", FALSE, width = "95%"),
                        selectInput("featGrouper", "Feature grouper", c("OpenMS", "XCMS"),
                                    "OpenMS", FALSE, width = "100%")
                    ),
                    fillRow(
                        height = 130,
                        fillCol(
                            flex = c(1, NA),
                            height = 110,
                            width = "95%",
                            fileSelect("suspectList", "suspectListButton", "Suspect list"),
                            textNote("Leave blank for no suspect screening (i.e. perform full non-target analysis)")
                        ),
                        fillCol(
                            flex = c(1, NA),
                            height = 130,
                            textInput("suspectAdduct", "Adduct", placeholder = "e.g. [M+H]+ or [M-H]-"),
                            textNote("Used for mass calculation. Can be left empty if ionized mass or adduct data information is available in the suspect list (\"mz\"/\"adduct\" columns).")
                        )
                    ),
                    hr(),
                    fillCol(
                        flex = NA,
                        height = 70,
                        strong("Post-Filtering of feature groups"),
                        textNote("Set below values to zero to disable a particular filter.")
                    ),
                    fillCol(
                        height = 325,
                        fillRow(
                            numericInput("preIntThr", "Pre-Intensity threshold", 1E2, 0, step = 100, width = "95%"),
                            numericInput("intThr", "Intensity threshold", 1E4, 0, step = 1000, width = "100%")
                        ),
                        fillRow(
                            numericInput("repAbundance", "Min. replicate abundance (relative)", 1, 0, 1.0, 0.1, width = "95%"),
                            numericInput("maxRepRSD", "Max. replicate intensity RSD", 0.75, 0, step = 0.1, width = "100%")
                        ),
                        fillRow(
                            numericInput("blankThr", "Min. blank threshold", 5, 0, step = 1, width = "95%"),
                            checkboxInput("removeBlanks", "Discard blanks after filtering", TRUE)
                        ),
                        rangeNumeric("retention", "retention time (s)", step = 10),
                        rangeNumeric("mz", "m/z", step = 10)
                    )
                )
            ),
            miniUI::miniTabPanel(
                "Annotation", icon = icon("chart-bar"),
                miniUI::miniContentPanel(
                    fillCol(
                        flex = NA,
                        
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
                            condition = "(input.formulaGen != \"\" && input.formulaGen != \"Bruker\") || input.compIdent != \"\"",
                            fillRow(
                                height = 90,
                                selectInput("peakListGen", "Peak list generator",
                                            c("mzR", "Bruker DataAnalysis" = "Bruker"),
                                            "mzR", multiple = FALSE, width = "95%"),
                                conditionalPanel(
                                    condition = "input.peakListGen == \"mzR\"",
                                    numericInput("precursorMzWindow", "MS/MS precursor m/z search window", 4, width = "100%"),
                                    textNote("The precursor m/z search window applied when finding MS/MS spectra. Set to zero for DIA experiments.")
                                )
                            )
                        ),
                        conditionalPanel(
                            condition = "input.formulaGen != \"\" || input.compIdent != \"\" || input.components != \"\"",
                            selectInput("polarity", "Polarity", c("positive", "negative"), "positive",
                                        multiple = FALSE, width = "100%")
                        ),
                        conditionalPanel(
                            condition = "input.suspectList != \"\" && (input.formulaGen != \"\" || input.compIdent != \"\")",
                            fillRow(
                                height = 65,
                                fillCol(
                                    strong("Suspect annotation"),
                                    checkboxInput("annotateSus", "Annotate suspects", TRUE, width = "100%"),
                                    conditionalPanel(
                                        condition = "input.annotateSus",
                                        checkboxInput("genIDLevelFile", "Generate template file with configurable identification levels",
                                                      TRUE, width = "100%")
                                    )
                                )
                            )
                        )
                    )
                )
            ),
            miniUI::miniTabPanel(
                "Reporting", icon = icon("file-medical-alt"), # scroll, paperclip, newspaper, file-medical-alt,
                miniUI::miniContentPanel(
                    fillCol(
                        checkboxGroupInput("report", "Report generation",
                                           c("CSV (text tables)" = "CSV", "PDF (basic plots)" = "PDF",
                                             "HTML (easy browsable plots, bit slower than PDF)" = "HTML"),
                                           c("CSV", "HTML"), width = "100%")
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
#'   \link[=analysis-information]{analysis information} will be written to a
#'   \file{.csv} file so that it can easily be modified afterwards.
#'
#' @param destPath Set destination path value to this value (useful for
#'   debugging). Set to \code{NULL} for a default value.
#'
#' @rdname GUI-utils
#' @export
newProject <- function(destPath = NULL)
{
    rstudioapi::verifyAvailable()

    # UNDONE: warning/message about empty groups

    # NOTE: disable column sorting so we don't have to worry about getting correct row index
    # (https://github.com/jrowen/rhandsontable/issues/257)
    # NOTE: preventOverflow should not be set when columnSorting = FALSE
    # https://github.com/handsontable/handsontable/issues/4303
    # NOTE: set selectionMode to range as only row series can currently be queried
    # (https://github.com/jrowen/rhandsontable/issues/313)
    hotOpts <- list(rowHeaderWidth = 40, readOnly = TRUE,
                    columnSorting = FALSE, sortIndicator = TRUE, selectCallback = TRUE,
                    currentRowClassName = "currentRow", stretchH = "all",
                    selectionMode = "range", outsideClickDeselects = FALSE,
                    contextMenu = FALSE, manualColumnResize = TRUE)

    server <- function(input, output, session)
    {
        rValues <- reactiveValues(analyses = data.table(analysis = character(0), format = character(0),
                                                        group = character(0), blank = character(0), path = character(0)))

        observeEvent(input$create, {
            if (input$destinationPath == "")
                rstudioapi::showDialog("Invalid destination", "Please select a destination path!", "")
            else if (input$outputScriptTo != "curFile" && input$scriptFile == "")
                rstudioapi::showDialog("No script file", "Please select a destination script file!", "")
            else if (input$generateAnaInfo %in% c("table", "script") && nrow(rValues$analyses) == 0)
                rstudioapi::showDialog("No analyses selected", "Please select some analyses!", "")
            else if (input$generateAnaInfo == "table" && file.exists(file.path(input$destinationPath, input$analysisTableFile)) &&
                     !rstudioapi::showQuestion("Analysis table file already exists",
                                               sprintf("Analysis table file already exists: '%s'.\nOverwrite?",
                                                       file.path(input$destinationPath, input$analysisTableFile)),
                                               "Yes", "No"))
            {}
            else if (input$outputScriptTo != "curFile" && file.exists(file.path(input$destinationPath, input$scriptFile)) &&
                     !rstudioapi::showQuestion("Script file already exists:",
                                               sprintf("Script file already exists: '%s'.\nOverwrite?",
                                                       file.path(input$destinationPath, input$scriptFile)),
                                               "Yes", "No"))
            {}
            else
            {
                doCreateProject(input, rValues$analyses)
                stopApp(TRUE)
            }
        })

        observeEvent(input$projectDestButton, {
            dest <- rstudioapi::selectDirectory("Select destination directory", path = input$destinationPath)
            if (!is.null(dest))
                updateTextInput(session, "destinationPath", value = dest)
        })

        observeEvent(input$analysesHot, {
            # HACK: input$analysesHot$params$maxRows: make sure we don't have empty table as hot_to_r errors otherwise
            if (input$analysesHot$params$maxRows > 0)
            {
                df <- rhandsontable::hot_to_r(input$analysesHot)
                rValues$analyses[, c("group", "blank") := .(df$group, df$blank)]
            }
        })

        observeEvent(input$addAnalysesDir, {
            anaDir <- rstudioapi::selectDirectory(path = "~/")
            if (!is.null(anaDir))
            {
                files <- listMSFiles(anaDir, MSFileFormats())

                if (length(files) > 0)
                {
                    dt <- data.table(path = dirname(files), analysis = simplifyAnalysisNames(files),
                                     group = "", blank = "")

                    fExts <- MSFileExtensions()
                    dt[, format := sapply(tools::file_ext(files), function(ext)
                    {
                        paste0(names(fExts)[sapply(fExts, function(e) ext %in% e)], collapse = "/")
                    })]


                    dt[, format := paste0(.SD$format, collapse = ", "), by = .(path, analysis)]
                    dt <- unique(dt, by = c("analysis", "path"))
                    setcolorder(dt, c("analysis", "format", "group", "blank", "path"))

                    rValues$analyses <- rbind(rValues$analyses, dt)
                }
            }
        })

        observeEvent(input$addAnalysesCSV, {
            csvFile <- rstudioapi::selectFile(path = "~/", filter = "csv files (*.csv)")
            if (!is.null(csvFile))
            {
                csvTab <- tryCatch(fread(csvFile, select = c("path", "analysis", "group", "blank"),
                                         colClasses = "character"),
                                   error = function(e) FALSE, warning = function(w) FALSE)
                if (is.logical(csvTab))
                    rstudioapi::showDialog("Error", "Failed to open/parse selected csv file!", "")
                else if (nrow(csvTab) > 0)
                {
                    msExts <- MSFileExtensions()
                    formats <- mapply(csvTab$analysis, csvTab$path, FUN = function(ana, path)
                    {
                        fp <- file.path(path, ana)
                        ret <- names(msExts)[(sapply(msExts, function(e) any(file.exists(paste0(fp, ".", e)))))]
                        ret <- paste0(ret, collapse = ", ")
                        return(if (length(ret) > 0) ret else "")
                    })

                    csvTab[, format := formats]
                    csvTab <- csvTab[nzchar(format)] # prune unknown files (might have been removed?)
                    rValues$analyses <- rbind(rValues$analyses, csvTab)
                }
            }
        })

        observeEvent(input$removeAnalyses, {
            rValues$analyses <- rValues$analyses[-seq.int(input$analysesHot_select$select$r,
                                                          input$analysesHot_select$select$r2)]
        })

        observeEvent(input$convAlgo, {
            from <- switch(input$convAlgo,
                           pwiz = MSFileFormats("pwiz"),
                           bruker = "bruker",
                           openms = MSFileFormats("openms"),
                           ""
                   )
            sel <- ""
            if (nzchar(input$convAlgo))
                sel <- MSFileFormats(input$convAlgo, input$convAlgo != "openms")[1]

            updateSelectInput(session, "convFrom", choices = from, selected = sel)
        })

        observeEvent(input$DAMethodButton, {
            dm <- rstudioapi::selectDirectory("Select DataAnalysis method")
            if (!is.null(dm))
            {
                if (!file.exists(file.path(dm, "DataAnalysis.method")))
                    rstudioapi::showDialog("Invalid DataAnalysis method", "Please select a valid DataAnalysis method!", "")
                else
                    updateTextInput(session, "DAMethod", value = dm)
            }
        })

        observeEvent(input$suspectListButton, {
            sl <- rstudioapi::selectFile("Select suspect list", filter = "csv files (*.csv)")
            if (!is.null(sl))
            {
                csvTab <- tryCatch(fread(sl), error = function(e) FALSE, warning = function(w) FALSE)
                cols <- names(csvTab)
                massCols <- c("mz", "neutralMass", "formula", "SMILES", "InChI")
                
                err <- NULL
                if (is.logical(csvTab))
                    err <- "Failed to open/parse selected csv file!"
                else if (nrow(csvTab) == 0)
                    err <- "The selected files seems to be empty."
                else if (!"name" %in% cols)
                    err <- "The selected file does not have a name column"
                else if (!any(massCols %in% cols))
                    err <- paste("The selected file should have either one of the columns:",
                                 paste(massCols, collapse = ", "))
                
                if (!is.null(err))                
                    rstudioapi::showDialog("Error", err, "")
                else
                    updateTextInput(session, "suspectList", value = sl)
            }
        })

        output$analysesHot <- rhandsontable::renderRHandsontable({
            hot <- do.call(rhandsontable::rhandsontable,
                           c(list(rValues$analyses, height = 250, maxRows = nrow(rValues$analyses)),
                             hotOpts)) %>%
                rhandsontable::hot_col(c("group", "blank"), readOnly = FALSE, type = "text")
            
            return(hot)
        })
    }

    runGadget(getNewProjectUI(destPath), server, viewer = dialogViewer("Create new project", width = 800, height = 600))
    # runGadget(getNewProjectUI(), server, viewer = paneViewer())
}
