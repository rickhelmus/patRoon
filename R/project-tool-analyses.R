newProjectAnalysesUI <- function(id)
{
    ns <- NS(id)
    
    anaTableCondition <- function(pol, wh)
    {
        ret <- "input.generateAnaInfo == \"table\""
        if (!missing(pol))
        {
            if (is.null(pol))
                ret <- paste(ret, "&& output.ionization != \"both\"")
            else
            {
                ret <- paste(ret, "&& output.ionization == \"both\"")
                if (pol != "both")
                    ret <- paste0(ret, " && input.currentSetStatic == \"", pol, "\"")
            }
        }
        return(ret)
    }
    
    anaDynCondition <- function(pol, wh)
    {
        ret <- "input.generateAnaInfo == \"dynamic\""
        if (is.null(pol))
            ret <- paste(ret, "&& output.ionization != \"both\"")
        else
        {
            ret <- paste(ret, "&& output.ionization == \"both\"")
            if (pol != "both")
                ret <- paste0(ret, " && input.currentSetDynamic == \"", pol, "\"")
        }
        return(ret)
    }
    
    getAnaTableFileUI <- function(pol)
    {
        # NOTE: pol may be NULL
        pf <- if (is.null(pol)) "" else if (pol == "positive") "Pos" else "Neg"
        
        fillRow(
            flex = NA,
            conditionalPanel(
                condition = "input.analysisTableFileType == \"CSV\"",
                ns = ns,
                textInput(ns(paste0("analysisTableFileCSV", pf)), "File name", width = "100%")
            ),
            conditionalPanel(
                condition = "input.analysisTableFileType == \"R\"",
                ns = ns,
                textInput(ns(paste0("analysisTableFileR", pf)), "File name", width = "100%")
            )
        )
    }
    
    genAnaDynUI <- function(pol)
    {
        pf <- if (is.null(pol)) "" else if (pol == "positive") "Pos" else "Neg"
        
        htmltools::tagList(
            fileSelect(ns(paste0("genAnaInfoDynRaw", pf)), ns("genAnaInfoDynRawButton"),
                       "raw analyses", "", placeholder = "Leave empty for none"),
            fileSelect(ns(paste0("genAnaInfoDynCentroid", pf)), ns("genAnaInfoDynCentroidButton"),
                       "centroided analyses", "", placeholder = "Leave empty for none"),
            fileSelect(ns(paste0("genAnaInfoDynIMS", pf)), ns("genAnaInfoDynIMSButton"),
                       "IMS analyses", "", placeholder = "Leave empty for none"),
            fileSelect(ns(paste0("genAnaInfoDynProfile", pf)), ns("genAnaInfoDynProfileButton"),
                       "profile analyses", "", placeholder = "Leave empty for none")
        )
    }
    
    htmltools::tagList(
        miniUI::miniContentPanel(
            fillCol(
                flex = NA,
                fillRow(
                    height = 180,
                    fillCol(
                        flex = NA,
                        radioButtons(ns("generateAnaInfo"), "Generate analysis information",
                                     c("None" = "none", "Analyses table" = "table",
                                       "Dynamically in script" = "dynamic", "Example data" = "example")),
                        conditionalPanel(
                            condition = anaTableCondition(),
                            ns = ns,
                            radioButtons(ns("analysisTableFileType"), "Save analyses table",
                                         c("as CSV file" = "CSV", "as R file" = "R", "embedded in script" = "embedded"),
                                         inline = TRUE),
                        ),
                        conditionalPanel(
                            condition = "input.generateAnaInfo == \"example\"",
                            ns = ns,
                            fillRow(
                                height = 30,
                                textNote("Make sure that the patRoonData package is installed.")
                            )
                        )
                    ),
                    fillCol(
                        flex = NA,
                        conditionalPanel(
                            condition = anaTableCondition("both"),
                            ns = ns,
                            radioButtons(ns("currentSetStatic"), "Define table for set", c("positive", "negative"),
                                         inline = TRUE)
                        ),
                        conditionalPanel(
                            condition = anaTableCondition(NULL),
                            ns = ns,
                            getAnaTableFileUI(NULL)
                        ),
                        conditionalPanel(
                            condition = anaTableCondition("positive"),
                            ns = ns,
                            getAnaTableFileUI("positive")
                        ),
                        conditionalPanel(
                            condition = anaTableCondition("negative"),
                            ns = ns,
                            getAnaTableFileUI("negative")
                        )
                    )
                ),
                conditionalPanel(
                    condition = "input.generateAnaInfo == \"dynamic\"",
                    ns = ns,
                    fillCol(
                        flex = NA,
                        conditionalPanel(
                            condition = anaDynCondition("both"),
                            ns = ns,
                            radioButtons(ns("currentSetDynamic"), "Define table for set", c("positive", "negative"),
                                         inline = TRUE)
                        ),
                        strong(em("Generate analysis information from the following directories")),
                        br(),
                        conditionalPanel(
                            condition = anaDynCondition(NULL),
                            ns = ns,
                            genAnaDynUI(NULL)
                        ),
                        conditionalPanel(
                            condition = anaDynCondition("positive"),
                            ns = ns,
                            genAnaDynUI("positive")
                        ),
                        conditionalPanel(
                            condition = anaDynCondition("negative"),
                            ns = ns,
                            genAnaDynUI("negative")
                        )
                    )
                ),
                conditionalPanel(
                    height = 30,
                    condition = paste(anaTableCondition(), "|| input.generateAnaInfo == \"dynamic\""),
                    ns = ns,
                    textNote(paste("Raw analyses are instrument files (.raw, .d etc);",
                                   "others are mzXML/mzML files.",
                                   "Exporting data is configured in Data pre-treatment."))
                ),
                conditionalPanel(
                    condition = anaTableCondition(),
                    ns = ns,
                    rhandsontable::rHandsontableOutput(ns("analysesHOT"))
                )
            )
        ),
        
        conditionalPanel(
            condition = anaTableCondition(),
            ns = ns,
            miniUI::miniButtonBlock(
                actionButton(ns("setAnaInfoFiles"), label = "Files", icon = icon("file-import"),
                             title = "Set analysis files"),
                shinyjs::disabled(actionButton(ns("removeAnaInfoRows"), label = "Remove", icon = icon("trash"),
                                               title = "Remove selected row(s)")),
                shinyjs::disabled(actionButton(ns("anaInfoRowsUp"), label = "Move up", icon = icon("arrow-up"),
                                               title = "Move selected row(s) up")),
                shinyjs::disabled(actionButton(ns("anaInfoRowsDown"), label = "Move down", icon = icon("arrow-down"),
                                               title = "Move selected row(s) down")),
                actionButton(ns("importAnaInfoCSV"), label = "Import CSV", icon = icon("file-csv"),
                             title = "Import previously generated analyses information from a csv file")
            )
        )
    )
}

newProjectAnalysesServer <- function(id, ionization, settings)
{
    ns <- NS(id)
    
    anaFilesToAnaInfo <- function(anaFiles)
    {
        ret <- dcast(anaFiles, analysis ~ type, value.var = "path")
        setnames(ret, getMSFileTypes(), paste0("path_", getMSFileTypes()), skip_absent = TRUE)
        for (ft in getMSFileTypes())
        {
            col <- paste0("path_", ft)
            if (is.null(ret[[col]]))
                set(ret, j = col, value = character(0))
        }
        return(ret)
    }
    
    anaInfoToAnaFiles <- function(anaInfo)
    {
        ret <- melt(anaInfo, id.vars = "analysis", measure.vars = getAnaInfoPathCols(anaInfo), variable.name = "type",
                    value.name = "path", na.rm = TRUE)
        ret[, type := gsub("^path_", "", type)]
        return(ret[])
    }
    
    emptyAnaTable <- function()
    {
        ai <- data.table(analysis = character(0), type = character(0), replicate = character(0), blank = character(0),
                         conc = numeric(0), norm_conc = numeric(0))
        ai[, paste0("path_", getMSFileTypes()) := character(0)]
        return(ai)
    }
    
    moduleServer(id, function(input, output, session)
    {
        rValues <- reactiveValues(triggerAnaInfoHOTUpdate = 0,
                                  anaInfoTabs = list(positive = emptyAnaTable(), negative = emptyAnaTable()),
                                  analysisFiles = data.table(analysis = character(), type = character(), path = character()))
        
        observeEvent(settings(), {
            updateRadioButtons(session, "generateAnaInfo", selected = settings()$generateAnaInfo)
            updateRadioButtons(session, "analysisTableFileType", selected = settings()$analysisTableFileType)
            updateTextInput(session, "analysisTableFileCSV", value = settings()$analysisTableFileCSV)
            updateTextInput(session, "analysisTableFileCSVPos", value = settings()$analysisTableFileCSVPos)
            updateTextInput(session, "analysisTableFileCSVNeg", value = settings()$analysisTableFileCSVNeg)
            updateTextInput(session, "analysisTableFileR", value = settings()$analysisTableFileR)
            updateTextInput(session, "analysisTableFileRPos", value = settings()$analysisTableFileRPos)
            updateTextInput(session, "analysisTableFileRNeg", value = settings()$analysisTableFileRNeg)
            updateTextInput(session, "genAnaInfoDynRaw", value = settings()$genAnaInfoDynRaw)
            updateTextInput(session, "genAnaInfoDynRawPos", value = settings()$genAnaInfoDynRawPos)
            updateTextInput(session, "genAnaInfoDynRawNeg", value = settings()$genAnaInfoDynRawNeg)
            updateTextInput(session, "genAnaInfoDynCentroid", value = settings()$genAnaInfoDynCentroid)
            updateTextInput(session, "genAnaInfoDynCentroidPos", value = settings()$genAnaInfoDynCentroidPos)
            updateTextInput(session, "genAnaInfoDynCentroidNeg", value = settings()$genAnaInfoDynCentroidNeg)
            updateTextInput(session, "genAnaInfoDynIMS", value = settings()$genAnaInfoDynIMS)
            updateTextInput(session, "genAnaInfoDynIMSPos", value = settings()$genAnaInfoDynIMSPos)
            updateTextInput(session, "genAnaInfoDynIMSNeg", value = settings()$genAnaInfoDynIMSNeg)
            updateTextInput(session, "genAnaInfoDynProfile", value = settings()$genAnaInfoDynProfile)
            updateTextInput(session, "genAnaInfoDynProfilePos", value = settings()$genAnaInfoDynProfilePos)
            updateTextInput(session, "genAnaInfoDynProfileNeg", value = settings()$genAnaInfoDynProfileNeg)
        })
        
        makeAnalysesHOT <- function()
        {
            rValues$triggerAnaInfoHOTUpdate
            
            dt <- isolate(copy(rValues$anaInfoTabs[[getCurPolarity()]]))
            pcols <- getAnaInfoPathCols(dt)
            dt[, type := {
                p <- unlist(mget(pcols))
                paste0(sub("^path_", "", names(p)[!is.na(p)]), collapse = ", ")
            }, by = .I]
            dt[, (pcols) := NULL]
            hot <- makeNewProjectHOT(dt, height = 250, maxRows = nrow(dt), columnSorting = FALSE) |>
                rhandsontable::hot_col(c("replicate", "blank"), readOnly = FALSE, type = "text") |>
                rhandsontable::hot_col(c("conc", "norm_conc"), readOnly = FALSE, type = "numeric")
            return(hot)
        }
        
        triggerAnaInfoHOTUpdate <- function() rValues$triggerAnaInfoHOTUpdate <- rValues$triggerAnaInfoHOTUpdate + 1
        
        makeAnalysisFilesHOT <- function()
        {
            tab <- copy(rValues$analysisFiles)
            tab[, format := mapply(analysis, path, type, FUN = function(a, p, t) {
                exts <- unique(MSFileExtensions()[getMSFileFormats(t)])
                paths <- file.path(p, paste0(a, ".", exts))
                exts <- exts[file.exists(paths)]
                if (length(exts) == 0)
                    "not found!"
                else
                    paste0(exts, collapse = ", ")
            })]
            setcolorder(tab, c("analysis", "type", "format", "path"))
            hot <- makeNewProjectHOT(tab, height = 300, maxRows = nrow(rValues$analysisFiles), rowHeaders = NULL)
            return(hot)
        }
        
        getCurPolarity <- function()
        {
            if (ionization() == "positive" || (ionization() == "both" && input$currentSetStatic == "positive"))
                return("positive")
            return("negative")
        }
        
        moveSelectedAnalyses <- function(dir)
        {
            # NOTE: assume selection is a block range
            sel <- input$analysesHOT_select$select$rAll
            pol <- getCurPolarity()
            rValues$anaInfoTabs[[pol]] <- moveSelectedTabRows(dir, sel, rValues$anaInfoTabs[[pol]])
            triggerAnaInfoHOTUpdate()
            moveHOTSel(session, dir, sel, ns("analysesHOT"))
        }
        
        observeEvent(input$analysesHOT, {
            # HACK: maxRows: make sure we don't have empty table as hot_to_r errors otherwise
            if (!is.null(input$analysesHOT) && input$analysesHOT$params$maxRows > 0)
            {
                # sync from table edits
                dt <- rhandsontable::hot_to_r(input$analysesHOT)
                pol <- getCurPolarity()
                ai <- rValues$anaInfoTabs[[pol]]
                mcols <- c("replicate", "blank", "conc", "norm_conc")
                ai[, (mcols) := dt[, mcols, with = FALSE]]
                rValues$anaInfoTabs[[pol]] <- ai
            }
            else
            {
                shinyjs::disable("removeAnaInfoRows")
                shinyjs::disable("anaInfoRowsUp")
                shinyjs::disable("anaInfoRowsDown")
            }
        })
        
        observeEvent(input$analysesHOT_select$select$r, {
            pol <- getCurPolarity()
            sel <- input$analysesHOT_select$select$rAll
            e <- nrow(rValues$anaInfoTabs[[pol]]) > 0 && length(sel) > 0
            shinyjs::toggleState("removeAnaInfoRows", condition = e)
            shinyjs::toggleState("anaInfoRowsUp", condition = e && min(sel) > 1)
            shinyjs::toggleState("anaInfoRowsDown", condition = e && max(sel) < nrow(rValues$anaInfoTabs[[pol]]))
        })
        
        observeEvent(input$setAnaInfoFiles, {
            pol <- getCurPolarity()
            rValues$analysisFiles <- if (nrow(rValues$anaInfoTabs[[pol]]) == 0)
                rValues$analysisFiles[0, ]
            else
                anaInfoToAnaFiles(rValues$anaInfoTabs[[pol]])
            
            showModal(modalDialog(
                title = "Select analyses",
                rhandsontable::rHandsontableOutput(ns("analysisFilesHOT"), height = 300, width = 700),
                easyClose = TRUE,
                fade = FALSE, # TRUE messes up HOT
                footer = tagList(
                    # use column() to (1) make sure selectInput only occupies a single row and (2) a button can
                    # be placed right next to it and (3) left alignment
                    column(
                        width = 2,
                        style = "padding-right: 0px;",
                        htmltools::tagAppendAttributes(
                            selectInput(ns("analysisFilesAddType"), NULL, c("raw", "centroid", "profile", "ims"),
                                        selectize = FALSE, width = "100%"),
                            style = "margin-bottom: 0px;")
                    ),
                    column(
                        width = 1,
                        style = "padding-left: 0px;",
                        title = "Add all analysis files with the selected type within a directory",
                        actionButton(ns("analysisFilesAdd"), label = NULL, icon = icon("plus"))
                    ),
                    column(
                        width = 2,
                        style = "padding-right: 0px;",
                        htmltools::tagAppendAttributes(
                            selectInput(ns("analysisFilesChangeType"), NULL, c("centroid", "profile", "ims"),
                                        selectize = FALSE, width = "100%"),
                            style = "margin-bottom: 0px;")
                    ),
                    column(
                        width = 1,
                        style = "padding-left: 0px;",
                        title = "Change type of selected analysis files",
                        shinyjs::disabled(actionButton(ns("analysisFilesChange"), label = NULL, icon = icon("pen-to-square")))
                    ),
                    column(
                        width = 1,
                        shinyjs::disabled(actionButton(ns("analysisFilesRemove"), label = NULL, icon = icon("trash")))
                    ),
                    modalButton("Cancel"),
                    actionButton(ns("analysisFilesOK"), "OK")
                )
            ))
        })
        
        observeEvent(input$removeAnaInfoRows, {
            sel <- input$analysesHOT_select$select$rAll
            pol <- getCurPolarity()
            rValues$anaInfoTabs[[pol]] <- rValues$anaInfoTabs[[pol]][-sel]
            triggerAnaInfoHOTUpdate()
        })
        
        observeEvent(input$anaInfoRowsUp, { moveSelectedAnalyses("up") })
        observeEvent(input$anaInfoRowsDown, { moveSelectedAnalyses("down") })
        
        observeEvent(input$analysisFilesHOT, {
            # HACK: maxRows: make sure we don't have empty table as hot_to_r errors otherwise
            if (is.null(input$analysesHOT) || input$analysisFilesHOT$params$maxRows == 0)
            {
                shinyjs::disable("analysisFilesChangeType")
                shinyjs::disable("analysisFilesChange")
                shinyjs::disable("analysisFilesRemove")
            }
        })
        observeEvent(input$analysisFilesHOT_select$select$r, {
            if (nrow(rValues$analysisFiles) > 0 && length(input$analysisFilesHOT_select$select$rAll) > 0)
            {
                sel <- input$analysisFilesHOT_select$select$rAll
                rawSelected <- "raw" %in% rValues$analysisFiles[sel, "type"]
                shinyjs::toggleState("analysisFilesChangeType", condition = !rawSelected)
                shinyjs::toggleState("analysisFilesChange", condition = !rawSelected)
                
                # Disallow IMS selection for mzXML files (not supported by format)
                mzXMLSelected <- "mzxml" %in% tolower(tools::file_ext(rValues$analysisFiles[sel ,"path"]))
                ch <- c("centroid", "profile", if (!mzXMLSelected) "ims")
                s <- rValues$analysisFiles[sel ,"type"]
                if (uniqueN(s) > 1)
                    s <- ch[1] # just default to first choice if selection isn't homogenous
                updateSelectInput(session, "analysisFilesChangeType", selected = s, choices = ch)
                
                shinyjs::enable("analysisFilesRemove")
            }
            else
            {
                shinyjs::disable("analysisFilesChange")
                shinyjs::disable("analysisFilesRemove")
            }
        })
        
        observeEvent(input$analysisFilesAdd, {
            anaDir <- rstudioapi::selectDirectory(path = "~/")
            if (!is.null(anaDir))
            {
                files <- listMSFiles(anaDir, getMSFileFormats(input$analysisFilesAddType))
                af <- data.table(analysis = tools::file_path_sans_ext(basename(files)),
                                 type = if (length(files) > 0) input$analysisFilesAddType else character(),
                                 path = dirname(files))
                af <- unique(af) # in case of multiple files with the same analysis name and path (eg mzXML+mzML)
                rValues$analysisFiles <- rbind(rValues$analysisFiles, af)
            }
        })
        
        observeEvent(input$analysisFilesChange, {
            sel <- input$analysisFilesHOT_select$select$rAll
            rValues$analysisFiles$type[sel] <- input$analysisFilesChangeType
        })
        
        observeEvent(input$analysisFilesRemove, {
            sel <- input$analysisFilesHOT_select$select$rAll
            rValues$analysisFiles <- rValues$analysisFiles[-sel]
        })
        
        observeEvent(input$analysisFilesOK, {
            removeModal()
            pol <- getCurPolarity()
            
            if (nrow(rValues$analysisFiles) > 0)
            {
                ai <- anaFilesToAnaInfo(rValues$analysisFiles)
                
                # remove the analyses for which all file paths were removed
                rValues$anaInfoTabs[[pol]] <- rValues$anaInfoTabs[[pol]][analysis %in% rValues$analysisFiles$analysis]
                
                # overlap: update paths
                aiOv <- ai[analysis %in% rValues$anaInfoTabs[[pol]]$analysis]
                if (nrow(aiOv) > 0)
                {
                    pcols <- getAnaInfoPathCols(rValues$anaInfoTabs[[pol]])
                    rValues$anaInfoTabs[[pol]][aiOv, (pcols) := mget(paste0("i.", pcols)), on = "analysis"]
                }
                
                # add analyses for new files                
                aiNew <- ai[!analysis %in% rValues$anaInfoTabs[[pol]]$analysis]
                if (nrow(aiNew) > 0)
                    rValues$anaInfoTabs[[pol]] <- rbind(rValues$anaInfoTabs[[pol]], aiNew, fill = TRUE)
            }
            else
                rValues$anaInfoTabs[[pol]] <- rValues$anaInfoTabs[[pol]][0]
            
            triggerAnaInfoHOTUpdate()
        })
        
        observeEvent(input$importAnaInfoCSV, {
            csvFile <- rstudioapi::selectFile(path = "~/", filter = "csv files (*.csv)")
            if (!is.null(csvFile))
            {
                ai <- tryCatch(assertAndPrepareAnaInfo(fread(csvFile)), error = function(e)
                {
                    rstudioapi::showDialog("Invalid CSV file",
                                           sprintf("The selected CSV file is not a valid analysis information file: '%s'", e))
                    return(FALSE)
                })
                
                if (!isFALSE(ai) && nrow(ai) > 0)
                {
                    pol <- getCurPolarity()
                    
                    # check for overlap
                    if (any(ai$analysis %in% rValues$anaInfoTabs[[pol]]$analysis))
                    {
                        overwrite <- rstudioapi::showQuestion("Overwrite analysis information",
                                                              "One or more analyses are already defined. Do you want to overwrite or keep the current information for these analyses?",
                                                              "Overwrite", "Keep", timeout = NULL)
                        if (overwrite) # UNDONE: try to keep original order?
                            rValues$anaInfoTabs[[pol]] <- rValues$anaInfoTabs[[pol]][!analysis %in% ai$analysis]
                        else
                            ai <- ai[!analysis %in% rValues$anaInfoTabs[[pol]]$analysis]
                    }
                    
                    rValues$anaInfoTabs[[pol]] <- rbind(rValues$anaInfoTabs[[pol]], ai, fill = TRUE)
                    triggerAnaInfoHOTUpdate()
                }
            }
        })
        
        doObserveGenAnaInfoDynSelDir <- function(textID, buttonID)
        {
            observeEvent(input[[buttonID]], {
                pf <- if (input$ionization != "both") "" else if (input$currentSetDynamic == "positive") "Pos" else "Neg"
                textID <- paste0(textID, pf)
                d <- rstudioapi::selectDirectory("Select directory", path = input[[textID]])
                if (!is.null(d))
                    updateTextInput(session, textID, value = d)
            })
        }
        
        doObserveGenAnaInfoDynSelDir("genAnaInfoDynRaw", "genAnaInfoDynRawButton")
        doObserveGenAnaInfoDynSelDir("genAnaInfoDynCentroid", "genAnaInfoDynCentroidButton")
        doObserveGenAnaInfoDynSelDir("genAnaInfoDynIMS", "genAnaInfoDynIMSButton")
        doObserveGenAnaInfoDynSelDir("genAnaInfoDynProfile", "genAnaInfoDynProfileButton")
        
        doObserveSelDir(input, session, "destinationPath", "destinationPathButton")
        
        output$analysesHOT <- rhandsontable::renderRHandsontable(makeAnalysesHOT())
        output$analysisFilesHOT <- rhandsontable::renderRHandsontable(makeAnalysisFilesHOT())
        
        output <- exportShinyOutputVal(output, "ionization", ionization)
        
        list(
            valid = reactive({
                ret <- TRUE
                if (input$generateAnaInfo == "table")
                {
                    for (pol in c("positive", "negative"))
                    {
                        if (nrow(anaInfoTabs[[pol]]) == 0 && input$ionization %in% c(pol, "both"))
                        {
                            ret <- list(title = "No analyses", msg = paste0("Please select some analyses for ", pol, " mode!"))
                            break
                        }
                    }
                }
                ret
            }),
            settings = reactive(list(
                generateAnaInfo = input$generateAnaInfo,
                analysisTableFileType = input$analysisTableFileType,
                analysisTableFileCSV = input$analysisTableFileCSV,
                analysisTableFileCSVPos = input$analysisTableFileCSVPos,
                analysisTableFileCSVNeg = input$analysisTableFileCSVNeg,
                analysisTableFileR = input$analysisTableFileR,
                analysisTableFileRPos = input$analysisTableFileRPos,
                analysisTableFileRNeg = input$analysisTableFileRNeg,
                genAnaInfoDynRaw = input$genAnaInfoDynRaw,
                genAnaInfoDynRawPos = input$genAnaInfoDynRawPos,
                genAnaInfoDynRawNeg = input$genAnaInfoDynRawNeg,
                genAnaInfoDynCentroid = input$genAnaInfoDynCentroid,
                genAnaInfoDynCentroidPos = input$genAnaInfoDynCentroidPos,
                genAnaInfoDynCentroidNeg = input$genAnaInfoDynCentroidNeg,
                genAnaInfoDynIMS = input$genAnaInfoDynIMS,
                genAnaInfoDynIMSPos = input$genAnaInfoDynIMSPos,
                genAnaInfoDynIMSNeg = input$genAnaInfoDynIMSNeg,
                genAnaInfoDynProfile = input$genAnaInfoDynProfile,
                genAnaInfoDynProfilePos = input$genAnaInfoDynProfilePos,
                genAnaInfoDynProfileNeg = input$genAnaInfoDynProfileNeg
            )),
            anaInfoTabs = reactive(rValues$anaInfoTabs)
        )
    })
}

defaultAnalysesSettings <- function()
{
    return(list(
        generateAnaInfo = "none",
        analysisTableFileType = "CSV",
        analysisTableFileCSV = "analyses.csv",
        analysisTableFileCSVPos = "analyses-pos.csv",
        analysisTableFileCSVNeg = "analyses-neg.csv",
        analysisTableFileR = "analyses.R",
        analysisTableFileRPos = "analyses-pos.R",
        analysisTableFileRNeg = "analyses-neg.R",
        genAnaInfoDynRaw = "",
        genAnaInfoDynRawPos = "",
        genAnaInfoDynRawNeg = "",
        genAnaInfoDynCentroid = "",
        genAnaInfoDynCentroidPos = "",
        genAnaInfoDynCentroidNeg = "",
        genAnaInfoDynIMS = "",
        genAnaInfoDynIMSPos = "",
        genAnaInfoDynIMSNeg = "",
        genAnaInfoDynProfile = "",
        genAnaInfoDynProfilePos = "",
        genAnaInfoDynProfileNeg = ""
    ))
}

upgradeAnalysesSettings <- function(settings)
{
    # NOTE: this updates from first file version
    return(modifyList(defaultAnalysesSettings(), list(
        analysisTableFileCSV = settings$analysisTableFile,
        analysisTableFileCSVPos = settings$analysisTableFilePos,
        analysisTableFileCSVNeg = settings$analysisTableFileNeg
    )))
}
