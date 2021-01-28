#' @include main.R
#' @include feature_groups.R

getUISettingsPath <- function()
{
    dirPath <- RUserDir("patRoon", "config")
    mkdirp(dirPath)
    return(file.path(dirPath, "check_features.yml"))
}

getUISettings <- function()
{
    path <- getUISettingsPath()
    if (!file.exists(path))
    {
        ret <- getDefaultUISettings()
        saveUISettings(ret)
    }
    else
        ret <- yaml::read_yaml(path, eval.expr = FALSE)
    ret <- ret[setdiff(names(ret), "version")]
    return(ret)
}

getDefaultUISettings <- function()
{
    return(list(retUnit = "sec", featQuantity = "intensity", fGroupQuantity = "average",
                fGroupColumns = c("retMZ", "estIDLevel", "totalScore"),
                featureColumns = c("retMZ", "quantity", "totalScore")))
}

saveUISettings <- function(settings)
{
    settings$version <- 1 # just store for now, in case if ever needed in the future
    yaml::write_yaml(settings, getUISettingsPath(), indent = 4)
}

#' @export
importCheckFeaturesSession <- function(sessionIn, sessionOut, fGroups, overWrite = FALSE)
{
    # UNDONE: docs
    
    checkmate::assertString(sessionIn, min.chars = 1)
    pathIn <- paste0(sessionIn, ".Rds")
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFileExists(pathIn, "r", .var.name = "session", add = ac)
    checkmate::assertString(sessionOut, min.chars = 1, add = ac)
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertFlag(overWrite, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(fGroups) == 0)
        stop("No feature groups, nothing to do...")
    
    pathOut <- paste0(sessionOut, ".Rds")
    if (file.exists(pathOut) && !overWrite)
        stop("Output session already exists. Set overWrite=TRUE to proceed anyway.")
    
    # settings import:
    # - if file has analyses or feature groups not present in target fGroups: remove
    # - default for any missing features/feature groups 
    
    gNames <- names(fGroups)
    settings <- readRDS(pathIn)
    otherGNames <- names(settings$enabledFeatures)
    commonGNames <- intersect(gNames, otherGNames)
    commonAnalyses <- intersect(analyses(fGroups), settings$enabledFeatures$analysis)
    
    if (length(commonGNames) == 0)
        warning("Imported session doesn't contain any relevant feature groups!")
    if (length(commonAnalyses) == 0)
        warning("Imported session doesn't contain any relevant analyses!")

    # only keep common
    settings$enabledFGroups <- intersect(gNames, settings$enabledFGroups)
    settings$enabledFeatures <- settings$enabledFeatures[, c("analysis", commonGNames), drop = FALSE]
    settings$enabledFeatures <- settings$enabledFeatures[settings$enabledFeatures$analysis %in% commonAnalyses, ,
                                                         drop = FALSE]
    
    # add missing
    missingFGroups <- setdiff(gNames, otherGNames)
    settings$enabledFGroups <- c(settings$enabledFGroups, missingFGroups)
    if (length(missingFGroups) > 0 && nrow(settings$enabledFeatures) > 0)
        settings$enabledFeatures[, missingFGroups] <- TRUE
    missingTbl <- data.frame(analysis = setdiff(analyses(fGroups), settings$enabledFeatures$analysis))
    if (nrow(missingTbl) > 0)
    {
        missingTbl[, gNames] <- TRUE
        if (nrow(settings$enabledFeatures) > 0)
            settings$enabledFeatures <- rbind(settings$enabledFeatures, missingTbl)
        else
            settings$enabledFeatures <- missingTbl
    }
    
    # match analysis order
    settings$enabledFeatures <- settings$enabledFeatures[match(analyses(fGroups), settings$enabledFeatures$analysis), ]
    
    saveRDS(settings, pathOut)
}

getCheckFeatsUI <- function(settings)
{
    showOpts <- c("Keep", "Don't keep")
    
    fillPage(
        tags$head(includeScript(system.file("js", "utils-checkUI.js", package = "patRoon"))),
        shinyjs::useShinyjs(),
        
        keys::useKeys(),
        keys::keysInput("keys", c("p", "n", "t")),
        
        title = "Check features tool",
        
        fillCol(
            flex = c(1, NA),
            
            fillCol(
                flex = c(NA, 1),
                fillRow(
                    flex = c(1, NA),
                    height = 40,
                    strong(style = "font-size: 200%; text-align: center;", textOutput("pageTitle")),
                    fillRow(
                        width = 130,
                        height = 30,
                        shinyjs::disabled(actionButton("saveSession", "Save session", icon("save")))
                    )
                ),
                fillCol(
                    plotOutput("plotChrom", width = "100%", height = "100%")
                )
            ),
            
            fillRow(
                height = 300,

                tabsetPanel(id = "tabs",
                    tabPanel("Feature groups", fillCol(
                        flex = c(NA, 1),
                        fillRow(
                            height = 40,
                            flex = c(NA, NA, NA, NA, NA, 1),
                            
                            fillCol(
                                width = 45,
                                actionButton("previousGroup", "", icon("arrow-left"))
                            ),
                            fillCol(
                                width = 45,
                                actionButton("nextGroup", "", icon("arrow-right"))
                            ),
                            fillCol(
                                width = 100,
                                actionButton("toggleGroup", "Toggle", icon("toggle-on"))
                            ),
                            fillCol(
                                width = 200,
                                div(style = "margin-top: 8px; margin-left: 5px",
                                    checkboxGroupInput("showWhat", NULL, showOpts, showOpts,
                                                       inline = TRUE))
                            ),
                            fillRow(
                                width = 250,
                                flex = c(NA, 1),
                                fillRow(
                                    width = 75,
                                    div(style = "margin-top: 8px", HTML("<strong>Plot mode</strong>"))
                                ),
                                selectInput("fGroupPlotMode", NULL,
                                            c("Top most group" = "topMost", "Top most replicates" = "topMostByRGroup",
                                              "All" = "all"), "topMost")
                            ),
                            div(style = "margin: 8px 10px 12px; font-size: small; text-align: right;",
                                HTML("<b>n</b>: next; <b>p</b>: previous; <b>t</b>: toggle"))
                        ),
                        fillRow(
                            tags$div(id = "fGroupsHotContainer", rhandsontable::rHandsontableOutput("fGroupsHot"))
                        )
                    )),
                    tabPanel("Features", fillCol(
                        flex = c(NA, 1),
                        
                        fillRow(
                            height = 40,
                            actionButton("resetFeatures", "Enable all features for all groups",
                                         icon("check-square"))
                        ),
                        
                        fillRow(
                            tags$div(id = "featuresHotContainer", rhandsontable::rHandsontableOutput("featuresHot"))
                        )
                    )),
                    tabPanel("Settings", fillCol(
                        flex = NA,
                        
                        fillRow(
                            height = 200,
                            
                            fillCol(
                                flex = NA,
                                radioButtons("retUnit", "Retention time unit", c("Seconds" = "sec", "Minutes" = "min"),
                                             settings$retUnit, inline = TRUE),
                                radioButtons("featQuantity", "Feature quantity", c("Peak intensity" = "intensity",
                                                                                   "Peak area" = "area"),
                                             settings$featQuantity, inline = TRUE),
                                radioButtons("fGroupQuantity", "Reported feature group quantities",
                                             c("Max" = "max", "Mean" = "average", "All" = "all"),
                                             settings$fGroupQuantity, inline = TRUE)
                            ),
                            fillRow(
                                checkboxGroupInput("fGroupColumns", "Feature groub table columns",
                                                   c("Retention time & m/z" = "retMZ",
                                                     "EIC preview" = "EICPreview",
                                                     "Suspect properties (RT, m/z, fragments)" = "suspProp",
                                                     "Estimated suspect identification level" = "estIDLevel",
                                                     "Other suspect annotations" = "suspOther",
                                                     "Total score" = "totalScore",
                                                     "Other scores" = "otherScores"),
                                                   settings$fGroupColumns)
                            ),
                            fillRow(
                                checkboxGroupInput("featureColumns", "Feature table columns",
                                                   c("Retention time & m/z" = "retMZ",
                                                     "Replicate group" = "rGroup",
                                                     "Blank" = "blank",
                                                     "Quantity" = "quantity",
                                                     "RT and m/z range" = "rtMZRange",
                                                     "Total score" = "totalScore",
                                                     "Other scores" = "otherScores"),
                                                   settings$featureColumns)
                            )
                        ),
                        fillRow(
                            height = 40,
                            flex = NA,
                            shinyjs::disabled(actionButton("saveApplySettings", "Save & Apply", icon = icon("save"))),
                            actionButton("resetSettings", "Reset defaults", icon = icon("redo"))
                        )
                    ))
                )
            )
        )
    )
}

#' @details \code{checkFeatures} is used to review chromatographic
#'   information for feature groups. This is especially useful to get a visual
#'   impression of the quality of detected features. In addition, this function
#'   may be used to remove unwanted (\emph{e.g.} outlier) features. Better
#'   performance is often obtained when an external browser is used to use this
#'   Shiny application. Furthermore, when a large \code{featureGroups} object is
#'   used it is recommended to limit the number of analyses/feature groups by
#'   subsetting the object.
#'
#' @param mzWindow Default \emph{m/z} window to be used for creating extracted
#'   ion chromatograms (EICs).
#' @param enabledFGroups A logical vector that states for each feature group
#'   whether it should be kept (\code{TRUE}) or not (\code{FALSE}). The order is
#'   the same as the \code{fGroups} parameter. If \code{NULL} then all feature
#'   groups are considered to be kept.
#'
#' @return \code{checkFeatures} returns a logical vector for all feature
#'   groups that were selected to be kept (\code{TRUE}) or not (\code{FALSE}).
#'   This result can be passed to the \code{enabledFGroups} parameter for
#'   subsequent calls to \code{checkFeatures} in order to restore the
#'   keep/not keep state from a previous call. To actually remove unwanted
#'   feature groups the object should be subset by the subsetting
#'   (\code{\link{[}}) operator to which the return value should be passed as
#'   the second parameter.
#'
#' @rdname GUI-utils
#' @aliases checkFeatures
#' @export
setMethod("checkFeatures", "featureGroups", function(fGroups, session, rtWindow, clearSession)
{
    # UNDONE: update docs
    
    if (length(fGroups) == 0)
        stop("No feature groups, nothing to check...")
    
    checkmate::assertFlag(clearSession)
    
    ac <- checkmate::makeAssertCollection()
    assertCheckFeaturesSession(session, fGroups, mustExist = FALSE, canClearSession = TRUE,
                               didClearSession = clearSession, add = ac)
    checkmate::assertNumber(rtWindow, finite = TRUE, lower = 0, add = ac)
    checkmate::reportAssertions(ac)
    
    sessionPath <- paste0(session, ".Rds")
    checkmate::assertPathForOutput(sessionPath, overwrite = TRUE, .var.name = "session")
    
    if (clearSession && file.exists(sessionPath))
        file.remove(sessionPath)
    
    anaInfo <- analysisInfo(fGroups)
    gNames <- names(fGroups)
    gCount <- length(fGroups)
    fTable <- featureTable(fGroups)
    ftind <- groupFeatIndex(fGroups)
    
    # UNDONE: make topMost/onlyPresent optional/interactive
    EICsTopMost <- getEICsForFGroups(fGroups, rtWindow, 0.001, topMost = 1, FALSE, onlyPresent = TRUE)
    EICsTopMostRG <- EICsAll <- NULL
    
    # format is in [[ana]][[fGroup]], since we only took top most intensive we can throw away the ana dimension
    EICPreviews <- Reduce(modifyList, EICsTopMost)
    EICPreviews <- Map(EICPreviews, names(EICPreviews), f = function(eic, grp)
    {
        anai <- which.max(fGroups[[grp]])
        return(eic[numGTE(eic$time, fTable[[anai]]$retmin[ftind[[grp]][anai]]) &
                   numLTE(eic$time, fTable[[anai]]$retmax[ftind[[grp]][anai]]), ])
    })
    
    hotOpts <- list(rowHeaderWidth = 40, readOnly = TRUE, disableVisualSelection = "area",
                    columnSorting = TRUE, sortIndicator = TRUE, selectCallback = TRUE,
                    currentRowClassName = "currentRow", stretchH = "all",
                    preventOverflow = "horizontal", multiSelect = FALSE,
                    outsideClickDeselects = FALSE, manualColumnResize = TRUE,
                    rowHeaders = NULL)
    
    settings <- getUISettings()
    
    curSession <- NULL
    sessionChanged <- FALSE
    if (file.exists(sessionPath))
        curSession <- readRDS(sessionPath)
    else
    {
        ef <- data.frame(analysis = anaInfo$analysis)
        ef[, gNames] <- TRUE
        curSession <- list(enabledFGroups = gNames, enabledFeatures = ef)
    }
    
    server <- function(input, output, session)
    {
        rValues <- reactiveValues(currentFGroup = gNames[1],
                                  triggerFGroupHotUpdate = 0,
                                  triggerFeatHotUpdate = 0,
                                  enabledFGroups = curSession$enabledFGroups,
                                  # NOTE: should be data.frame not data.table, as Shiny doesn't register changes with the latter
                                  enabledFeatures = curSession$enabledFeatures,
                                  settings = settings,
                                  fGroupPlotMode = "topMost")
        
        getCurFGroupRow <- function()
        {
            tbl <- rhandsontable::hot_to_r(input$fGroupsHot)
            tblRow <- match(rValues$currentFGroup, tbl$group)
            if (is.na(tblRow))
                warning("Cannot find fgroup row!")
            return(tblRow)
        }
        
        updateFGroupRow <- function(new)
        {
            old <- rValues$currentFGroup
            rValues$currentFGroup <- new
            
            # update feature selection if needed
            if (!isTRUE(all.equal(rValues$enabledFeatures[[old]],
                                  rValues$enabledFeatures[[new]])))
                rValues$triggerFeatHotUpdate <- rValues$triggerFeatHotUpdate + 1
        }
        
        advanceFGroupSelection <- function(dir)
        {
            tbl <- rhandsontable::hot_to_r(input$fGroupsHot)
            tblRow <- getCurFGroupRow()
            if (!is.na(tblRow))
            {
                tblRow <- tblRow + dir
                if (tblRow < 1)
                    tblRow <- 1
                else if (tblRow > nrow(tbl))
                    tblRow <- nrow(tbl)
                
                updateFGroupRow(tbl$group[tblRow])
                session$sendCustomMessage("selectFGroupRow", tblRow)
            }
        }
        
        toggleFGroup <- function()
        {
            tblRow <- getCurFGroupRow()
            if (!is.na(tblRow))
                session$sendCustomMessage("toggleFGroupRow", tblRow)
        }
        
        reAddHOT <- function(name)
        {
            # BUG: avoid errors after adding/removing columns when column sorting is enabled.
            # work-around from https://github.com/jrowen/rhandsontable/issues/303
            removeUI(selector = paste0("#", name))
            insertUI(selector = paste0("#", name, "Container"), where = "afterEnd", ui = rhandsontable::rHandsontableOutput(name))
        }
        
        getInputSettings <- function()
        {
            return(sapply(names(rValues$settings), function(col) input[[col]], simplify = FALSE))
        }
        
        doApplySettings <- function(settings)
        {
            curSettings <- rValues$settings
            rValues$settings <- settings
            if (!isTRUE(all.equal(curSettings$fGroupColumns, settings$fGroupColumns)) ||
                curSettings$featQuantity != settings$featQuantity ||
                curSettings$fGroupQuantity != settings$fGroupQuantity)
                reAddHOT("fGroupsHot")
            if (!isTRUE(all.equal(curSettings$featureColumns, settings$featureColumns)))
                reAddHOT("featuresHot")
            saveUISettings(settings)
        }
        
        syncInputSettings <- function()
        {
            for (s in c("retUnit", "featQuantity", "fGroupQuantity"))
                updateRadioButtons(session, s, selected = rValues$settings[[s]])
            for (s in c("fGroupColumns", "featureColumns"))
                updateCheckboxGroupInput(session, s, selected = rValues$settings[[s]])
        }
        
        setSessionChanged <- function(changed)
        {
            if (changed != sessionChanged)
            {
                sessionChanged <<- changed
                shinyjs::toggleState("saveSession", condition = changed)
                session$sendCustomMessage("setSessionChanged", changed)
            }
        }
        
        fGroupData <- reactive({
            getScores <- any(c("totalScore", "otherScores") %in% rValues$settings$fGroupColumns)
            args <- list(fGroups, areas = rValues$settings$featQuantity == "area",
                         average = rValues$settings$fGroupQuantity == "average",
                         qualities = if (getScores) "score" else FALSE)
            if (isScreening(fGroups) && any(c("suspProp", "estIDLevel", "suspOther") %in% rValues$settings$fGroupColumns))
                args <- c(args, list(collapseSuspects = NULL))
            gData <- do.call(as.data.table, args)
            
            if (rValues$settings$fGroupQuantity == "max")
            {
                gData[, max_intensity := rowSums(.SD) / nrow(anaInfo), .SDcols = anaInfo$analysis]
                gData <- gData[, (anaInfo$analysis) := NULL]
            }
            
            if ("EICPreview" %in% rValues$settings$fGroupColumns)
            {
                gData[, EIC := sapply(gNames, function(g)
                {
                    jsonlite::toJSON(list(values = EICPreviews[[g]]$intensity, xvalues = EICPreviews[[g]]$time,
                                          options = list(type = "line", height = 50)))
                })]
                setcolorder(gData, c("group", "EIC"))
            }

            if (!"retMZ" %in% rValues$settings$fGroupColumns)
                gData[, c("ret", "mz") := NULL]
            else if (rValues$settings$retUnit == "min")
                gData[, ret := ret / 60]
            
            if (isScreening(fGroups))
            {
                if (!"suspProp" %in% rValues$settings$fGroupColumns)
                    gData[, (intersect(names(gData), c("susp_rt", "susp_mz", "fragments_mz", "fragments_formula"))) := NULL]
                if (!"estIDLevel" %in% rValues$settings$fGroupColumns && !is.null(gData[["estIDLevel"]]))
                    gData[, "estIDLevel" := NULL]
                if (!"suspOther" %in% rValues$settings$fGroupColumns)
                    gData[, (intersect(names(gData),
                                       c("suspFormRank", "suspCompRank", "annSimForm", "annSimComp", "annSimBoth",
                                         "maxFrags", "maxFragMatches", "maxFragMatchesRel"))) := NULL]
            }
            
            if (getScores && hasFGroupScores(fGroups))
            {
                # scores were added by as.data.table(). Remove those we don't want.
                if (!"totalScore" %in% rValues$settings$fGroupColumns)
                    gData[, totalScore := NULL]
                if (!"otherScores" %in% rValues$settings$fGroupColumns)
                    gData[, (c(featureScoreNames(), featureGroupScoreNames())) := NULL]
            }
            
            return(gData)
        })
        
        featureData <- reactive({
            fti <- ftind[[rValues$currentFGroup]]
            ft <- fTable[fti != 0]; ai <- anaInfo[fti != 0, ]; fti <- fti[fti != 0]
            feat <- rbindlist(Map(ft, fti, f = function(f, i) f[i]))
            
            divret <- if (rValues$settings$retUnit == "min") 60 else 1
            
            fData <- data.table(analysis = ai$analysis)
            if ("retMZ" %in% rValues$settings$featureColumns)
                fData[, c("ret", "mz") := .(feat$ret / divret, feat$mz)]
            if ("rGroup" %in% rValues$settings$featureColumns)
                fData[, replicate_group := ai$group]
            if ("blank" %in% rValues$settings$featureColumns)
                fData[, blank := ai$blank]
            if ("quantity" %in% rValues$settings$featureColumns)
                fData[, quantity := if (rValues$settings$featQuantity == "intensity") feat$intensity else feat$area]
            if ("rtMZRange" %in% rValues$settings$featureColumns)
                fData[, c("retmin", "retmax", "mzmin", "mzmax") := .(feat$retmin / divret, feat$retmax / divret,
                                                                     feat$mzmin, feat$mzmax)]
            if (hasFGroupScores(fGroups))
            {
                if ("otherScores" %in% rValues$settings$featureColumns)
                    fData[, (featureScoreNames()) := feat[, featureScoreNames(), with = FALSE]]
                if ("totalScore" %in% rValues$settings$featureColumns)
                    fData[, totalScore := feat$totalScore]
            }

            return(fData)
        })
        
        observeEvent(input$keys, {
            if (input$tabs == "Feature groups")
            {
                switch(input$keys,
                       p = advanceFGroupSelection(-1),
                       n = advanceFGroupSelection(1),
                       t = toggleFGroup())
            }
        })
        
        observeEvent(input$saveSession, {
            saveRDS(list(enabledFGroups = rValues$enabledFGroups,
                         enabledFeatures = rValues$enabledFeatures,
                         version = 1), sessionPath)
            setSessionChanged(FALSE)
        })
        
        observeEvent(input$tabs, {
            if (input$tabs != "Settings")
            {
                set <- getInputSettings()
                if (!isTRUE(all.equal(set, rValues$settings)))
                {
                    showModal(modalDialog("Settings were changed but not yet applied. Do you want to do so now?",
                                          footer = tagList(actionButton("saveAndApplyModal", "Save & Apply", icon = icon("save")),
                                                           actionButton("discardModal", "Discard changes", icon = icon("window-close")))))
                }
                
                if (input$tabs == "Feature groups")
                    rValues$fGroupPlotMode <- input$fGroupPlotMode
                else
                    rValues$fGroupPlotMode <- "all"
            }
        })
        
        observeEvent(input$saveAndApplyModal, {
            doApplySettings(getInputSettings())
            shinyjs::disable("saveApplySettings")
            removeModal()
        })

        observeEvent(input$discardModal, {
            syncInputSettings()
            shinyjs::disable("saveApplySettings")
            removeModal()
        })
        
        observeEvent(input$previousGroup, {
            advanceFGroupSelection(-1)
        })
        
        observeEvent(input$nextGroup, {
            advanceFGroupSelection(1)
        })
        
        observeEvent(input$toggleGroup, {
            toggleFGroup()
        })
        
        observeEvent(input$fGroupPlotMode, {
            if ((input$fGroupPlotMode == "topMostByRGroup" && is.null(EICsTopMostRG)) ||
                (input$fGroupPlotMode == "all" && is.null(EICsAll)))
            {
                not <- showNotification("Loading EICs...", duration = NULL, closeButton = FALSE, type = "message")
                if (input$fGroupPlotMode == "topMostByRGroup")
                    EICsTopMostRG <<- getEICsForFGroups(fGroups, rtWindow, 0.001, 1, TRUE, TRUE)
                else
                    EICsAll <<- getEICsForFGroups(fGroups, rtWindow, 0.001, NULL, FALSE, TRUE)
                removeNotification(not)
            }
            rValues$fGroupPlotMode <- input$fGroupPlotMode
        })
        
        observeEvent(input$fGroupsHot, {
            # HACK: input$fGroupsHot$params$maxRows: make sure we don't have empty table as hot_to_r errors otherwise
            if (input$fGroupsHot$params$maxRows > 0)
            {
                tbl <- rhandsontable::hot_to_r(input$fGroupsHot)
                keep <- tbl[keep == TRUE]$group
                notkeep <- tbl[keep == FALSE]$group
                oldeg <- rValues$enabledFGroups
                rValues$enabledFGroups <- setdiff(union(rValues$enabledFGroups, keep), notkeep)
                egChanged <- !isTRUE(all.equal(oldeg, rValues$enabledFGroups))
                if (egChanged)
                    setSessionChanged(TRUE)
                
                selr <- input$fGroupsHot_select$select$rAll[1]
                
                # HACK: init selection
                if (is.null(selr))
                    session$sendCustomMessage("selectFGroupRow", 1)
                
                # filters are active?
                if (length(input$showWhat) < 2)
                {
                    # update table
                    if (egChanged)
                        rValues$triggerFGroupHotUpdate <- rValues$triggerFGroupHotUpdate + 1
                    else
                    {
                        # update selection after table update was triggered
                        if (!is.null(selr))
                            updateFGroupRow(tbl$group[selr])
                    }
                }
            }
        })

        observeEvent(input$fGroupsHot_select$select$r, {
            tbl <- rhandsontable::hot_to_r(input$fGroupsHot)
            updateFGroupRow(tbl$group[input$fGroupsHot_select$select$rAll[1]])
        })
        
        observeEvent(input$enableAllGroups, {
            rValues$enabledFGroups <- gNames
            rValues$triggerFGroupHotUpdate <- rValues$triggerFGroupHotUpdate + 1
            setSessionChanged(TRUE)
        })
        observeEvent(input$disableAllGroups, {
            rValues$enabledFGroups <- character()
            rValues$triggerFGroupHotUpdate <- rValues$triggerFGroupHotUpdate + 1
            setSessionChanged(TRUE)
        })
        
        observeEvent(input$resetFeatures, {
            rValues$enabledFeatures[, gNames] <-TRUE
            rValues$triggerFeatHotUpdate <- rValues$triggerFeatHotUpdate + 1
            setSessionChanged(TRUE)
        })
        
        observeEvent(input$featuresHot, {
            # HACK: input$featuresHot$params$maxRows: make sure we don't have empty table as hot_to_r errors otherwise
            if (input$featuresHot$params$maxRows > 0)
            {
                tbl <- rhandsontable::hot_to_r(input$featuresHot)
                oldef <- rValues$enabledFeatures[[rValues$currentFGroup]]
                rValues$enabledFeatures[match(tbl$analysis, rValues$enabledFeatures$analysis),
                                        rValues$currentFGroup] <- tbl$keep
                if (!isTRUE(all.equal(oldef, rValues$enabledFeatures[[rValues$currentFGroup]])))
                    setSessionChanged(TRUE)
            }
        })
        
        observeEvent(input$enableAllFeatures, {
            rValues$enabledFeatures[, rValues$currentFGroup] <- TRUE
            rValues$triggerFeatHotUpdate <- rValues$triggerFeatHotUpdate + 1
            setSessionChanged(TRUE)
        })
        observeEvent(input$disableAllFeatures, {
            rValues$enabledFeatures[, rValues$currentFGroup] <- FALSE
            rValues$triggerFeatHotUpdate <- rValues$triggerFeatHotUpdate + 1
            setSessionChanged(TRUE)
        })
        
        observeEvent({ input$retUnit; input$featQuantity; input$fGroupQuantity; input$fGroupColumns; input$featureColumns }, {
            shinyjs::toggleState("saveApplySettings", !isTRUE(all.equal(getInputSettings(), rValues$settings)))
        })
        
        observeEvent(input$saveApplySettings, {
            doApplySettings(getInputSettings())
            shinyjs::disable("saveApplySettings")
        })
        
        observeEvent(input$resetSettings, {
            doApplySettings(getDefaultUISettings())
            syncInputSettings()
            shinyjs::disable("saveApplySettings")
        })
        
        output$pageTitle <- renderText({
            sprintf("Group %s (%d/%d)", rValues$currentFGroup,
                    match(rValues$currentFGroup, gNames), gCount)
        })
        
        output$plotChrom <- renderPlot({
            EICs <- switch(rValues$fGroupPlotMode,
                topMost = EICsTopMost,
                topMostByRGroup = EICsTopMostRG,
                all = EICsAll
            )
            
            fg <- fGroups[, rValues$currentFGroup]
            if (rValues$fGroupPlotMode == "all") # UNDONE: also for rGroups top most somehow?
                fg <- fg[rValues$enabledFeatures[[rValues$currentFGroup]]]
            
            withr::with_par(list(mar = c(4, 4, 0.1, 1), cex = 1.5), {
                plotChroms(fg, EICs = EICs, colourBy = "rGroups", showPeakArea = TRUE,
                           showFGroupRect = FALSE, title = "",
                           topMost = if (rValues$fGroupPlotMode == "all") NULL else 1,
                           topMostByRGroup = rValues$fGroupPlotMode == "topMostByRGroup",
                           retMin = rValues$settings$retUnit == "min")
            })
        })
        
        output$fGroupsHot <- rhandsontable::renderRHandsontable({
            not <- showNotification("Updating feature group table...", duration = NULL, closeButton = FALSE, type = "message")
            rValues$triggerFGroupHotUpdate
            
            gData <- fGroupData()
            gData[, keep := group %in% isolate(rValues$enabledFGroups)]
            setcolorder(gData, c("group", "keep"))
            
            if (!"Keep" %in% input$showWhat)
                gData <- gData[keep == FALSE, ]
            if (!"Don't keep" %in% input$showWhat)
                gData <- gData[keep == TRUE, ]
            
            hot <- do.call(rhandsontable::rhandsontable,
                           c(list(gData, width = NULL, height = 200, maxRows = nrow(gData)),
                             hotOpts)) %>%
                # rhandsontable::hot_cols(valign = "htMiddle", fixedColumnsLeft = 2) %>%
                rhandsontable::hot_cols(valign = "htMiddle", fixedColumnsLeft = 2) %>%
                rhandsontable::hot_col("keep", readOnly = FALSE, halign = "htCenter") %>%
                rhandsontable::hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE, customOpts = list(
                    enableAll = list(
                        name = "Enable all",
                        callback = htmlwidgets::JS(
                            'function(key, options) { Shiny.onInputChange("enableAllGroups", Math.random()); }'
                        )
                    ),
                    disableAll = list(
                        name = "Disable all",
                        callback = htmlwidgets::JS(
                            'function(key, options) { Shiny.onInputChange("disableAllGroups", Math.random()); }'
                        )
                    )
                )) %>%
                # BUG: table is messed up after tab switch
                # work-around from https://github.com/jrowen/rhandsontable/issues/366
                htmlwidgets::onRender("function(el, x){
                  var hot = this.hot;
                  $('a[data-value=\"Feature groups\"').on('click', function(){
                    setTimeout(function(){hot.render();}, 0);
                  });
                }")
            
            if (!is.null(gData[["EIC"]]))
            {
                hot <- rhandsontable::hot_col(hot, "EIC", renderer = htmlwidgets::JS("renderSparkline")) %>%
                    rhandsontable::hot_rows(rowHeights = 50)
            }
            
            # HACK: disable some default context options
            hot$x$contextMenu$items[c("undo", "redo", "alignment")] <- NULL
            
            # BUG: table is messed up after tab switch
            # work-around from https://github.com/jrowen/rhandsontable/issues/366
            outputOptions(output, "fGroupsHot", suspendWhenHidden = FALSE)
            
            removeNotification(not)
            return(hot)
        })
        
        output$featuresHot <- rhandsontable::renderRHandsontable({
            rValues$triggerFeatHotUpdate
            
            fData <- featureData()
            isolate(fData[, keep := rValues$enabledFeatures[match(analysis, rValues$enabledFeatures$analysis),
                          rValues$currentFGroup]])
            setcolorder(fData, c("analysis", "keep"))
            
            hot <- do.call(rhandsontable::rhandsontable,
                           c(list(fData, height = 200, maxRows = nrow(fData)), hotOpts)) %>%
                rhandsontable::hot_col("keep", readOnly = FALSE, halign = "htCenter") %>%
                rhandsontable::hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE,
                                                customOpts = list(
                    enableAll = list(
                        name = "Enable all",
                        callback = htmlwidgets::JS(
                            'function(key, options) { Shiny.onInputChange("enableAllFeatures", Math.random()); }'
                        )
                    ),
                    disableAll = list(
                        name = "Disable all",
                        callback = htmlwidgets::JS(
                            'function(key, options) { Shiny.onInputChange("disableAllFeatures", Math.random()); }'
                        )
                    )
                )) %>%
                # BUG: table is messed up after tab switch
                # work-around from https://github.com/jrowen/rhandsontable/issues/366
                htmlwidgets::onRender("function(el, x){
                  var hot = this.hot;
                  $('a[data-value=\"Features\"').on('click', function(){
                    setTimeout(function(){hot.render();}, 0);
                  });
                }")
            
            # HACK: disable some default context options
            hot$x$contextMenu$items[c("undo", "redo", "alignment")] <- NULL

            # BUG: table is messed up after tab switch
            # work-around from https://github.com/jrowen/rhandsontable/issues/366
            outputOptions(output, "featuresHot", suspendWhenHidden = FALSE)
            
            return(hot)
        })
    }
    
    runApp(shinyApp(getCheckFeatsUI(settings), server))
})
