#' @include main.R

checkUIInterface <- setRefClass("checkUIInterface", contains = "VIRTUAL",
                                fields = list(primarySelections = "character", curSession = "list"))

checkUIInterface$methods(

    resetSecondaryUITitle = function() stop("VIRTUAL"),
    settingsTabUI = function() stop("VIRTUAL"),
    
    primaryTab = function() stop("VIRTUAL"),
    secondaryTab = function() stop("VIRTUAL"),
    UISettingsFileName = function() stop("VIRTUAL"),
    
    primarySettingsChanged = function(cur, new) stop("VIRTUAL"),
    secondarySettingsChanged = function(cur, new) stop("VIRTUAL"),
    syncInputSettings = function(session, settings) stop("VIRTUAL"),
    
    primaryTableData = function(rValues) stop("VIRTUAL"),
    secondaryTableData = function(rValues) stop("VIRTUAL"),
    
    doObserveEvents = function(input) NULL,
    
    plotChrom = function(rValues) stop("VIRTUAL")
)

getUISettingsPath <- function(fileName)
{
    dirPath <- RUserDir("patRoon", "config")
    mkdirp(dirPath)
    return(file.path(dirPath, fileName))
}

getUISettings <- function(fileName)
{
    path <- getUISettingsPath(fileName)
    if (!file.exists(path))
    {
        ret <- getDefaultUISettings()
        saveUISettings(ret)
    }
    else
        ret <- yaml::read_yaml(path, eval.expr = FALSE)
    return(ret)
}

getDefaultUISettings <- function()
{
    return(list(retUnit = "sec", featQuantity = "intensity", fGroupQuantity = "average",
                fGroupColumns = c("retMZ", "estIDLevel", "totalScore"),
                featureColumns = c("retMZ", "quantity", "totalScore")))
}

saveUISettings <- function(fileName, settings)
{
    yaml::write_yaml(settings, getUISettingsPath(fileName), indent = 4)
}

# UNDONE
if (FALSE) {
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
}

getCheckUI <- function(title, primaryTab, secondaryTab, resetSecondaryTitle, settingsTab)
{
    showOpts <- c("Keep", "Don't keep")
    
    fillPage(
        tags$head(includeScript(system.file("js", "utils-checkUI.js", package = "patRoon"))),
        shinyjs::useShinyjs(),
        
        keys::useKeys(),
        keys::keysInput("keys", c("p", "n", "t")),
        
        title = title,
        
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
                
                tabsetPanel(id = "tabs", tabPanel(primaryTab, fillCol(
                    flex = c(NA, 1),
                    fillRow(
                        height = 40,
                        flex = c(NA, NA, NA, NA, NA, 1),
                        
                        fillCol(
                            width = 45,
                            actionButton("previousRow", "", icon("arrow-left"))
                        ),
                        fillCol(
                            width = 45,
                            actionButton("nextRow", "", icon("arrow-right"))
                        ),
                        fillCol(
                            width = 100,
                            actionButton("toggleRow", "Toggle", icon("toggle-on"))
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
                            # UNDONE: disable for components somehow
                            selectInput("fGroupPlotMode", NULL,
                                        c("Top most group" = "topMost", "Top most replicates" = "topMostByRGroup",
                                          "All" = "all"), "topMost")
                        ),
                        div(style = "margin: 8px 10px 12px; font-size: small; text-align: right;",
                            HTML("<b>n</b>: next; <b>p</b>: previous; <b>t</b>: toggle"))
                    ),
                    fillRow(
                        tags$div(id = "primaryHotContainer", rhandsontable::rHandsontableOutput("primaryHot"))
                    )
                )),
                tabPanel(secondaryTab, fillCol(
                    flex = c(NA, 1),
                    
                    fillRow(
                        height = 40,
                        actionButton("resetAllSecondary", resetSecondaryTitle, icon("check-square"))
                    ),
                    
                    fillRow(
                        tags$div(id = "secondaryHotContainer", rhandsontable::rHandsontableOutput("secondaryHot"))
                    )
                )),
                tabPanel("Settings", fillCol(
                    flex = NA,
                    
                    settingsTab,
                    
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

runCheckUI <- function(UIInterface)
{
    hotOpts <- list(rowHeaderWidth = 40, readOnly = TRUE, disableVisualSelection = "area",
                    columnSorting = TRUE, sortIndicator = TRUE, selectCallback = TRUE,
                    currentRowClassName = "currentRow", stretchH = "all",
                    preventOverflow = "horizontal", multiSelect = FALSE,
                    outsideClickDeselects = FALSE, manualColumnResize = TRUE,
                    rowHeaders = NULL)
    
    settings <- getUISettings(UIInterface$UISettingsFileName())
    
    hotRenderJSFunc <- function(tab) sprintf("function(el, x) {
                  var hot = this.hot;
                  $('a[data-value=\"%s\"').on('click', function(){
                      setTimeout(function() { hot.render(); }, 0);
                  });
                }", tab)
    
    primarySelections <- UIInterface$primarySelections
    curSession <- UIInterface$curSession
    
    server <- function(input, output, session)
    {
        rValues <- reactiveValues(currentPrimSel = primarySelections[1],
                                  triggerPrimaryHotUpdate = 0,
                                  triggerSecondaryHotUpdate = 0,
                                  primarySelections = curSession$primarySelections,
                                  # NOTE: should be data.frame not data.table, as Shiny doesn't register changes with the latter
                                  secondarySelections = curSession$secondarySelections,
                                  settings = settings)
                                  # UNDONE fGroupPlotMode = "topMost")
        
        getCurPrimaryRow <- function()
        {
            tbl <- rhandsontable::hot_to_r(input$primaryHot)
            tblRow <- match(rValues$currentPrimSel, tbl$name)
            if (is.na(tblRow))
                warning("Cannot find fgroup row!")
            return(tblRow)
        }
        
        updatePrimaryRow <- function(new)
        {
            old <- rValues$currentPrimSel
            rValues$currentPrimSel <- new
            
            # update feature selection if needed
            if (!isTRUE(all.equal(rValues$secondarySelections[[old]],
                                  rValues$secondarySelections[[new]])))
                rValues$triggerSecondaryHotUpdate <- rValues$triggerSecondaryHotUpdate + 1
        }
        
        advancePrimaryRow <- function(dir)
        {
            tbl <- rhandsontable::hot_to_r(input$primaryHot)
            tblRow <- getCurPrimaryRow()
            if (!is.na(tblRow))
            {
                tblRow <- tblRow + dir
                if (tblRow < 1)
                    tblRow <- 1
                else if (tblRow > nrow(tbl))
                    tblRow <- nrow(tbl)
                
                updatePrimaryRow(tbl$group[tblRow])
                session$sendCustomMessage("selectPrimaryRow", tblRow)
            }
        }
        
        togglePrimaryRow <- function()
        {
            tblRow <- getCurPrimaryRow()
            if (!is.na(tblRow))
                session$sendCustomMessage("togglePrimaryRow", tblRow)
        }
        
        reAddHOT <- function(name)
        {
            # BUG: avoid errors after adding/removing columns when column sorting is enabled.
            # work-around from https://github.com/jrowen/rhandsontable/issues/303
            removeUI(selector = paste0("#", name))
            insertUI(selector = paste0("#", name, "Container"), where = "afterEnd",
                     ui = rhandsontable::rHandsontableOutput(name))
        }
        
        getInputSettings <- function()
        {
            return(sapply(names(rValues$settings), function(col) input[[col]], simplify = FALSE))
        }
        
        doApplySettings <- function(settings)
        {
            curSettings <- rValues$settings
            rValues$settings <- settings
            
            if (UIInterface$primarySettingsChanged(curSettings, settings))
                reAddHOT("primaryHot")
            if (UIInterface$secondarySettingsChanged(curSettings, settings))
                reAddHOT("secondaryHot")
            
            saveUISettings(settings)
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
        
        primaryTable <- reactive({ UIInterface$primaryTableData(rValues) })
        secondaryTable <- reactive({ UIInterface$secondaryTableData(rValues) })
        
        observeEvent(input$keys, {
            if (input$tabs == "Feature groups")
            {
                switch(input$keys,
                       p = advancePrimaryRow(-1),
                       n = advancePrimaryRow(1),
                       t = togglePrimaryRow())
            }
        })
        
        observeEvent(input$saveSession, {
            saveRDS(list(primarySelections = rValues$primarySelections,
                         secondarySelections = rValues$secondarySelections), sessionPath)
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
            }
        })
        
        observeEvent(input$saveAndApplyModal, {
            doApplySettings(getInputSettings())
            shinyjs::disable("saveApplySettings")
            removeModal()
        })
        
        observeEvent(input$discardModal, {
            UIInterface$syncInputSettings(session, rValues$settings)
            shinyjs::disable("saveApplySettings")
            removeModal()
        })
        
        observeEvent(input$previousRow, {
            advancePrimaryRow(-1)
        })
        
        observeEvent(input$nextRow, {
            advancePrimaryRow(1)
        })
        
        observeEvent(input$toggleRow, {
            togglePrimaryRow()
        })
        
        observeEvent(input$primaryHot, {
            # HACK: input$primaryHot$params$maxRows: make sure we don't have empty table as hot_to_r errors otherwise
            if (input$primaryHot$params$maxRows > 0)
            {
                tbl <- rhandsontable::hot_to_r(input$primaryHot)
                keep <- tbl[keep == TRUE]$name
                notkeep <- tbl[keep == FALSE]$name
                oldsel <- rValues$primarySelections
                rValues$primarySelections <- setdiff(union(rValues$primarySelections, keep), notkeep)
                selChanged <- !isTRUE(all.equal(oldsel, rValues$primarySelections))
                if (selChanged)
                    setSessionChanged(TRUE)
                
                selr <- input$primaryHot_select$select$rAll[1]
                
                # HACK: init selection
                if (is.null(selr))
                    session$sendCustomMessage("selectPrimaryRow", 1)
                
                # filters are active?
                if (length(input$showWhat) < 2)
                {
                    # update table
                    if (selChanged)
                        rValues$triggerPrimaryHotUpdate <- rValues$triggerPrimaryHotUpdate + 1
                    else
                    {
                        # update selection after table update was triggered
                        if (!is.null(selr))
                            updatePrimaryRow(tbl$group[selr])
                    }
                }
            }
        })
        
        observeEvent(input$primaryHot_select$select$r, {
            tbl <- rhandsontable::hot_to_r(input$primaryHot)
            updatePrimaryRow(tbl$group[input$primaryHot_select$select$rAll[1]])
        })
        
        observeEvent(input$enableAllPrimary, {
            rValues$primarySelections <- primarySelections
            rValues$triggerPrimaryHotUpdate <- rValues$triggerPrimaryHotUpdate + 1
            setSessionChanged(TRUE)
        })
        observeEvent(input$disableAllPrimary, {
            rValues$primarySelections <- character()
            rValues$triggerPrimaryHotUpdate <- rValues$triggerPrimaryHotUpdate + 1
            setSessionChanged(TRUE)
        })
        
        observeEvent(input$resetAllSecondary, {
            rValues$secondarySelections[, primarySelections] <-TRUE
            rValues$triggerSecondaryHotUpdate <- rValues$triggerSecondaryHotUpdate + 1
            setSessionChanged(TRUE)
        })
        
        observeEvent(input$secondaryHot, {
            # HACK: input$secondaryHot$params$maxRows: make sure we don't have empty table as hot_to_r errors otherwise
            if (input$secondaryHot$params$maxRows > 0)
            {
                tbl <- rhandsontable::hot_to_r(input$secondaryHot)
                oldsel <- rValues$secondarySelections[[rValues$currentPrimSel]]
                rValues$secondarySelections[match(tbl$name, rValues$secondarySelections$name),
                                            rValues$currentPrimSel] <- tbl$keep
                if (!isTRUE(all.equal(oldsel, rValues$secondarySelections[[rValues$currentPrimSel]])))
                    setSessionChanged(TRUE)
            }
        })
        
        observeEvent(input$enableAllSecondary, {
            rValues$secondarySelections[, rValues$currentPrimSel] <- TRUE
            rValues$triggerSecondaryHotUpdate <- rValues$triggerSecondaryHotUpdate + 1
            setSessionChanged(TRUE)
        })
        observeEvent(input$disableAllSecondary, {
            rValues$secondarySelections[, rValues$currentPrimSel] <- FALSE
            rValues$triggerSecondaryHotUpdate <- rValues$triggerSecondaryHotUpdate + 1
            setSessionChanged(TRUE)
        })
        
        # UNDONE
        if (FALSE){
        observeEvent({ input$retUnit; input$featQuantity; input$fGroupQuantity; input$fGroupColumns; input$featureColumns }, {
            shinyjs::toggleState("saveApplySettings", !isTRUE(all.equal(getInputSettings(), rValues$settings)))
        })}
        
        observeEvent(input$saveApplySettings, {
            doApplySettings(getInputSettings())
        })
        
        observeEvent(input$resetSettings, {
            doApplySettings(getDefaultUISettings())
            UIInterface$syncInputSettings(session, rValues$settings)
            shinyjs::disable("saveApplySettings")
        })

        UIInterface$doObserveEvents(input)
        
        output$pageTitle <- renderText({
            sprintf("Selection %s (%d/%d)", rValues$currentPrimSel,
                    match(rValues$currentPrimSel, primarySelections), length(primarySelections))
        })
        
        output$plotChrom <- renderPlot({ UIInterface$plotChrom(rValues) })

        output$primaryHot <- rhandsontable::renderRHandsontable({
            not <- showNotification(sprintf("Updating %s table...", primaryTabName), duration = NULL,
                                    closeButton = FALSE, type = "message")
            rValues$triggerPrimaryHotUpdate
            
            tabData <- primaryTable()
            tabData[, keep := group %in% isolate(rValues$primarySelections)]
            setcolorder(tabData, c("name", "keep"))
            
            if (!"Keep" %in% input$showWhat)
                tabData <- tabData[keep == FALSE, ]
            if (!"Don't keep" %in% input$showWhat)
                tabData <- tabData[keep == TRUE, ]
            
            hot <- do.call(rhandsontable::rhandsontable,
                           c(list(tabData, width = NULL, height = 200, maxRows = nrow(tabData)),
                             hotOpts)) %>%
                # rhandsontable::hot_cols(valign = "htMiddle", fixedColumnsLeft = 2) %>%
                rhandsontable::hot_cols(valign = "htMiddle", fixedColumnsLeft = 2) %>%
                rhandsontable::hot_col("keep", readOnly = FALSE, halign = "htCenter") %>%
                rhandsontable::hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE, customOpts = list(
                    enableAll = list(
                        name = "Enable all",
                        callback = htmlwidgets::JS(
                            'function(key, options) { Shiny.onInputChange("enableAllPrimary", Math.random()); }'
                        )
                    ),
                    disableAll = list(
                        name = "Disable all",
                        callback = htmlwidgets::JS(
                            'function(key, options) { Shiny.onInputChange("disableAllPrimary", Math.random()); }'
                        )
                    )
                )) %>%
                # BUG: table is messed up after tab switch
                # work-around from https://github.com/jrowen/rhandsontable/issues/366
                htmlwidgets::onRender(hotRenderJSFunc(primaryTabName))
            
            if (!is.null(tabData[["EIC"]]))
            {
                hot <- rhandsontable::hot_col(hot, "EIC", renderer = htmlwidgets::JS("renderSparkline")) %>%
                    rhandsontable::hot_rows(rowHeights = 50)
            }
            
            # HACK: disable some default context options
            hot$x$contextMenu$items[c("undo", "redo", "alignment")] <- NULL
            
            # BUG: table is messed up after tab switch
            # work-around from https://github.com/jrowen/rhandsontable/issues/366
            outputOptions(output, "primaryHot", suspendWhenHidden = FALSE)
            
            removeNotification(not)
            return(hot)
        })
        
        output$secondaryHot <- rhandsontable::renderRHandsontable({
            rValues$triggerSecondaryHotUpdate
            
            tabData <- secondaryTable()
            isolate(tabData[, keep := rValues$secondarySelections[match(analysis, rValues$secondarySelections$analysis),
                                                                  rValues$currentPrimSel]])
            setcolorder(tabData, c("name", "keep"))
            
            hot <- do.call(rhandsontable::rhandsontable,
                           c(list(tabData, height = 200, maxRows = nrow(tabData)), hotOpts)) %>%
                rhandsontable::hot_col("keep", readOnly = FALSE, halign = "htCenter") %>%
                rhandsontable::hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE, customOpts = list(
                    enableAll = list(
                        name = "Enable all",
                        callback = htmlwidgets::JS(
                            'function(key, options) { Shiny.onInputChange("enableAllSecondary", Math.random()); }'
                        )
                    ),
                    disableAll = list(
                        name = "Disable all",
                        callback = htmlwidgets::JS(
                            'function(key, options) { Shiny.onInputChange("disableAllSecondary", Math.random()); }'
                        )
                    )
                )) %>%
                # BUG: table is messed up after tab switch
                # work-around from https://github.com/jrowen/rhandsontable/issues/366
                htmlwidgets::onRender(hotRenderJSFunc(secondaryTabName))
            
            # HACK: disable some default context options
            hot$x$contextMenu$items[c("undo", "redo", "alignment")] <- NULL
            
            # BUG: table is messed up after tab switch
            # work-around from https://github.com/jrowen/rhandsontable/issues/366
            outputOptions(output, "secondaryHot", suspendWhenHidden = FALSE)
            
            return(hot)
        })
    }
    
    runApp(shinyApp(getCheckUI(settings, UIInterface$primaryTab(), UIInterface$secondaryTab(),
                               UIInterface$resetSecondaryUITitle(), UIInterface$settingsTabUI()), server))
}
