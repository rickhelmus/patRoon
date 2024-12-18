# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R

checkUIInterface <- setRefClass("checkUIInterface", contains = "VIRTUAL",
                                fields = list(primarySelections = "character", curSession = "list",
                                              session = "character"))

checkUIInterface$methods(

    UITitle = function() stop("VIRTUAL"),
    resetSecondaryUITitle = function() stop("VIRTUAL"),
    primaryTabTopUI = function() "",
    settingsTabUI = function(settings) stop("VIRTUAL"),
    
    primaryTab = function() stop("VIRTUAL"),
    secondaryTab = function() stop("VIRTUAL"),
    
    defaultUISettings = function() stop("VIRTUAL"),
    UISettingsFileName = function() stop("VIRTUAL"),
    
    getSecondarySelections = function(primSel) stop("VIRTUAL"),
    getSecondarySelectionsFromTab = function(tab) stop("VIRTUAL"),
    
    init = function(rValues) rValues,
    
    settingsChangedExpression = function(input) stop("VIRTUAL"),
    primarySettingsChanged = function(cur, new) stop("VIRTUAL"),
    secondarySettingsChanged = function(cur, new) stop("VIRTUAL"),
    syncInputSettings = function(session, settings) stop("VIRTUAL"),
    
    primaryTableData = function(rValues) stop("VIRTUAL"),
    secondaryTableData = function(rValues) stop("VIRTUAL"),
    
    doObserveEvents = function(input, rValues) NULL,
    
    plotMain = function(input, rValues) stop("VIRTUAL"),
    
    saveSession = function(s) stop("VIRTUAL")
)

getUISettingsPath <- function(fileName)
{
    dirPath <- RUserDir("patRoon", "config")
    mkdirp(dirPath)
    return(file.path(dirPath, fileName))
}

getUISettings <- function(fileName, default)
{
    path <- getUISettingsPath(fileName)
    if (!file.exists(path))
    {
        ret <- default
        saveUISettings(fileName, ret)
    }
    else
        ret <- readYAML(path)
    ret <- ret[setdiff(names(ret), "version")]
    return(ret)
}

saveUISettings <- function(fileName, settings)
{
    settings$version <- 1 # just store for now, in case if ever needed in the future
    writeYAML(settings, getUISettingsPath(fileName))
}

readCheckSession <- function(session, type)
{
    ret <- readYAML(session)
    if (is.null(ret[["type"]]))
        ret$type <- "unknown" # for error below
    if (ret$type != type)
        stop("The specified session file type is incorrect: ", ret$type)
    
    if (length(ret$removeFully) == 0)
        ret$removeFully <- character() # YAML converts to empty list
    
    return(ret)
}

saveCheckSession <- function(session, path, fGroups, type)
{
    session$version <- 1
    session$type <- type
    
    gInfo <- groupInfo(fGroups)
    session$featureGroups <- sapply(names(fGroups),
                                    function(grp) list(ret = gInfo[grp, "rts"], mz = gInfo[grp, "mzs"]),
                                    simplify = FALSE)
    
    writeYAML(session, path)
}

importCheckUISessionGroups <- function(oldSession, fGroups, rtWindow, mzWindow)
{
    
    gInfo <- groupInfo(fGroups)
    gInfoDT <- as.data.table(gInfo[, c("rts", "mzs")])
    setnames(gInfoDT, c("ret", "mz")) # equalize column names between old/new tables
    gInfoDT[, group := rownames(gInfo)]
    
    oldGroupTab <- rbindlist(oldSession$featureGroups, idcol = "group")
    
    warnTol <- FALSE
    newGroups <- setNames(lapply(split(oldGroupTab, seq_len(nrow(oldGroupTab))), function(ogtr)
    {
        gi <- gInfoDT[numLTE(abs(ret - ogtr$ret), rtWindow) & numLTE(abs(mz - ogtr$mz), mzWindow)]
        if (nrow(gi) == 0)
        {
            printf("Could not find any matching feature groups for old group %s\n", ogtr$group)
            warnTol <<- TRUE
            return(NULL)
        }
        else if (nrow(gi) > 1)
        {
            printf("Old group %s matched to multiple new groups: %s\n", ogtr$group, paste0(gi$group, collapse = ", "))
            warnTol <<- TRUE
        }
        return(gi)
    }), oldGroupTab$group)
    
    if (warnTol)
        printf(paste("NOTE: You may consider tweaking the retention and/or m/z tolerances by setting",
                     "the rtWindow/mzWindow arguments\n"))
    
    return(rbindlist(pruneList(newGroups), idcol = "oldGroup"))
}

getCheckUI <- function(UIInterface, settings)
{
    showOpts <- c("Keep", "Don't keep")
    
    fillPage(
        tags$head(includeScript(system.file("js", "utils-checkUI.js", package = "patRoon"))),
        shinyjs::useShinyjs(),
        
        keys::useKeys(),
        keys::keysInput("keys", c("p", "n", "t")),
        
        title = UIInterface$UITitle(),
        
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
                    plotOutput("plotMain", width = "100%", height = "100%")
                )
            ),
            
            fillRow(
                height = 300,
                
                tabsetPanel(id = "tabs",
                            tabPanel(UIInterface$primaryTab(), fillCol(
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
                                    UIInterface$primaryTabTopUI(),
                                    div(style = "margin: 8px 10px 12px; font-size: small; text-align: right;",
                                        HTML("<b>n</b>: next; <b>p</b>: previous; <b>t</b>: toggle"))
                                ),
                                fillRow(
                                    tags$div(id = "primaryHotContainer",
                                             rhandsontable::rHandsontableOutput("primaryHot"))
                                )
                            )),
                            tabPanel(UIInterface$secondaryTab(), fillCol(
                                flex = c(NA, 1),
                                
                                fillRow(
                                    height = 40,
                                    actionButton("resetAllSecondary", UIInterface$resetSecondaryUITitle(),
                                                 icon("check-square"))
                                ),
                                
                                fillRow(
                                    tags$div(id = "secondaryHotContainer",
                                             rhandsontable::rHandsontableOutput("secondaryHot"))
                                )
                            )),
                            tabPanel("Settings", fillCol(
                                flex = NA,
                                
                                UIInterface$settingsTabUI(settings),
                                
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
    
    hotRenderJSFunc <- function(tab) sprintf("function(el, x) {
                  var hot = this.hot;
                  $('a[data-value=\"%s\"').on('click', function(){
                      setTimeout(function() { hot.render(); }, 0);
                  });
                }", tab)
    
    settings <- getUISettings(UIInterface$UISettingsFileName(), UIInterface$defaultUISettings())
    sessionChanged <- FALSE
    
    server <- function(input, output, session)
    {
        rValues <- reactiveValues(currentPrimSel = UIInterface$primarySelections[1],
                                  triggerPrimaryHotUpdate = 0,
                                  triggerSecondaryHotUpdate = 0,
                                  removeFully = UIInterface$curSession$removeFully,
                                  # NOTE: should be data.frame not data.table, as Shiny doesn't register changes with the latter
                                  removePartially = UIInterface$curSession$removePartially,
                                  settings = settings)
        
        rValues <- UIInterface$init(rValues)
        
        getCurPrimaryRow <- function()
        {
            tbl <- rhandsontable::hot_to_r(input$primaryHot)
            tblRow <- match(rValues$currentPrimSel, tbl[[1]])
            if (is.na(tblRow))
                warning("Cannot find row!")
            return(tblRow)
        }
        
        updatePrimaryRow <- function(new)
        {
            old <- rValues$currentPrimSel
            rValues$currentPrimSel <- new
            
            # update secondary selection if needed
            rpo <- rValues$removePartially[[old]]; rpn <- rValues$removePartially[[new]]
            if (is.null(rpo) != is.null(rpn) || (!is.null(rpo) && !setequal(rpo, rpn)))
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
                
                updatePrimaryRow(tbl[[1]][tblRow])
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
            
            saveUISettings(UIInterface$UISettingsFileName(), settings)
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
            if (input$tabs == UIInterface$primaryTab())
            {
                switch(input$keys,
                       p = advancePrimaryRow(-1),
                       n = advancePrimaryRow(1),
                       t = togglePrimaryRow())
            }
        })
        
        observeEvent(input$saveSession, {
            UIInterface$saveSession(list(removeFully = rValues$removeFully, removePartially = rValues$removePartially))
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
                keep <- tbl[keep == TRUE][[1]]
                notkeep <- tbl[keep == FALSE][[1]]
                oldsel <- rValues$removeFully
                rValues$removeFully <- setdiff(union(rValues$removeFully, notkeep), keep)
                selChanged <- !setequal(oldsel, rValues$removeFully)
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
                            updatePrimaryRow(tbl[[1]][selr])
                    }
                }
            }
        })
        
        observeEvent(input$primaryHot_select$select$r, {
            tbl <- rhandsontable::hot_to_r(input$primaryHot)
            updatePrimaryRow(tbl[[1]][input$primaryHot_select$select$rAll[1]])
        })
        
        observeEvent(input$enableAllPrimary, {
            rValues$removeFully <- character()
            rValues$triggerPrimaryHotUpdate <- rValues$triggerPrimaryHotUpdate + 1
            setSessionChanged(TRUE)
        })
        observeEvent(input$disableAllPrimary, {
            rValues$removeFully <- UIInterface$primarySelections
            rValues$triggerPrimaryHotUpdate <- rValues$triggerPrimaryHotUpdate + 1
            setSessionChanged(TRUE)
        })
        
        observeEvent(input$resetAllSecondary, {
            rValues$removePartially <- list()
            rValues$triggerSecondaryHotUpdate <- rValues$triggerSecondaryHotUpdate + 1
            setSessionChanged(TRUE)
        })
        
        observeEvent(input$secondaryHot, {
            # HACK: input$secondaryHot$params$maxRows: make sure we don't have empty table as hot_to_r errors otherwise
            if (input$secondaryHot$params$maxRows > 0)
            {
                tbl <- rhandsontable::hot_to_r(input$secondaryHot)
                oldsel <- rValues$removePartially[[rValues$currentPrimSel]]
                newsel <- UIInterface$getSecondarySelectionsFromTab(tbl)[!tbl$keep]
                if (length(newsel) == 0)
                    newsel <- NULL
                rValues$removePartially[[rValues$currentPrimSel]] <- newsel
                if (is.null(oldsel) != is.null(newsel) ||
                    (!is.null(oldsel) && !setequal(oldsel, newsel)))
                    setSessionChanged(TRUE)
            }
        })
        
        observeEvent(input$enableAllSecondary, {
            rValues$removePartially[[rValues$currentPrimSel]] <- NULL
            rValues$triggerSecondaryHotUpdate <- rValues$triggerSecondaryHotUpdate + 1
            setSessionChanged(TRUE)
        })
        observeEvent(input$disableAllSecondary, {
            rValues$removePartially[[rValues$currentPrimSel]] <-
                UIInterface$getSecondarySelections(rValues$currentPrimSel)
            rValues$triggerSecondaryHotUpdate <- rValues$triggerSecondaryHotUpdate + 1
            setSessionChanged(TRUE)
        })
        
        observeEvent(UIInterface$settingsChangedExpression(input), {
            shinyjs::toggleState("saveApplySettings", !isTRUE(all.equal(getInputSettings(), rValues$settings)))
        })
        
        observeEvent(input$saveApplySettings, {
            doApplySettings(getInputSettings())
            shinyjs::disable("saveApplySettings")
        })
        
        observeEvent(input$resetSettings, {
            doApplySettings(UIInterface$defaultUISettings())
            UIInterface$syncInputSettings(session, rValues$settings)
            shinyjs::disable("saveApplySettings")
        })

        UIInterface$doObserveEvents(input, rValues)
        
        output$pageTitle <- renderText({
            sprintf("%s (%d/%d)", rValues$currentPrimSel,
                    match(rValues$currentPrimSel, UIInterface$primarySelections), length(UIInterface$primarySelections))
        })
        
        output$plotMain <- renderPlot({ UIInterface$plotMain(input, rValues) })

        output$primaryHot <- rhandsontable::renderRHandsontable({
            not <- showNotification(sprintf("Updating %s table...", UIInterface$primaryTab()), duration = NULL,
                                    closeButton = FALSE, type = "message")
            rValues$triggerPrimaryHotUpdate
            
            tabData <- primaryTable()
            tabData[, keep := !.SD[[1]] %in% isolate(rValues$removeFully)]
            setcolorder(tabData, c(names(tabData)[1], "keep"))
            
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
                htmlwidgets::onRender(hotRenderJSFunc(UIInterface$primaryTab()))
            
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
            isolate(rp <- rValues$removePartially[[rValues$currentPrimSel]])
            isolate(tabData[, keep := is.null(rp) |
                                !UIInterface$getSecondarySelections(rValues$currentPrimSel) %in% rp])
            setcolorder(tabData, c(names(tabData)[1], "keep"))
            
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
                htmlwidgets::onRender(hotRenderJSFunc(UIInterface$secondaryTab()))
            
            # HACK: disable some default context options
            hot$x$contextMenu$items[c("undo", "redo", "alignment")] <- NULL
            
            # BUG: table is messed up after tab switch
            # work-around from https://github.com/jrowen/rhandsontable/issues/366
            outputOptions(output, "secondaryHot", suspendWhenHidden = FALSE)
            
            return(hot)
        })
    }
    
    runApp(shinyApp(getCheckUI(UIInterface, settings), server))
}
