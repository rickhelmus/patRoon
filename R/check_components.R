#' @include main.R
#' @include components.R
#' @include feature_groups.R
#' @include check_ui.R

# UNDONE: can't make components field type components?
checkComponentsInterface <- setRefClass("checkComponentsInterface", contains = "checkUIInterface",
                                        fields = c(components = "ANY", fGroups = "featureGroups", EICs = "list"))

checkComponentsInterface$methods(
    
    UITitle = function() "Check components tool",
    resetSecondaryUITitle = function() "Enable all feature groups for all components",
    primaryTabTopUI = function()
    {
        fillRow(
            width = 150,
            div(style = "margin-top: 8px; margin-left: 0px",
                checkboxGroupInput("plotSpec", NULL, c("Plot spectrum" = "plotSpec"), inline = TRUE))
        )        
    },
    settingsTabUI = function(settings)
    {
        fillRow(
            height = 200,
            fillCol(
                flex = NA,
                radioButtons("retUnit", "Retention time unit", c("Seconds" = "sec", "Minutes" = "min"),
                             settings$retUnit, inline = TRUE)
            )
        )
    },
    
    primaryTab = function() "Components",
    secondaryTab = function() "Feature groups",
    
    defaultUISettings = function()
    {
        return(list(retUnit = "sec"))
    },
    UISettingsFileName = function() "check_components.yml",
    
    settingsChangedExpression = function(input)
    {
        input$retUnit
    },
    primarySettingsChanged = function(cur, new) FALSE,
    secondarySettingsChanged = function(cur, new) FALSE,
    syncInputSettings = function(session, settings)
    {
        updateRadioButtons(session, "retUnit", selected = rValues$settings[["retUnit"]])
    },
    
    primaryTableData = function(rValues)
    {
        # make sure that right copy() is called... https://stackoverflow.com/a/16566247
        tab <- data.table::copy(componentInfo(components))
        if (rValues$settings$retUnit == "min")
        {
            for (col in c("cmp_ret", "cmp_retsd", "ret_increment", "ret_min", "ret_max", "ret_range"))
            {
                if (!is.null(tab[[col]]))
                    set(tab, j = col, value = tab[[col]] / 60)
            }
        }
        
        setcolorder(tab, "name")
        if (!is.null(tab[["links"]]))
            tab[, links := NULL]
        
        return(tab)
    },
    secondaryTableData = function(rValues)
    {
        tab <- data.table::copy(components[[rValues$currentPrimSel]])
        if (rValues$settings$retUnit == "min")
            tab[, ret := ret / 60]
        
        # NOTE: group names (ie secondary names) should be first column
        setcolorder(tab, "group")
        
        return(tab)
    },
    
    plotMain = function(input, rValues)
    {
        enFGroups <- rValues$secondarySelections$name[rValues$secondarySelections[[rValues$currentPrimSel]]]
        cmp <- components[, enFGroups]
        
        withr::with_par(list(mar = c(4, 4, 0.1, 1), cex = 1.5), {
            if ("plotSpec" %in% input$plotSpec)
            {
                scr <- split.screen(c(2, 1))
                screen(scr[1])
            }
            
            plotChroms(cmp, index = rValues$currentPrimSel, fGroups = fGroups, EICs = EICs,
                       title = "", retMin = rValues$settings$retUnit == "min")

            if ("plotSpec" %in% input$plotSpec)
            {
                screen(scr[2])
                plotSpectrum(cmp, index = rValues$currentPrimSel)
                
                close.screen(scr)
            }
        })
    }
)

#' @export
importCheckComponentsSession <- function(sessionIn, sessionOut, components, overWrite = FALSE)
{
    # UNDONE: docs
    
    aapply(checkmate::assertString, . ~ sessionIn + sessionOut, min.chars = 1)
    pathIn <- getCheckSessionPath(sessionIn, "components"); pathOut <- getCheckSessionPath(sessionOut, "components")
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFileExists(pathIn, "r", .var.name = "session", add = ac)
    checkmate::assertPathForOutput(pathOut, overwrite = TRUE, .var.name = "sessionOut", add = ac)
    checkmate::assertClass(components, "components", add = ac)
    checkmate::assertFlag(overWrite, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(components) == 0)
        stop("No components, nothing to do...")
    
    importCheckUISession(pathIn, pathOut, "components", "feature groups", names(components), groupNames(components),
                         overWrite = overWrite)
}

#' @rdname GUI-utils
#' @aliases checkComponents
#' @export
setMethod("checkComponents", "components", function(components, fGroups, session, rtWindow, clearSession)
{
    # UNDONE: update docs
    
    checkmate::assertClass(fGroups, "featureGroups") # do first so we can sync

    if (length(components) == 0 || length(fGroups) == 0)
        stop("No components or feature groups, nothing to check...")
    
    if (!all(groupNames(components) %in% names(fGroups)))
        stop("The components object contains results for feature groups not present in the given featureGroups object. ",
             "You may need to sync the components object, eg: components <- components[, names(fGroups)]")
    
    ac <- checkmate::makeAssertCollection()
    assertCheckComponentsSession(session, components, mustExist = FALSE, canClearSession = TRUE,
                                 didClearSession = clearSession, add = ac)
    checkmate::assertNumber(rtWindow, finite = TRUE, lower = 0, add = ac)
    checkmate::assertFlag(clearSession, add = ac)
    checkmate::reportAssertions(ac)
    
    sessionPath <- getCheckSessionPath(session, "components")
    checkmate::assertPathForOutput(sessionPath, overwrite = TRUE, .var.name = "session")

    if (clearSession && file.exists(sessionPath))
        file.remove(sessionPath)
    
    fGroups <- fGroups[, groupNames(components)] # remove any fGroups not in components
    
    cmpNames <- names(components)
    
    EICs <- getEICsForFGroups(fGroups, rtWindow = rtWindow, mzExpWindow = 0.001, topMost = 1,
                              topMostByRGroup = FALSE, onlyPresent = TRUE)
    
    curSession <- NULL
    if (file.exists(sessionPath))
        curSession <- readRDS(sessionPath)
    else
    {
        eg <- data.frame(name = names(fGroups))
        eg[, cmpNames] <- TRUE
        curSession <- list(primarySelections = cmpNames, secondarySelections = eg)
    }
    
    int <- checkComponentsInterface$new(components = components, fGroups = fGroups, EICs = EICs,
                                        primarySelections = cmpNames, curSession = curSession,
                                        sessionPath = sessionPath)
    return(runCheckUI(int))
})
