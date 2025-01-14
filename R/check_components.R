# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

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
    
    getSecondarySelections = function(primSel) components[[primSel]]$group,
    getSecondarySelectionsFromTab = function(tab) tab[["group"]],

    init = function(rValues)
    {
        cNames <- names(components)
        extraComps <- !all(curSession$removeFully %in% cNames) ||
            (length(curSession$removePartially) > 0 && !all(names(curSession$removePartially) %in% cNames))
        allSessionFGroups <- unlist(curSession$removePartially)
        extraFGroups <- length(allSessionFGroups) > 0 && !all(allSessionFGroups %in% groupNames(components))
        
        if (extraComps || extraFGroups)
        {
            extraWhat <- if (extraComps && extraFGroups) "components and feature groups" else if (extraComps) "components" else "feature groups"
            showModal(modalDialog(title = "Session data",
                                  easyClose = TRUE,
                                  paste(sprintf("Some additional selection data for %s is present in the loaded session.",
                                                extraWhat),
                                        "This can occur if the components object was e.g. subset or filtered",
                                        "or the session was made for another components object.",
                                        "In the latter case you probably want to use importCheckComponentsSession() first.",
                                        "When you save the session now the additional selection data will be removed!")))
        }
        
        return(rValues)
    },
    
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
        rp <- rValues$removePartially[[rValues$currentPrimSel]]
        cmp <- if (!is.null(rp)) delete(components, j = rp) else components
        
        bg <- if (rValues$currentPrimSel %in% rValues$removeFully) RColorBrewer::brewer.pal(9, "Reds")[[1]] else "white"
        withr::with_par(list(mar = c(4, 4, 0.1, 1), cex = 1.5, bg = bg), {
            if (!rValues$currentPrimSel %in% names(cmp))
            {
                # may happen if all fGroups are disabled
                noDataPlot()
            }
            else
            {
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
            }
        })
    },
    
    saveSession = function(s)
    {
        # remove oldd selections if present, eg now removed due to subsetting the components
        s$removeFully <- intersect(s$removeFully, names(components))
        s$removePartially <- s$removePartially[names(s$removePartially) %in% names(components)]
        sessionGrps <- character()
        if (length(s$removePartially) > 0)
            sessionGrps <- unlist(s$removePartially)
        saveCheckSession(s, session, fGroups[, sessionGrps], "components")
    }
)

#' @details \code{checkComponents} is used to review components and their feature groups contained within. A typical use
#'   case is to verify that peaks from features that were annotated as related adducts and/or isotopes are correctly
#'   aligned.
#'
#' @param components The \code{\link{components}} to be checked.
#'
#' @note \code{checkComponents}: Some componentization algorithms (\emph{e.g.} \code{\link{generateComponentsNontarget}}
#'   and \code{\link{generateComponentsTPs}}) may output components where the same feature group in a component is
#'   present multiple times, for instance, when multiple TPs are matched to the same feature group. If such a feature
#'   group is selected for removal, then \emph{all} of its result in the component will be marked for removal.
#'
#' @rdname check-GUI
#' @aliases checkComponents
#' @export
setMethod("checkComponents", "components", function(components, fGroups, session, EICParams, clearSession)
{
    checkmate::assertClass(fGroups, "featureGroups") # do first so we can sync
    checkmate::assertFlag(clearSession)
    
    if (length(components) == 0 || length(fGroups) == 0)
        stop("No components or feature groups, nothing to check...")
    
    if (!all(groupNames(components) %in% names(fGroups)))
        stop("The components object contains results for feature groups not present in the given featureGroups object. ",
             "You may need to sync the components object, eg: components <- components[, names(fGroups)]")
    
    ac <- checkmate::makeAssertCollection()
    assertCheckSession(session, mustExist = FALSE, add = ac)
    assertEICParams(EICParams, add = ac)
    checkmate::reportAssertions(ac)
    
    if (clearSession && file.exists(session))
        file.remove(session)
    
    fGroups <- fGroups[, groupNames(components)] # remove any fGroups not in components
    
    cmpNames <- names(components)
    
    EICs <- getEICsForFGroups(fGroups, EICParams = modifyList(EICParams, list(topMost = 1, topMostByRGroup = FALSE)))
    
    curSession <- NULL
    if (file.exists(session))
        curSession <- readCheckSession(session, "components")
    else
        curSession <- list(removeFully = character(), removePartially = list())
    
    int <- checkComponentsInterface$new(components = collapseComponents(components), fGroups = fGroups, EICs = EICs,
                                        primarySelections = cmpNames, curSession = curSession,
                                        session = session)
    return(runCheckUI(int))
})
