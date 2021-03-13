#' @include main.R
#' @include feature_groups.R
#' @include check_ui.R


checkFeaturesInterface <- setRefClass("checkFeaturesInterface", contains = "checkUIInterface",
                                      fields = c(fGroups = "featureGroups", EICsTopMost = "list",
                                                 EICsTopMostRG = "list", EICsAll = "list",
                                                 EICPreviews = "list"))

checkFeaturesInterface$methods(
    
    UITitle = function() "Check features tool",
    resetSecondaryUITitle = function() "Enable all features for all groups",
    primaryTabTopUI = function()
    {
        fillRow(
            width = 250,
            flex = c(NA, 1),
            fillRow(
                width = 75,
                div(style = "margin-top: 8px", HTML("<strong>Plot mode</strong>"))
            ),
            selectInput("fGroupPlotMode", NULL,
                        c("Top most group" = "topMost",
                          "Top most replicates" = "topMostByRGroup",
                          "All" = "all"), "topMost")
        )
    },
    settingsTabUI = function(settings)
    {
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
        )
    },
    
    primaryTab = function() "Feature groups",
    secondaryTab = function() "Features",
    
    defaultUISettings = function()
    {
        return(list(retUnit = "sec", featQuantity = "intensity", fGroupQuantity = "average",
                    fGroupColumns = c("retMZ", "estIDLevel", "totalScore"),
                    featureColumns = c("retMZ", "quantity", "totalScore")))
    },
    UISettingsFileName = function() "check_features.yml",
    
    getSecondarySelections = function(primSel)
    {
        fti <- groupFeatIndex(fGroups)[[primSel]]
        return(analyses(fGroups)[fti != 0])
    },
    
    init = function(rValues)
    {
        gNames <- names(fGroups)
        extraFGroups <- !all(curSession$removeFully %in% gNames) ||
            (length(curSession$removePartially) > 0 && !all(names(curSession$removePartially) %in% gNames))
        allSessionAnas <- unlist(curSession$removePartially)
        extraAnas <- length(allSessionAnas) > 0 && !all(allSessionAnas %in% analyses(fGroups))
        
        if (extraFGroups || extraAnas)
        {
            extraWhat <- if (extraFGroups && extraAnas) "feature groups and analyses" else if (extraFGroups) "feature groups" else "analyses"
            showModal(modalDialog(title = "Session data",
                                  easyClose = TRUE,
                                  paste(sprintf("Some additional selection data for %s is present in the loaded session.",
                                                extraWhat),
                                        "This can occur if the feature groups object was e.g. subset or filtered",
                                        "or the session was made for another feature groups object.",
                                        "In the latter case you probably want to use importCheckFeaturesSession() first.",
                                        "When you save the session now the additional selection data will be removed!")))
        }
        
        rValues$fGroupPlotMode <- "topMost"
        return(rValues)
    },
    
    settingsChangedExpression = function(input)
    {
        input$retUnit; input$featQuantity; input$fGroupQuantity; input$fGroupColumns; input$featureColumns
    },
    primarySettingsChanged = function(cur, new)
    {
        return(!isTRUE(all.equal(cur$fGroupColumns, new$fGroupColumns)) ||
                   cur$featQuantity != new$featQuantity ||
                   cur$fGroupQuantity != new$fGroupQuantity)
    },
    secondarySettingsChanged = function(cur, new) !isTRUE(all.equal(cur$featureColumns, new$featureColumns)),
    syncInputSettings = function(session, settings)
    {
        for (s in c("retUnit", "featQuantity", "fGroupQuantity"))
            updateRadioButtons(session, s, selected = settings[[s]])
        for (s in c("fGroupColumns", "featureColumns"))
            updateCheckboxGroupInput(session, s, selected = settings[[s]])
    },
    
    primaryTableData = function(rValues)
    {
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
            gData[, EIC := sapply(names(fGroups), function(g)
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
    },
    secondaryTableData = function(rValues)
    {
        fti <- groupFeatIndex(fGroups)[[rValues$currentPrimSel]]
        ft <- featureTable(fGroups)[fti != 0]; ai <- analysisInfo(fGroups)[fti != 0, ]; fti <- fti[fti != 0]
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
    },
    
    doObserveEvents = function(input, rValues)
    {
        observeEvent(input$tabs, {
            if (input$tabs != "Settings")
            {
                if (input$tabs == primaryTab())
                    rValues$fGroupPlotMode <- input$fGroupPlotMode
                else
                    rValues$fGroupPlotMode <- "all"
            }
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
    },
    
    plotMain = function(input, rValues)
    {
        EICs <- switch(rValues$fGroupPlotMode,
                       topMost = EICsTopMost,
                       topMostByRGroup = EICsTopMostRG,
                       all = EICsAll
        )
        if (length(EICs) == 0)
            EICs <- NULL # not (yet) loaded, in this case plotChroms() will make its own but EICs must be NULL
        
        fg <- fGroups[, rValues$currentPrimSel]
        if (rValues$fGroupPlotMode == "all") # UNDONE: also for rGroups top most somehow?
        {
            rp <- rValues$removePartially[[rValues$currentPrimSel]]
            if (!is.null(rp))
                fg <- fg[setdiff(getSecondarySelections(rValues$currentPrimSel), rp)]
        }
        
        withr::with_par(list(mar = c(4, 4, 0.1, 1), cex = 1.5), {
            plotChroms(fg, EICs = EICs, colourBy = "rGroups", showPeakArea = TRUE,
                       showFGroupRect = FALSE, title = "",
                       topMost = if (rValues$fGroupPlotMode == "all") NULL else 1,
                       topMostByRGroup = rValues$fGroupPlotMode == "topMostByRGroup",
                       retMin = rValues$settings$retUnit == "min")
        })
    },
    
    saveSession = function(s)
    {
        sessionGrps <- s$removeFully
        if (length(s$removePartially) > 0)
            sessionGrps <- union(sessionGrps, names(s$removePartially))
        saveCheckSession(s, session, fGroups[, sessionGrps], "featureGroups")
    }
)

#' @export
importCheckFeaturesSession <- function(sessionIn, sessionOut, fGroups, rtWindow = 6, mzWindow = 0.002,
                                       overWrite = FALSE)
{
    # UNDONE: docs
    
    ac <- checkmate::makeAssertCollection()
    assertCheckSession(sessionIn, mustExist = TRUE, add = ac)
    assertCheckSession(sessionOut, mustExist = FALSE, add = ac)
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    aapply(checkmate::assertNumber, . ~ rtWindow + mzWindow, lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::assertFlag(overWrite, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(fGroups) == 0)
    {
        printf("No feature groups, nothing to do...\n")
        invisible(return(NULL))
    }
    
    # UNDONE: handle overWrite
    
    oldSession <- readCheckSession(sessionIn, "featureGroups")
    
    if (length(oldSession$removeFully) == 0 && length(oldSession$removePartially) == 0)
    {
        printf("Old session is empty, nothing to do...\n")
        return(invisible(NULL))
    }
    
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
    
    newGroupsTab <- rbindlist(pruneList(newGroups), idcol = "oldGroup")

    removeFully <- newGroupsTab[oldGroup %chin% oldSession$removeFully]$group
    
    rmpwh <- which(newGroupsTab$oldGroup %chin% names(oldSession$removePartially))
    removePartially <- oldSession$removePartially[newGroupsTab$oldGroup[rmpwh]]
    names(removePartially) <- newGroupsTab$group[rmpwh]
    
    saveCheckSession(list(removeFully = removeFully, removePartially = removePartially), sessionOut,
                     fGroups[, newGroupsTab$group], "featureGroups")
    
    invisible(NULL)
    # importCheckUISession(pathIn, pathOut, "feature groups", "analyses", names(fGroups),
    #                      analyses(fGroups), overWrite = overWrite)
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
    
    ac <- checkmate::makeAssertCollection()
    assertCheckSession(session, mustExist = FALSE, add = ac)
    checkmate::assertNumber(rtWindow, finite = TRUE, lower = 0, add = ac)
    checkmate::assertFlag(clearSession, add = ac)
    checkmate::reportAssertions(ac)
    
    if (clearSession && file.exists(session))
        file.remove(session)
    
    gNames <- names(fGroups)
    fTable <- featureTable(fGroups)
    ftind <- groupFeatIndex(fGroups)
    
    EICsTopMost <- getEICsForFGroups(fGroups, rtWindow, 0.001, topMost = 1, topMostByRGroup = FALSE,
                                     onlyPresent = TRUE)
    EICsTopMostRG <- EICsAll <- list()
    
    # format is in [[ana]][[fGroup]], since we only took top most intensive we can throw away the ana dimension
    EICPreviews <- Reduce(modifyList, EICsTopMost)
    EICPreviews <- Map(EICPreviews, names(EICPreviews), f = function(eic, grp)
    {
        anai <- which.max(fGroups[[grp]])
        return(eic[numGTE(eic$time, fTable[[anai]]$retmin[ftind[[grp]][anai]]) &
                       numLTE(eic$time, fTable[[anai]]$retmax[ftind[[grp]][anai]]), ])
    })
    
    curSession <- NULL
    if (file.exists(session))
        curSession <- readCheckSession(session, "featureGroups")
    else
        curSession <- list(removeFully = character(), removePartially = list())
    
    int <- checkFeaturesInterface$new(fGroups = fGroups, EICsTopMost = EICsTopMost,
                                      EICsTopMostRG = EICsTopMostRG, EICsAll = EICsAll,
                                      EICPreviews = EICPreviews, primarySelections = gNames,
                                      curSession = curSession, session = session)
    return(runCheckUI(int))
})

convertQualitiesToMCData <- function(fGroups)
{
    if (!hasFGroupScores(fGroups))
        stop("No feature qualities were calculated. Please run calculatePeakQualities() first.")
    
    ret <- copy(groupQualities(fGroups))
    
    hasNA <- unique(unlist(lapply(ret, function(x) which(is.na(x)))))
    if (length(hasNA) > 0)
    {
        warning("The following feature groups have one or more NA peak qualities and will be omitted: ",
                paste0(ret$group[hasNA], collapse = ", "))
        ret <- ret[-hasNA]
    }
    
    ret[, EICNo := match(group, names(fGroups))]
    setcolorder(ret, "EICNo")
    qcols <- c(featureQualityNames(), featureGroupQualityNames())
    setnames(ret, qcols, paste0(qcols, "_mean"))
    ret[, group := NULL][]
    return(ret)
}

getMCTrainData <- function(fGroups, session)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups")
    assertCheckSession(session, mustExist = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    session <- readCheckSession(session, "featureGroups")
    ret <- convertQualitiesToMCData(fGroups)
    gNames <- names(fGroups)
    ret[, Class := fifelse(gNames[EICNo] %chin% session$removeFully, "BAD", "GOOD")]
    
    return(as.data.frame(ret))
}

predictCheckFeaturesSession <- function(fGroups, session, model = NULL, overWrite = FALSE)
{
    checkPackage("MetaClean")
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups")
    checkmate::assertClass(model, "train", null.ok = TRUE, add = ac)
    assertCheckSession(session, mustExist = FALSE, add = ac)
    checkmate::assertFlag(overWrite, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(fGroups) == 0)
        stop("No feature groups, nothing to do...")
    
    if (is.null(model))
    {
        if (!requireNamespace("MetaCleanData", quietly = TRUE))
            stop("Please install MetaCleanData to use the default example model")
        data("example_model", package = "MetaCleanData", envir = environment())
        model <- example_model
    }
    
    if (file.exists(session) && !overWrite)
        stop("Output session already exists. Set overWrite=TRUE to proceed anyway.")
    
    testd <- convertQualitiesToMCData(fGroups)
    preds <- MetaClean::getPredicitons(model = model, testData = testd, eicColumn = "EICNo")
    
    gNames <- names(fGroups)
    # UNDONE: when is it GOOD/BAD or Pass/Fail?
    rmf <- gNames[preds[preds$Pred_Class %in% c("GOOD", "Pass"), "EIC"]]
    saveCheckSession(list(removeFully = rmf, removePartially = list()), session, fGroups[, rmf], "featureGroups")
    
    invisible(NULL)
}
