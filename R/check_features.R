#' @include main.R
#' @include feature_groups.R
#' @include check_ui.R


checkFeaturesInterface <- setRefClass("checkFeaturesInterface", contains = "checkUIInterface",
                                      fields = c(fGroups = "featureGroups", EICParams = "list",
                                                 EICsTopMost = "list", EICsTopMostRG = "list", EICsAll = "list",
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
                                     "Ion annotations" = "ionAnn",
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
    getSecondarySelectionsFromTab = function(tab) tab[["analysis"]],
    
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
        gData <- do.call(as.data.table, args)
        
        if (rValues$settings$fGroupQuantity == "max")
        {
            anaInfo <- analysisInfo(fGroups)
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
        
        if (!"ionAnn" %in% rValues$settings$fGroupColumns)
        {
            if (!is.null(gData[["adduct"]]))
                gData[, adduct := NULL]
            if (!is.null(gData[["neutralMass"]]))
                gData[, neutralMass := NULL]
        }
        
        if (getScores && hasFGroupScores(fGroups))
        {
            # scores were added by as.data.table(). Remove those we don't want.
            if (!"totalScore" %in% rValues$settings$fGroupColumns)
                gData[, totalScore := NULL]
            if (!"otherScores" %in% rValues$settings$fGroupColumns)
                gData[, (featureQualityNames(scores = TRUE, totScore = FALSE)) := NULL]
        }
        
        return(gData)
    },
    secondaryTableData = function(rValues)
    {
        fti <- groupFeatIndex(fGroups)[[rValues$currentPrimSel]]
        ft <- featureTable(fGroups)[fti != 0]; ai <- analysisInfo(fGroups)[fti != 0]; fti <- fti[fti != 0]
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
            {
                fs <- featureQualityNames(group = FALSE, scores = TRUE, totScore = FALSE)
                fData[, (fs) := feat[, fs, with = FALSE]]
            }
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
            if ((input$fGroupPlotMode == "topMostByRGroup" && length(EICsTopMostRG) == 0) ||
                (input$fGroupPlotMode == "all" && length(EICsAll) == 0))
            {
                not <- showNotification("Loading EICs...", duration = NULL, closeButton = FALSE, type = "message")
                if (input$fGroupPlotMode == "topMostByRGroup")
                    EICsTopMostRG <<- getEICsForFGroups(fGroups,
                                                        EICParams = modifyList(EICParams, list(topMost = 1,
                                                                                               topMostByRGroup = TRUE)))
                else
                    EICsAll <<- getEICsForFGroups(fGroups, EICParams = modifyList(EICParams, list(topMost = NULL),
                                                                                  keep.null = TRUE))
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
        
        ep <- getDefEICParams(topMost = if (rValues$fGroupPlotMode == "all") NULL else 1,
                              topMostByRGroup = rValues$fGroupPlotMode == "topMostByRGroup")
        withr::with_par(list(mar = c(4, 4, 0.1, 1), cex = 1.5), {
            plotChroms(fg, EICs = EICs, colourBy = "rGroups", showPeakArea = TRUE, EICParams = ep,
                       showFGroupRect = FALSE, title = "", retMin = rValues$settings$retUnit == "min")
        })
    },
    
    saveSession = function(s)
    {
        # remove oldd selections if present, eg now removed due to subsetting fGroups
        s$removeFully <- intersect(s$removeFully, names(fGroups))
        s$removePartially <- s$removePartially[names(s$removePartially) %in% names(fGroups)]
        sessionGrps <- s$removeFully
        if (length(s$removePartially) > 0)
            sessionGrps <- union(sessionGrps, names(s$removePartially))
        saveCheckSession(s, session, fGroups[, sessionGrps], "featureGroups")
    }
)

#' @details \code{importCheckFeaturesSession} is used to import a session file that was generated from a different
#'   \code{\link{featureGroups}} object. This is useful to avoid re-doing manual interpretation of chromatographic peaks
#'   when, for instance, feature group data is re-created with different parameters.
#'
#' @param sessionIn,sessionOut The file names for the input and output sessions.
#' @param rtWindow The retention time window (seconds) used to relate 'old' with 'new' feature groups.
#' @param mzWindow The \emph{m/z} window (in Da) used to relate 'old' with 'new' feature groups.
#' @param overWrite Set to \code{TRUE} to overwrite the output session file if it already exists. If \code{FALSE}, the
#'   function will stop with an error message.
#' 
#' @rdname check-GUI
#' @export
importCheckFeaturesSession <- function(sessionIn, sessionOut, fGroups, rtWindow = 6, mzWindow = 0.002,
                                       overWrite = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    assertCheckSession(sessionIn, mustExist = TRUE, add = ac)
    assertCheckSession(sessionOut, mustExist = FALSE, add = ac)
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    aapply(checkmate::assertNumber, . ~ rtWindow + mzWindow, lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::assertFlag(overWrite, add = ac)
    checkmate::reportAssertions(ac)
    
    if (file.exists(sessionOut) && !overWrite)
        stop("Output session already exists. Set overWrite=TRUE to proceed anyway.")
    
    if (length(fGroups) == 0)
    {
        printf("No feature groups, nothing to do...\n")
        return(invisible(NULL))
    }
    
    oldSession <- readCheckSession(sessionIn, "featureGroups")
    
    if (length(oldSession$removeFully) == 0 && length(oldSession$removePartially) == 0)
    {
        printf("Old session is empty, nothing to do...\n")
        return(invisible(NULL))
    }
    
    newGroupsTab <- importCheckUISessionGroups(oldSession, fGroups, rtWindow, mzWindow)
    
    if (nrow(newGroupsTab) == 0)
    {
        printf("Nothing could be matched, nothing to do...\n")
        return(invisible(NULL))
    }
    
    removeFully <- newGroupsTab[oldGroup %chin% oldSession$removeFully]$group
    
    rmpwh <- which(newGroupsTab$oldGroup %chin% names(oldSession$removePartially))
    removePartially <- oldSession$removePartially[newGroupsTab$oldGroup[rmpwh]]
    names(removePartially) <- newGroupsTab$group[rmpwh]
    
    saveCheckSession(list(removeFully = removeFully, removePartially = removePartially), sessionOut,
                     fGroups[, newGroupsTab$group], "featureGroups")
    
    invisible(NULL)
}

#' @details \code{checkFeatures} is used to review chromatographic information for feature groups. Its main purpose is
#'   to assist in reviewing the quality of detected feature (groups) and easily select unwanted data such as features
#'   with poor peak shapes or noise.
#' @rdname check-GUI
#' @aliases checkFeatures
#' @export
setMethod("checkFeatures", "featureGroups", function(fGroups, session, EICParams, clearSession)
{
    if (length(fGroups) == 0)
        stop("No feature groups, nothing to check...")
    
    ac <- checkmate::makeAssertCollection()
    assertCheckSession(session, mustExist = FALSE, add = ac)
    assertEICParams(EICParams, add = ac)
    checkmate::assertFlag(clearSession, add = ac)
    checkmate::reportAssertions(ac)
    
    if (clearSession && file.exists(session))
        file.remove(session)
    
    gNames <- names(fGroups)
    fTable <- featureTable(fGroups)
    ftind <- groupFeatIndex(fGroups)
    
    EICsTopMost <- getEICsForFGroups(fGroups, EICParams = modifyList(EICParams, list(topMost = 1,
                                                                                     topMostByRGroup = FALSE)))
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
    
    int <- checkFeaturesInterface$new(fGroups = fGroups, EICParams = EICParams, EICsTopMost = EICsTopMost,
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
    qcols <- featureQualityNames()
    setnames(ret, qcols, paste0(qcols, "_mean"))
    ret[, group := NULL][]
    return(ret)
}

#' @details \code{getMCTrainData} converts a session created by \code{checkFeatures} to a \code{data.frame} that can be
#'   used by the \pkg{MetaClean} to train a new model. The output format is comparable to that from
#'   \code{\link[MetaClean]{getPeakQualityMetrics}}.
#' @note \code{getMCTrainData} only uses session data for selected feature groups. Selected features for removal are
#'   ignored, as this is not supported by \pkg{MetaClean}.
#' @rdname check-GUI
#' @export
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

#' @details \code{predictCheckFeaturesSession} Uses ML data from \pkg{MetaClean} to predict the quality (Pass/Fail) of
#'   feature group data, and converts this to a session which can be reviewed with \code{checkFeatures} and used to
#'   remove unwanted feature groups by \code{\link[=filter,featureGroups-method]{filter}}.
#' @param model The model that was created with \pkg{MetaClean} and that should be used to predict pass/fail data. If
#'   \code{NULL}, the example model of the \pkg{MetaCleanData} package is used.
#' @rdname check-GUI
#' @export
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
    rmf <- gNames[preds[preds$Pred_Class %in% c("BAD", "Fail"), "EIC"]]
    saveCheckSession(list(removeFully = rmf, removePartially = list()), session, fGroups[, rmf], "featureGroups")
    
    invisible(NULL)
}
