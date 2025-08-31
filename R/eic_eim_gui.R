#' Launch EIC/EIM GUI
#'
#' Launches a Shiny app for generating and plotting EICs and EIMs.
#'
#' @param obj The object to use for analysis information or features/featureGroups.
#' @param suspects For data.frame method, an optional suspect list data.frame.
#' @return A Shiny app object.
#' @export
setMethod("launchEICGUI", "data.frame", function(obj, suspects = NULL, adduct = NULL)
{
    if (!is.null(adduct))
        adduct <- checkAndToAdduct(adduct)
    
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(obj, add = ac)
    assertSuspectList(suspects, needsAdduct = is.null(adduct), skipInvalid = TRUE, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!is.null(suspects))
    {
        suspects <- prepareSuspectList(suspects, adduct = adduct, skipInvalid = TRUE, checkDesc = TRUE,
                                       prefCalcChemProps = TRUE, neutralChemProps = FALSE)
        
        if (!is.null(suspects[["mobility"]]) || !is.null(suspects[["CCS"]]))
        {
            suspects[, mobility_susp := selectFromSuspAdductCol(suspects, "mobility", data.table(), if (!is.null(adduct)) as.character(adduct))]
            suspects <- expandSuspMobilities(suspects)
        }
    }
    gui <- createEICGUI(obj, analysisInfo, suspects)
    shiny::shinyApp(gui$ui, gui$server)
})

setMethod("launchEICGUI", "features", function(obj)
{
    analysisInfo <- analysisInfo(obj)
    gui <- createEICGUI(obj, analysisInfo)
    shiny::shinyApp(gui$ui, gui$server)
})

setMethod("launchEICGUI", "featureGroups", function(obj)
{
    analysisInfo <- analysisInfo(obj)
    gui <- createEICGUI(obj, analysisInfo)
    shiny::shinyApp(gui$ui, gui$server)
})

setMethod("launchEICGUI", "featureGroupsScreening", function(obj)
{
    analysisInfo <- analysisInfo(obj)
    suspects <- screenInfo(obj)
    gui <- createEICGUI(obj, analysisInfo, suspects)
    shiny::shinyApp(gui$ui, gui$server)
})

createEICGUI <- function(obj, analysisInfo, suspects = NULL)
{
    featureChoices <- NULL
    suspectChoices <- NULL
    if (is(obj, "features"))
    {
        ft <- featureTable(obj)[[analysisInfo$analysis[1]]]
        featureChoices <- seq_len(nrow(ft))
        names(featureChoices) <- sapply(seq_len(nrow(ft)), function(i) {
            f <- ft[i, ]
            label <- sprintf("%s (mz: %.5f, rt: %.2f", f$ID, f$mz, f$ret)
            if ("mobility" %in% names(ft) && !is.na(f$mobility)) label <- paste(label, sprintf(", mob: %.3f)", f$mobility))
            else label <- paste(label, ")")
            label
        })
    }
    else if (is(obj, "featureGroups"))
    {
        ft <- featureTable(obj)[[analysisInfo$analysis[1]]]
        groups <- unique(ft$group)
        featureChoices <- seq_along(groups)
        names(featureChoices) <- sapply(seq_along(groups), function(i) {
            g <- groups[i]
            gf <- ft[group == g][1]  # Use first feature in group
            label <- sprintf("%s (mz: %.5f, rt: %.2f", g, gf$mz, gf$ret)
            if ("mobility" %in% names(gf) && !is.na(gf$mobility)) label <- paste(label, sprintf(", mob: %.3f)", gf$mobility))
            else label <- paste(label, ")")
            label
        })
    }
    if (!is.null(suspects))
    {
        suspectChoices <- seq_len(nrow(suspects))
        names(suspectChoices) <- sapply(seq_len(nrow(suspects)), function(i) {
            s <- suspects[i, ]
            label <- s$name
            hasMZ <- "mz" %in% names(suspects) && !is.na(s$mz)
            hasRT <- "rt" %in% names(suspects) && !is.na(s$rt)
            hasMob <- "mobility" %in% names(suspects) && !is.na(s$mobility)
            hasCCS <- "CCS" %in% names(suspects) && !is.na(s$CCS)
            hasGroup <- is(obj, "featureGroupsScreening") && "group" %in% names(suspects) && !is.na(s$group)
            if (hasMZ || hasRT || hasMob || hasCCS || hasGroup)
            {
                details <- c()
                if (hasMZ) details <- c(details, sprintf("mz: %.5f", s$mz))
                if (hasRT) details <- c(details, sprintf("rt: %.2f", s$rt))
                if (hasMob) details <- c(details, sprintf("mob: %.3f", s$mobility))
                if (hasCCS) details <- c(details, sprintf("CCS: %.1f", s$CCS))
                if (hasGroup) details <- c(details, sprintf("group: %s", s$group))
                label <- paste(label, sprintf("(%s)", paste(details, collapse = ", ")))
            }
            label
        })
    }
    
    ui <- shiny::fluidPage(
        shiny::titlePanel("EIC/EIM Viewer"),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::selectInput("analysis", "Analysis", choices = analysisInfo$analysis, selected = analysisInfo$analysis[1]),
                if (!is.null(featureChoices))
                    shiny::selectInput("feature", "Select Feature", choices = featureChoices),
                if (!is.null(suspectChoices))
                    shiny::selectInput("suspect", "Select Suspect", choices = suspectChoices),
                shiny::radioButtons("mzMode", "MZ Specification", choices = c("Range", "Center +/- Window", "Neutral Mass +/- MZ Window", "Formula +/- MZ Window"), selected = if (is(obj, "features") || is(obj, "featureGroups")) "Center +/- Window" else "Range", inline = TRUE),
                shiny::conditionalPanel(
                    condition = "input.mzMode == 'Range'",
                    shiny::fluidRow(
                        shiny::column(6, shiny::numericInput("mzmin", "MZ Min", value = 100, min = 0)),
                        shiny::column(6, shiny::numericInput("mzmax", "MZ Max", value = 200, min = 0))
                    )
                ),
                shiny::conditionalPanel(
                    condition = "input.mzMode == 'Center +/- Window'",
                    shiny::fluidRow(
                        shiny::column(6, shiny::numericInput("mzcenter", "MZ Center", value = 150, min = 0)),
                        shiny::column(6, shiny::numericInput("mzwindow", "MZ Window", value = 0.005, min = 0))
                    )
                ),
                shiny::conditionalPanel(
                    condition = "input.mzMode == 'Neutral Mass +/- MZ Window'",
                    shiny::fluidRow(
                        shiny::column(4, shiny::numericInput("neutralMass", "Neutral Mass", value = 100, min = 0)),
                        shiny::column(4, shiny::selectInput("adductNM", "Adduct", choices = c("[M+H]+", "[M-H]-", "[M+Na]+", "[M+K]+", "[2M+H]+", "[2M-H]-", "[2M+Na]+", "[2M+K]+"))),
                        shiny::column(4, shiny::numericInput("mzWindowNM", "MZ Window", value = 0.005, min = 0))
                    )
                ),
                shiny::conditionalPanel(
                    condition = "input.mzMode == 'Formula +/- MZ Window'",
                    shiny::fluidRow(
                        shiny::column(4, shiny::textInput("formula", "Formula", value = "C6H12O6")),
                        shiny::column(4, shiny::selectInput("adductF", "Adduct", choices = c("[M+H]+", "[M-H]-", "[M+Na]+", "[M+K]+", "[2M+H]+", "[2M-H]-", "[2M+Na]+", "[2M+K]+"))),
                        shiny::column(4, shiny::numericInput("mzWindowF", "MZ Window", value = 0.005, min = 0))
                    )
                ),
                shiny::radioButtons("rtMode", "RT Specification", choices = c("Disabled", "Range", "Value +/- Window"), selected = if (is(obj, "features") || is(obj, "featureGroups")) "Value +/- Window" else "Disabled", inline = TRUE),
                shiny::conditionalPanel(
                    condition = "input.rtMode == 'Range'",
                    shiny::fluidRow(
                        shiny::column(6, shiny::numericInput("retmin", "RT Min", value = 0, min = 0)),
                        shiny::column(6, shiny::numericInput("retmax", "RT Max", value = 10, min = 0))
                    )
                ),
                shiny::conditionalPanel(
                    condition = "input.rtMode == 'Value +/- Window'",
                    shiny::fluidRow(
                        shiny::column(6, shiny::numericInput("retvalue", "RT Value", value = 5, min = 0)),
                        shiny::column(6, shiny::numericInput("retwindow", "RT Window", value = 0.5, min = 0))
                    )
                ),
                shiny::radioButtons("mobMode", "Mobility Specification", choices = c("Disabled", "Range", "Value +/- Window"), selected = if (is(obj, "features") || is(obj, "featureGroups")) "Value +/- Window" else "Disabled", inline = TRUE),
                shiny::conditionalPanel(
                    condition = "input.mobMode == 'Range'",
                    shiny::fluidRow(
                        shiny::column(6, shiny::numericInput("mobmin", "Mobility Min", value = 0, min = 0, step = 0.01)),
                        shiny::column(6, shiny::numericInput("mobmax", "Mobility Max", value = 1.3, min = 0, step = 0.01))
                    )
                ),
                shiny::conditionalPanel(
                    condition = "input.mobMode == 'Value +/- Window'",
                    shiny::fluidRow(
                        shiny::column(6, shiny::numericInput("mobvalue", "Mobility Value", value = 0.5, min = 0, step = 0.01)),
                        shiny::column(6, shiny::numericInput("mobwindow", "Mobility Window", value = 0.1, min = 0, step = 0.01))
                    )
                ),
                shiny::checkboxInput("plotEIM", "Also plot EIM", value = FALSE),
                shiny::conditionalPanel(
                    condition = "input.plotEIM == true",
                    shiny::checkboxInput("enableEIMSmoothing", "Enable EIM Smoothing", value = FALSE),
                    shiny::conditionalPanel(
                        condition = "input.enableEIMSmoothing == true",
                        shiny::numericInput("sgLength", "Smoothing Length", value = 5, min = 3, step = 2)
                    )
                ),
                shiny::checkboxInput("autoGenerate", "Auto-generate plots on input change", value = TRUE),
                shiny::conditionalPanel(
                    condition = "input.autoGenerate == false",
                    shiny::actionButton("generate", "Generate plot(s)")
                )
            ),
            shiny::mainPanel(
                shiny::plotOutput("eicPlot"),
                shiny::conditionalPanel(
                    condition = "input.plotEIM == true",
                    shiny::plotOutput("eimPlot")
                )
            )
        )
    )
    
    server <- function(input, output, session)
    {
        if (!is.null(featureChoices))
        {
            observeEvent(input$analysis, {
                if (is(obj, "features"))
                {
                    ft <- featureTable(obj)[[input$analysis]]
                    featureChoices <- seq_len(nrow(ft))
                    names(featureChoices) <- sapply(seq_len(nrow(ft)), function(i) {
                        f <- ft[i, ]
                        label <- sprintf("%s (mz: %.5f, rt: %.2f", f$ID, f$mz, f$ret)
                        if ("mobility" %in% names(ft) && !is.na(f$mobility)) label <- paste(label, sprintf(", mob: %.3f)", f$mobility))
                        else label <- paste(label, ")")
                        label
                    })
                    updateSelectInput(session, "feature", choices = featureChoices)
                }
                else if (is(obj, "featureGroups"))
                {
                    ft <- featureTable(obj)[[input$analysis]]
                    groups <- unique(ft$group)
                    featureChoices <- seq_along(groups)
                    names(featureChoices) <- sapply(seq_along(groups), function(i) {
                        g <- groups[i]
                        gf <- ft[group == g][1]  # Use first feature in group
                        label <- sprintf("%s (mz: %.5f, rt: %.2f", g, gf$mz, gf$ret)
                        if ("mobility" %in% names(gf) && !is.na(gf$mobility)) label <- paste(label, sprintf(", mob: %.3f)", gf$mobility))
                        else label <- paste(label, ")")
                        label
                    })
                    updateSelectInput(session, "feature", choices = featureChoices)
                }
                # For featureGroups, no change needed
            })
            observeEvent(input$feature, {
                if (is(obj, "features"))
                {
                    ft <- featureTable(obj)[[input$analysis]]
                    feat <- ft[as.integer(input$feature)]
                }
                else if (is(obj, "featureGroups"))
                {
                    ft <- featureTable(obj)[[input$analysis]]
                    feat <- ft[group == names(obj)[as.integer(input$feature)]][1]  # Use first feature in the group
                }
                if (nrow(feat) == 0)
                    return()
                # Update MZ
                if (length(feat$mz) > 0 && !is.na(feat$mz))
                {
                    updateNumericInput(session, "mzcenter", value = round(feat$mz, 5))
                    updateNumericInput(session, "mzwindow", value = round((feat$mzmax - feat$mzmin) / 2, 5))
                    updateNumericInput(session, "mzmin", value = round(feat$mzmin, 5))
                    updateNumericInput(session, "mzmax", value = round(feat$mzmax, 5))
                }
                # Update RT
                if (length(feat$ret) > 0 && !is.na(feat$ret))
                {
                    updateNumericInput(session, "retvalue", value = round(feat$ret, 2))
                    updateNumericInput(session, "retwindow", value = round((feat$retmax - feat$retmin) / 2, 2))
                    updateNumericInput(session, "retmin", value = round(feat$retmin, 2))
                    updateNumericInput(session, "retmax", value = round(feat$retmax, 2))
                }
                # Update Mobility
                print(feat)
                if (length(feat$mobility) > 0 && !is.na(feat$mobility) && length(feat$mobmin) > 0 && !is.na(feat$mobmin))
                {
                    updateNumericInput(session, "mobvalue", value = round(feat$mobility, 3))
                    updateNumericInput(session, "mobwindow", value = round((feat$mobmax - feat$mobmin) / 2, 3))
                    updateNumericInput(session, "mobmin", value = round(feat$mobmin, 3))
                    updateNumericInput(session, "mobmax", value = round(feat$mobmax, 3))
                }
            })
        }
        if (!is.null(suspects))
        {
            observeEvent(input$suspect, {
                susp <- suspects[as.integer(input$suspect)]
                if (nrow(susp) == 0)
                    return()
                # For featureGroupsScreening, use the corresponding feature
                if (is(obj, "featureGroupsScreening") && "group" %in% names(suspects) && !is.na(susp$group))
                {
                    ft <- featureTable(obj)[[input$analysis]]
                    feat <- ft[group == susp$group][1]
                    if (nrow(feat) > 0)
                    {
                        # Update from feature
                        if (length(feat$mz) > 0 && !is.na(feat$mz))
                        {
                            updateNumericInput(session, "mzcenter", value = round(feat$mz, 5))
                            updateNumericInput(session, "mzwindow", value = round((feat$mzmax - feat$mzmin) / 2, 5))
                            updateNumericInput(session, "mzmin", value = round(feat$mzmin, 5))
                            updateNumericInput(session, "mzmax", value = round(feat$mzmax, 5))
                        }
                        if (length(feat$ret) > 0 && !is.na(feat$ret))
                        {
                            updateNumericInput(session, "retvalue", value = round(feat$ret, 2))
                            updateNumericInput(session, "retwindow", value = round((feat$retmax - feat$retmin) / 2, 2))
                            updateNumericInput(session, "retmin", value = round(feat$retmin, 2))
                            updateNumericInput(session, "retmax", value = round(feat$retmax, 2))
                        }
                        if (length(feat$mobility) > 0 && !is.na(feat$mobility) && length(feat$mobmin) > 0 && !is.na(feat$mobmin))
                        {
                            updateNumericInput(session, "mobvalue", value = round(feat$mobility, 3))
                            updateNumericInput(session, "mobwindow", value = round((feat$mobmax - feat$mobmin) / 2, 3))
                            updateNumericInput(session, "mobmin", value = round(feat$mobmin, 3))
                            updateNumericInput(session, "mobmax", value = round(feat$mobmax, 3))
                        }
                    }
                }
                else
                {
                    # Update from suspect
                    if (length(susp$mz) > 0 && !is.na(susp$mz))
                    {
                        updateNumericInput(session, "mzcenter", value = round(susp$mz, 5))
                        updateNumericInput(session, "mzwindow", value = 0.005)
                        updateNumericInput(session, "mzmin", value = round(susp$mz - 0.005, 5))
                        updateNumericInput(session, "mzmax", value = round(susp$mz + 0.005, 5))
                    }
                    if (length(susp$rt) > 0 && !is.na(susp$rt))
                    {
                        updateNumericInput(session, "retvalue", value = round(susp$rt, 2))
                        updateNumericInput(session, "retwindow", value = 0.5)
                        updateNumericInput(session, "retmin", value = round(susp$rt - 0.5, 2))
                        updateNumericInput(session, "retmax", value = round(susp$rt + 0.5, 2))
                    }
                    if (length(susp$mobility) > 0 && !is.na(susp$mobility))
                    {
                        updateNumericInput(session, "mobvalue", value = round(susp$mobility, 3))
                        updateNumericInput(session, "mobwindow", value = 0.1)
                        updateNumericInput(session, "mobmin", value = round(susp$mobility - 0.1, 3))
                        updateNumericInput(session, "mobmax", value = round(susp$mobility + 0.1, 3))
                    }
                }
                # Update Adduct and Formula from suspect
                if (length(susp$adduct) > 0 && !is.na(susp$adduct))
                {
                    updateSelectInput(session, "adductNM", selected = susp$adduct)
                    updateSelectInput(session, "adductF", selected = susp$adduct)
                }
                if (length(susp$formula) > 0 && !is.na(susp$formula))
                {
                    updateTextInput(session, "formula", value = susp$formula)
                }
            })
        }
        
        eicData <- shiny::reactive({
            req(input$mzMode, input$rtMode, input$mobMode, input$analysis)
            if (!input$autoGenerate)
                req(input$generate)
            
            mzmin <- mzmax <- NULL
            if (input$mzMode == "Range")
            {
                mzmin <- input$mzmin
                mzmax <- input$mzmax
            }
            else if (input$mzMode == "Center +/- Window")
            {
                mzmin <- input$mzcenter - input$mzwindow
                mzmax <- input$mzcenter + input$mzwindow
            }
            else if (input$mzMode == "Neutral Mass +/- MZ Window")
            {
                neutralMass <- input$neutralMass
                mz <- tryCatch(calculateMasses(neutralMass, as.adduct(input$adductNM), "mz"), error = function(e)
                {
                    shiny::showNotification("Error calculating MZ from neutral mass and adduct", type = "error")
                    NA
                })
                if (!is.na(mz))
                {
                    mzmin <- mz - input$mzWindowNM
                    mzmax <- mz + input$mzWindowNM
                }
            }
            else if (input$mzMode == "Formula +/- MZ Window")
            {
                neutralMass <- tryCatch(getFormulaMass(input$formula), error = function(e)
                {
                    shiny::showNotification("Error parsing formula", type = "error")
                    NA
                })
                if (!is.na(neutralMass))
                {
                    mz <- tryCatch(calculateMasses(neutralMass, as.adduct(input$adductF), "mz"), error = function(e)
                    {
                        shiny::showNotification("Error calculating MZ from formula and adduct", type = "error")
                        NA
                    })
                    if (!is.na(mz))
                    {
                        mzmin <- mz - input$mzWindowF
                        mzmax <- mz + input$mzWindowF
                    }
                }
            }

            if (is.null(mzmin) || is.null(mzmax))
            {
                mzmin <- 100
                mzmax <- 200
            }
            ranges <- data.table::data.table(
                mzmin = mzmin,
                mzmax = mzmax,
                retmin = if (input$rtMode == "Disabled") 0 else if (input$rtMode == "Range") input$retmin else input$retvalue - input$retwindow,
                retmax = if (input$rtMode == "Disabled") 0 else if (input$rtMode == "Range") input$retmax else input$retvalue + input$retwindow
            )
            if (input$mobMode != "Disabled")
            {
                ranges$mobmin <- if (input$mobMode == "Range") input$mobmin else input$mobvalue - input$mobwindow
                ranges$mobmax <- if (input$mobMode == "Range") input$mobmax else input$mobvalue + input$mobwindow
            }
            rngList <- list(ranges)
            names(rngList) <- input$analysis
            getEICs(analysisInfo[analysis == input$analysis], rngList, gapFactor = 3, output = "pad")
        })
        
        eimData <- shiny::reactive({
            req(input$plotEIM, input$mzMode, input$rtMode, input$mobMode, input$analysis)
            if (!input$autoGenerate)
                req(input$generate)
            
            if (!input$plotEIM)
            {
                return(NULL)
            }
            mzmin <- mzmax <- NULL
            if (input$mzMode == "Range")
            {
                mzmin <- input$mzmin
                mzmax <- input$mzmax
            }
            else if (input$mzMode == "Center +/- Window")
            {
                mzmin <- input$mzcenter - input$mzwindow
                mzmax <- input$mzcenter + input$mzwindow
            }
            else if (input$mzMode == "Neutral Mass +/- MZ Window")
            {
                neutralMass <- input$neutralMass
                mz <- tryCatch(calculateMasses(neutralMass, as.adduct(input$adductNM), "mz"), error = function(e)
                {
                    shiny::showNotification("Error calculating MZ from neutral mass and adduct", type = "error")
                    NA
                })
                if (!is.na(mz))
                {
                    mzmin <- mz - input$mzWindowNM
                    mzmax <- mz + input$mzWindowNM
                }
            }
            else if (input$mzMode == "Formula +/- MZ Window")
            {
                neutralMass <- tryCatch(getFormulaMass(input$formula), error = function(e)
                {
                    shiny::showNotification("Error parsing formula", type = "error")
                    NA
                })
                if (!is.na(neutralMass))
                {
                    mz <- tryCatch(calculateMasses(neutralMass, as.adduct(input$adductF), "mz"), error = function(e)
                    {
                        shiny::showNotification("Error calculating MZ from formula and adduct", type = "error")
                        NA
                    })
                    if (!is.na(mz))
                    {
                        mzmin <- mz - input$mzWindowF
                        mzmax <- mz + input$mzWindowF
                    }
                }
            }
            if (is.null(mzmin) || is.null(mzmax))
            {
                mzmin <- 100
                mzmax <- 200
            }
            eimInfo <- data.table::data.table(
                mzmin = mzmin,
                mzmax = mzmax,
                retmin = if (input$rtMode == "Disabled") 0 else if (input$rtMode == "Range") input$retmin else input$retvalue - input$retwindow,
                retmax = if (input$rtMode == "Disabled") 0 else if (input$rtMode == "Range") input$retmax else input$retvalue + input$retwindow,
                mobmin = if (input$mobMode == "Disabled") 0 else if (input$mobMode == "Range") input$mobmin else input$mobvalue - input$mobwindow,
                mobmax = if (input$mobMode == "Disabled") 0 else if (input$mobMode == "Range") input$mobmax else input$mobvalue + input$mobwindow
            )
            eimList <- list(eimInfo)
            names(eimList) <- input$analysis
            sgLength <- if (input$enableEIMSmoothing) input$sgLength else 0
            shiny::validate(shiny::need(!input$enableEIMSmoothing || (sgLength %% 2 == 1 && sgLength >= 3), "Smoothing length must be odd and >= 3"))
            patRoon:::doGetEIMs(analysisInfo, eimList, minIntensity = 25, sgOrder = 3, sgLength = sgLength)
        })
        
        output$eicPlot <- shiny::renderPlot(
        {
            req(eicData())
            eic <- eicData()[[1]][[1]]  # Assuming single analysis and single range
            if (is.null(eic) || nrow(eic) == 0)
            {
                plot.new()
                text(0.5, 0.5, "No EIC data", cex = 2)
                return()
            }
            plot(eic[, "time"], eic[, "intensity"], type = "l", xlab = "Retention Time", ylab = "Intensity",
                 main = "Extracted Ion Chromatogram")
        })
        
        output$eimPlot <- shiny::renderPlot(
        {
            req(eimData())
            eim <- eimData()[[1]][[1]]  # Assuming single analysis and single range
            if (is.null(eim) || nrow(eim) == 0)
            {
                plot.new()
                text(0.5, 0.5, "No EIM data", cex = 2)
                return()
            }
            plot(eim[, "mobility"], eim[, "intensity"], type = "l", xlab = "Mobility", ylab = "Intensity",
                 main = "Extracted Ion Mobilogram")
        })
    }
    
    return(list(ui = ui, server = server))
}
