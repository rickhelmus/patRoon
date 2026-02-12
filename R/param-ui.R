#' Shiny UI for configuring param objects
#'
#' This module provides a user interface for viewing and editing the \code{data}
#' slot of \code{param} objects. It automatically generates appropriate input
#' widgets based on the type of each parameter.
#'
#' @param obj A \code{param} object to edit. If \code{NULL}, a default
#'   \code{FeatureOpenMSParam} object will be created.
#'
#' @return A list containing \code{ui} and \code{server} components that can be
#'   passed to \code{shiny::shinyApp()}.
#'
#' @examples
#' \dontrun{
#' library(shiny)
#' paramObj <- FeatureOpenMSParam()
#' gui <- createParamUI(paramObj)
#' shinyApp(gui$ui, gui$server)
#' }
#'
#' @export
createParamUI <- function(obj = NULL)
{
    if (is.null(obj))
        obj <- FeatureOpenMSParam()
    
    # Store initial object for reset functionality
    initObj <- obj
    
    ui <- shiny::fluidPage(
        shiny::titlePanel("Parameter Configuration"),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                width = 3,
                shiny::h4("Actions"),
                shiny::fileInput("loadParam", "Load from R Data", accept = c(".rda", ".RData")),
                shiny::downloadButton("saveParam", "Save to R Data"),
                shiny::hr(),
                shiny::selectInput("exportFormat", "Export Format", 
                                   choices = c("R" = "R", "JSON" = "json", "YAML" = "yaml")),
                shiny::downloadButton("exportParam", "Export"),
                shiny::hr(),
                shiny::actionButton("resetParam", "Reset to Defaults", class = "btn-warning"),
                shiny::hr(),
                shiny::verbatimTextOutput("paramSummary"),
                shiny::hr(),
                shiny::uiOutput("validationStatus")
            ),
            shiny::mainPanel(
                width = 9,
                shiny::div(
                    style = "overflow-y: auto; max-height: 85vh; padding-right: 15px;",
                    shiny::h3("Parameter Values"),
                    shiny::uiOutput("paramInputs"),
                    shiny::hr(),
                    shiny::h4("Raw Data Preview"),
                    shiny::verbatimTextOutput("dataPreview")
                )
            )
        )
    )
    
    server <- function(input, output, session)
    {
        # Reactive value to store the current param object
        paramRV <- shiny::reactiveVal(obj)
        
        # Store initial object for reset
        initParam <- initObj
        
        # Get data from current param object
        getData <- shiny::reactive({
            as.list(paramRV())
        })
        
        # Generate dynamic UI for each parameter
        output$paramInputs <- shiny::renderUI({
            data <- getData()
            paramNames <- names(data)
            
            if (length(paramNames) == 0)
                return(shiny::helpText("No parameters to display."))
            
            # Create a list to hold all input widgets
            widgets <- lapply(paramNames, function(nm) {
                val <- data[[nm]]
                id <- paste0("param_", nm)
                
                # Determine widget type based on value
                widget <- createParamInput(id, nm, val)
                
                # Wrap in a div with some styling
                shiny::div(
                    style = "margin-bottom: 15px;",
                    widget
                )
            })
            
            do.call(shiny::tagList, widgets)
        })
        
        # Observer to update param object when inputs change
        shiny::observe({
            data <- getData()
            paramNames <- names(data)
            
            if (length(paramNames) == 0)
                return()
            
            # Collect values from inputs
            newData <- lapply(paramNames, function(nm) {
                id <- paste0("param_", nm)
                val <- data[[nm]]
                
                # Get value based on type
                if (shiny::isTruthy(input[[id]])) {
                    if (is.logical(val)) {
                        return(input[[id]])
                    } else if (is.numeric(val)) {
                        return(input[[id]])
                    } else if (is.character(val)) {
                        return(input[[id]])
                    } else {
                        # Try to parse as JSON for complex types
                        tryCatch(
                            jsonlite::fromJSON(input[[id]]),
                            error = function(e) val
                        )
                    }
                }
                return(val)
            })
            names(newData) <- paramNames
            
            # Update param object by modifying the data slot directly
            currParam <- paramRV()
            currData <- currParam@data
            for (nm in paramNames) {
                currData[[nm]] <- newData[[nm]]
            }
            currParam@data <- currData
            paramRV(currParam)
        })
        
        # Display param summary
        output$paramSummary <- shiny::renderPrint({
            p <- paramRV()
            cat(sprintf("Class: %s\n", p@name))
            cat(sprintf("Base name: %s\n", p@baseName))
            cat(sprintf("Version: %s\n", p@version))
            cat(sprintf("Date: %s\n", format(p@date)))
        })
        
        # Display raw data preview
        output$dataPreview <- shiny::renderPrint({
            utils::str(as.list(paramRV()))
        })
        
        # Validation status
        output$validationStatus <- shiny::renderUI({
            p <- paramRV()
            valid <- tryCatch(
                {
                    validObject(p)
                    TRUE
                },
                error = function(e) e$message
            )
            
            if (isTRUE(valid)) {
                shiny::div(
                    style = "color: green; font-weight: bold;",
                    shiny::icon("check-circle"),
                    "Valid"
                )
            } else {
                shiny::div(
                    style = "color: red;",
                    shiny::icon("exclamation-triangle"),
                    "Invalid: ", valid
                )
            }
        })
        
        # Load param from file
        shiny::observeEvent(input$loadParam, {
            req(input$loadParam)
            
            tryCatch({
                env <- new.env()
                load(input$loadParam$datapath, envir = env)
                
                # Find the param object in the loaded environment
                objNames <- ls(env)
                paramObj <- NULL
                for (nm in objNames) {
                    obj <- get(nm, envir = env)
                    if (is(obj, "param")) {
                        paramObj <- obj
                        break
                    }
                }
                
                if (!is.null(paramObj)) {
                    paramRV(paramObj)
                    shiny::showNotification("Parameters loaded successfully", type = "message")
                } else {
                    shiny::showNotification("No param object found in file", type = "error")
                }
            }, error = function(e) {
                shiny::showNotification(paste("Error loading file:", e$message), type = "error")
            })
        })
        
        # Save param to file
        output$saveParam <- shiny::downloadHandler(
            filename = function() {
                paste0("param_", paramRV()@baseName, ".rda")
            },
            content = function(file) {
                paramObj <- paramRV()
                save(paramObj, file = file)
            }
        )
        
        # Export param in various formats
        output$exportParam <- shiny::downloadHandler(
            filename = function() {
                ext <- switch(input$exportFormat,
                              "R" = ".R",
                              "json" = ".json",
                              "yaml" = ".yaml",
                              ".txt")
                paste0("param_", paramRV()@baseName, ext)
            },
            content = function(file) {
                export(paramRV(), input$exportFormat, file)
            }
        )
        
        # Reset to defaults
        shiny::observeEvent(input$resetParam, {
            shiny::showModal(shiny::modalDialog(
                title = "Confirm Reset",
                "Are you sure you want to reset all parameters to their default values?",
                footer = shiny::tagList(
                    shiny::modalButton("Cancel"),
                    shiny::actionButton("confirmReset", "Reset", class = "btn-warning")
                )
            ))
        })
        
        shiny::observeEvent(input$confirmReset, {
            shiny::removeModal()
            # Create new default object based on class
            currClass <- class(initParam)
            newObj <- do.call(currClass[[1]], list())
            paramRV(newObj)
            shiny::showNotification("Parameters reset to defaults", type = "message")
        })
    }
    
    return(list(ui = ui, server = server))
}

#' Create appropriate input widget for a parameter value
#'
#' @param id Input ID
#' @param name Parameter name (label)
#' @param val Current value
#' @return A Shiny input widget
#' @noRd
createParamInput <- function(id, name, val)
{
    # NULL values -> text area for JSON
    if (is.null(val)) {
        return(shiny::textAreaInput(id, name, value = "null", rows = 2))
    }
    
    # Logical/boolean -> checkbox
    if (is.logical(val) && length(val) == 1) {
        return(shiny::checkboxInput(id, name, value = val))
    }
    
    # Integer (positive count) -> numeric input with integer constraint
    if (is.numeric(val) && length(val) == 1 && 
        name %in% c("traceTermOutliers", "minTraceLength", "maxTraceLength")) {
        return(shiny::numericInput(id, name, value = val, step = 1))
    }
    
    # Numeric scalar -> numeric input
    if (is.numeric(val) && length(val) == 1) {
        # Determine appropriate step based on value magnitude
        step <- if (abs(val) >= 100) 1 else if (abs(val) >= 1) 0.1 else 0.001
        return(shiny::numericInput(id, name, value = val, step = step))
    }
    
    # Numeric vector length 2 -> two numeric inputs
    if (is.numeric(val) && length(val) == 2) {
        return(shiny::fluidRow(
            shiny::column(6, shiny::numericInput(paste0(id, "_min"), 
paste0(name, " (min)"), value = val[1])),
            shiny::column(6, shiny::numericInput(paste0(id, "_max"), 
paste0(name, " (max)"), value = val[2]))
        ))
    }
    
    # Character choices -> select input
    if (is.character(val) && length(val) == 1) {
        choices <- inferChoices(name, val)
        if (!is.null(choices)) {
            return(shiny::selectInput(id, name, choices = choices, selected = val))
        }
        # Free text fallback for single character
        return(shiny::textInput(id, name, value = val))
    }
    
    # Character vector -> select input with multiple selection
    if (is.character(val) && length(val) > 0) {
        return(shiny::selectInput(id, name, choices = val, selected = val, multiple = TRUE))
    }
    
    # List -> text area for JSON editing
    if (is.list(val)) {
        jsonVal <- tryCatch(
            jsonlite::toJSON(val, auto_unbox = TRUE, pretty = TRUE),
            error = function(e) "{}"
        )
        return(shiny::tagList(
            shiny::strong(name),
            shiny::textAreaInput(id, NULL, value = jsonVal, rows = 4, width = "100%")
        ))
    }
    
    # Fallback: text area for any other type
    shiny::textAreaInput(id, name, value = paste(as.character(val), collapse = ", "), rows = 2)
}

#' Infer choice options for known parameters
#'
#' @param name Parameter name
#' @param val Current value
#' @return Character vector of choices or NULL
#' @noRd
inferChoices <- function(name, val)
{
    choices <- switch(name,
                      "traceTermCriterion" = c("sample_rate", "outliers"),
                      "widthFiltering" = c("fixed", "adaptive"),
                      "isotopeFilteringModel" = c("metabolites (5% RMS)", "metabolites (10% RMS)", 
                                                  "lipids (5% RMS)", "lipids (10% RMS)"),
                      NULL
    )
    return(choices)
}

#' Launch the parameter configuration UI
#'
#' @param obj A \code{param} object to edit. If \code{NULL}, a default
#'   \code{FeatureOpenMSParam} object will be created.
#'
#' @return This function launches a Shiny app and does not return.
#'
#' @examples
#' \dontrun{
#' # Launch with default parameters
#' configureParams()
#'
#' # Launch with existing parameters
#' params <- FeatureOpenMSParam()
#' configureParams(params)
#' }
#'
#' @export
configureParams <- function(obj = NULL)
{
    gui <- createParamUI(obj)
    shiny::shinyApp(gui$ui, gui$server)
}
