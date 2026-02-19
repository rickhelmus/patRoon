#' Shiny UI for configuring param objects
#'
#' This module provides a user interface for viewing and editing the \code{data}
#' slot of \code{param} objects. It automatically generates appropriate input
#' widgets based on the parameter definitions.
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
        
        # Get definitions from current param object
        getDefinitions <- shiny::reactive({
            paramRV()@definitions
        })
        
        # Get data from current param object
        getData <- shiny::reactive({
            as.list(paramRV())
        })
        
        # Generate dynamic UI for each parameter
        output$paramInputs <- shiny::renderUI({
            definitions <- getDefinitions()
            data <- getData()
            paramNames <- names(definitions)
            
            if (length(paramNames) == 0)
                return(shiny::helpText("No parameters to display."))
            
            # Create a list to hold all input widgets
            widgets <- lapply(paramNames, function(nm) {
                def <- definitions[[nm]]
                val <- data[[nm]]
                id <- paste0("param_", nm)
                
                # Create widget based on definition
                widget <- createParamInput(id, nm, val, def)
                
                # Wrap in a div with some styling
                shiny::div(
                    style = "margin-bottom: 15px; padding-bottom: 10px; border-bottom: 1px solid #eee;",
                    widget
                )
            })
            
            do.call(shiny::tagList, widgets)
        })
        
        # Observer to update param object when inputs change
        shiny::observe({
            definitions <- getDefinitions()
            paramNames <- names(definitions)
            
            if (length(paramNames) == 0)
                return()
            
            # Collect values from inputs
            newData <- lapply(paramNames, function(nm) {
                def <- definitions[[nm]]
                id <- paste0("param_", nm)
                
                # Get value based on definition type
                if (shiny::isTruthy(input[[id]]) || def$type == "flag") {
                    return(getInputValue(id, def, input))
                }
                return(paramRV()@data[[nm]])
            })
            names(newData) <- paramNames
            
            # Update param object using $<- operator (triggers validation)
            currParam <- paramRV()
            for (nm in paramNames) {
                # Use $<- for proper validation
                currParam[[nm]] <- newData[[nm]]
            }
            paramRV(currParam)
        })
        
        # Display param summary
        output$paramSummary <- shiny::renderPrint({
            p <- paramRV()
            cat(sprintf("Class: %s\n", p@name))
            cat(sprintf("Base name: %s\n", p@baseName))
            cat(sprintf("Description: %s\n", p@description))
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
            shiny::req(input$loadParam)
            
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

#' Create appropriate input widget for a parameter based on its definition
#'
#' @param id Input ID
#' @param name Parameter name (label)
#' @param val Current value
#' @param def Parameter definition list
#' @return A Shiny input widget
#' @noRd
createParamInput <- function(id, name, val, def)
{
    # Create label with description
    label <- shiny::tagList(
        shiny::strong(name),
        if (!is.null(def$description) && nchar(def$description) > 0) {
            shiny::tags$small(
                style = "display: block; color: #666; font-weight: normal;",
                def$description
            )
        }
    )
    
    # Handle based on definition type
    if (def$type == "flag") {
        return(shiny::checkboxInput(id, label, value = isTRUE(val)))
    }
    
    if (def$type == "choice") {
        return(shiny::selectInput(
            id, label, 
            choices = def$choices, 
            selected = val
        ))
    }
    
    if (def$type == "count") {
        return(shiny::numericInput(
            id, label,
            value = val,
            min = if (isTRUE(def$positive)) 0 else NA,
            step = 1
        ))
    }
    
    if (def$type == "number") {
        # Determine bounds
        min_val <- if (!is.null(def$lower)) def$lower else if (isTRUE(def$positive)) 0 else NA
        max_val <- if (!is.null(def$upper)) def$upper else NA
        
        # Determine step based on value magnitude
        step <- if (!is.null(val) && abs(val) >= 100) 1 else if (!is.null(val) && abs(val) >= 1) 0.1 else 0.001
        
        return(shiny::numericInput(
            id, label,
            value = val,
            min = min_val,
            max = max_val,
            step = step
        ))
    }
    
    if (def$type == "list") {
        # For list types, use text area for JSON editing
        jsonVal <- tryCatch(
            jsonlite::toJSON(val, auto_unbox = TRUE, pretty = TRUE),
            error = function(e) if (is.null(val)) "null" else "{}"
        )
        return(shiny::tagList(
            label,
            shiny::textAreaInput(id, NULL, value = jsonVal, rows = 4, width = "100%")
        ))
    }
    
    # Fallback for unknown types: text input
    shiny::textInput(id, label, value = as.character(val))
}

#' Get input value based on definition type
#'
#' @param id Input ID
#' @param def Parameter definition list
#' @param input Shiny input object
#' @return The parsed value
#' @noRd
getInputValue <- function(id, def, input)
{
    rawVal <- input[[id]]
    
    if (def$type == "flag") {
        return(isTRUE(rawVal))
    }
    
    if (def$type == "choice") {
        return(rawVal)
    }
    
    if (def$type == "count") {
        return(as.integer(rawVal))
    }
    
    if (def$type == "number") {
        return(as.numeric(rawVal))
    }
    
    if (def$type == "list") {
        if (is.null(rawVal) || rawVal == "" || rawVal == "null") {
            return(NULL)
        }
        return(tryCatch(
            jsonlite::fromJSON(rawVal),
            error = function(e) NULL
        ))
    }
    
    # Fallback
    return(rawVal)
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