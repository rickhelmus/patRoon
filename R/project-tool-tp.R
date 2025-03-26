getTPGenInputs <- function(isLib)
{
    ret <- c("Suspect list" = "suspects",
             "Screening results" = "screening")
    if (isLib)
        ret <- c(ret, "All library parents" = "all")
    return(ret)
}

newProjectTPUI <- function(id)
{
    ns <- NS(id)
    
    miniUI::miniContentPanel(
        fillCol(
            flex = NA,
            checkboxInput(ns("doTPs"), "Perform transformation product screening"),
            conditionalPanel(
                condition = "input.doTPs",
                ns = ns,
                selectInput(ns("TPGen"), "TP algorithm", c("BioTransformer", "CTS", "Library", "Logic")),
                conditionalPanel(
                    condition = "input.TPGen != \"Logic\"",
                    ns = ns,
                    selectInput(ns("TPGenInput"), "Parent input", getTPGenInputs(FALSE)),
                    conditionalPanel(
                        condition = "input.TPGenInput == \"suspects\"",
                        ns = ns,
                        fileSelect(ns("TPSuspectList"), ns("TPSuspButton"), "Parent suspect list",
                                   placeholder = "Please specify parent suspect list")
                    ),
                    checkboxInput(ns("TPDoMFDB"), "Generate TP MetFrag database", TRUE)
                )
            )
        )
    )
}

newProjectTPServer <- function(id, hasSuspects, settings)
{
    moduleServer(id, function(input, output, session)
    {
        observeEvent(settings(), {
            updateCheckboxInput(session, "doTPs", value = settings()$doTPs)
            updateSelectInput(session, "TPGen", selected = settings()$TPGen)
            updateSelectInput(session, "TPGenInput", selected = settings()$TPGenInput)
            updateTextInput(session, "TPSuspectList", value = settings()$TPSuspectList)
            updateCheckboxInput(session, "TPDoMFDB", value = settings()$TPDoMFDB)
        })
        
        observeEvent(input$TPGen, {
            updateSelectInput(inputId = "TPGenInput", choices = getTPGenInputs(input$TPGen == "Library"))
        })
        
        observeEvent(input$TPGenInput, {
            if (input$TPGenInput == "screening" && !hasSuspects())
            {
                rstudioapi::showDialog("Enable suspect screening",
                                       "This requires a workflow with suspect screening. Please configure in the Features tab.", "")
                updateSelectInput(inputId = "TPGenInput", selected = "suspects")
            }
        })
        
        observeEvent(input$TPSuspButton, selectSuspList(session, "TPSuspectList"))
        
        list(
            valid = reactive({
                if (input$doTPs && input$TPGenInput == "suspects" && !nzchar(input$TPSuspectList))
                    list(title = "No suspect list", msg = "Please select a suspect list!")
                else
                    TRUE
            }),
            settings = reactive(list(
                doTPs = input$doTPs,
                TPGen = input$TPGen,
                TPGenInput = input$TPGenInput,
                TPSuspectList = input$TPSuspectList,
                TPDoMFDB = input$TPDoMFDB
            ))
        )
    })
}
