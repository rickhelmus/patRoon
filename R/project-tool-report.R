newProjectReportUI <- function(id)
{
    ns <- NS(id)
    
    miniUI::miniContentPanel(
        fillRow(
            checkboxGroupInput(ns("reportGen"), "Report generation", c("HTML reports" = "HTML",
                                                                       "Legacy interface" = "legacy"),
                               "HTML", width = "100%"),
            conditionalPanel(
                condition = "input.reportGen.includes('legacy')",
                ns = ns,
                checkboxGroupInput(ns("reportLegacy"), "Legacy report formats", c("CSV", "PDF"),
                                   "CSV", width = "100%")
            )
        )
    )
}

newProjectReportServer <- function(id, settings)
{
    moduleServer(id, function(input, output, session)
    {
        list(
            reportGen = reactive(input$reportGen),
            reportLegacy = reactive(input$reportLegacy)
        )
    })
}
