# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

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
        observeEvent(settings(), {
            updateCheckboxGroupInput(session, "reportGen", selected = settings()$reportGen)
            updateCheckboxGroupInput(session, "reportLegacy", selected = settings()$reportLegacy)
        })
        
        list(
            valid = reactive({
                if ("legacy" %in% input$reportGen && length(input$reportLegacy) == 0)
                    list(title = "Invalid report settings", msg = "Please select at least one report format!")
                else
                    TRUE
            }),
            settings = reactive(list(
                reportGen = input$reportGen,
                reportLegacy = input$reportLegacy
            ))
        )
    })
}

defaultReportSettings <- function()
{
    return(list(reportGen = "HTML", reportLegacy = "CSV"))
}

upgradeReportSettings <- function(settings)
{
    # NOTE: this updates from first file version
    ret <- modifyList(defaultReportSettings(), settings[c("reportGen", "reportLegacy")])
    ret$reportLegacy <- setdiff(ret$reportLegacy, "HTML")
    return(ret)
}
