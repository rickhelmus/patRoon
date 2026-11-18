# Throwaway Shiny GUI to review advanced flags without editing package source.
# Run from the package root: Rscript tools/param-advanced-review.R

for (pkg in c("shiny", "jsonlite")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Install %s first: install.packages('%s')", pkg, pkg))
    }
}

files <- c(
    "R/param.R",
    sort(list.files("R", pattern = "-param\\.R$", full.names = TRUE))
)
files <- unique(normalizePath(files, winslash = "/", mustWork = TRUE))

extract_file <- function(file) {
    lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
    starts <- grep("paramConfigDefsFact\\s*\\(\\s*list\\s*\\(", lines)
    if (!length(starts)) {
        return(NULL)
    }
    result <- list()
    for (ln in starts) {
        rest <- paste(lines[ln:length(lines)], collapse = "\n")
        call_start <- regexpr(
            "paramConfigDefsFact\\s*\\(\\s*list\\s*\\(",
            rest,
            perl = TRUE
        )[[1]]
        rest <- substring(rest, call_start)
        chars <- strsplit(rest, "", fixed = TRUE)[[1]]
        depth <- 0L
        quote <- ""
        esc <- FALSE
        end <- NA_integer_
        for (i in seq_along(chars)) {
            ch <- chars[[i]]
            if (nzchar(quote)) {
                if (esc) {
                    esc <- FALSE
                } else if (ch == "\\") {
                    esc <- TRUE
                } else if (ch == quote) {
                    quote <- ""
                }
                next
            }
            if (ch %in% c("'", '"', "`")) {
                quote <- ch
                next
            }
            if (ch == "(") {
                depth <- depth + 1L
            }
            if (ch == ")") {
                depth <- depth - 1L
                if (depth == 0L) {
                    end <- i
                    break
                }
            }
        }
        if (is.na(end)) {
            next
        }
        call <- substring(rest, 1L, end)
        tail <- paste(lines[ln:min(length(lines), ln + 180L)], collapse = "\n")
        cm <- regexpr(
            'setClass\\("[^\"]+",\\s*contains\\s*=\\s*"param"',
            tail,
            perl = TRUE
        )
        cls <- if (cm[[1]] > 0) {
            sub('^setClass\\("([^\"]+)".*$', '\\1', regmatches(tail, cm))
        } else {
            basename(file)
        }
        places <- gregexpr(
            "[A-Za-z][A-Za-z0-9._]*\\s*=\\s*list\\s*\\(",
            call,
            perl = TRUE
        )[[1]]
        if (places[[1]] < 0) {
            next
        }
        for (pos in places) {
            prefix <- substring(call, 1L, pos - 1L)
            d <- 0L
            q <- ""
            e <- FALSE
            for (ch in strsplit(prefix, "", fixed = TRUE)[[1]]) {
                if (nzchar(q)) {
                    if (e) {
                        e <- FALSE
                    } else if (ch == "\\") {
                        e <- TRUE
                    } else if (ch == q) {
                        q <- ""
                    }
                    next
                }
                if (ch %in% c("'", '"', "`")) {
                    q <- ch
                } else if (ch == "(") {
                    d <- d + 1L
                } else if (ch == ")") {
                    d <- d - 1L
                }
            }
            if (d != 2L) {
                next
            }
            entry <- substring(call, pos)
            name <- sub("^([A-Za-z][A-Za-z0-9._]*).*", "\\1", entry)
            lm <- regexpr("list\\s*\\(", entry, perl = TRUE)
            body <- substring(entry, attr(lm, "match.length") + 1L)
            bc <- strsplit(body, "", fixed = TRUE)[[1]]
            d <- 1L
            q <- ""
            e <- FALSE
            close <- length(bc)
            for (i in seq_along(bc)) {
                ch <- bc[[i]]
                if (nzchar(q)) {
                    if (e) {
                        e <- FALSE
                    } else if (ch == "\\") {
                        e <- TRUE
                    } else if (ch == q) {
                        q <- ""
                    }
                    next
                }
                if (ch %in% c("'", '"', "`")) {
                    q <- ch
                } else if (ch == "(") {
                    d <- d + 1L
                } else if (ch == ")") {
                    d <- d - 1L
                    if (d == 0L) {
                        close <- i
                        break
                    }
                }
            }
            body <- paste(bc[seq_len(max(0L, close - 1L))], collapse = "")
            adv <- regmatches(
                body,
                regexpr(
                    "(?m)^\\s*advanced\\s*=\\s*(TRUE|FALSE)",
                    body,
                    perl = TRUE
                )
            )
            result[[length(result) + 1L]] <- data.frame(
                file = substring(
                    file,
                    nchar(normalizePath(".", winslash = "/")) + 2L
                ),
                class = cls,
                parameter = name,
                advanced = length(adv) > 0L && grepl("TRUE$", trimws(adv)),
                stringsAsFactors = FALSE
            )
        }
    }
    if (length(result)) do.call(rbind, result) else NULL
}

found <- Filter(Negate(is.null), lapply(files, extract_file))
if (!length(found)) {
    stop("No parameter definitions found in R/param.R or R/*-param.R")
}
params <- do.call(rbind, found)
params <- params[
    !duplicated(params[c("file", "class", "parameter")]),
    ,
    drop = FALSE
]
params$key <- paste(params$file, params$class, params$parameter, sep = "::")
params$input_id <- paste0("adv_", seq_len(nrow(params)))

ui <- shiny::fluidPage(
    shiny::titlePanel("Parameter advanced-flag review"),
    shiny::p(
        "Choose which definitions should be advanced. This tool does not change R source files."
    ),
    shiny::fluidRow(
        shiny::column(
            4,
            shiny::textInput("search", "Search class, parameter, or file")
        ),
        shiny::column(
            4,
            shiny::selectInput(
                "file",
                "Source file",
                c("All files", sort(unique(params$file)))
            )
        ),
        shiny::column(
            4,
            shiny::selectInput(
                "status",
                "Filter by selection",
                c("All", "Advanced", "Regular")
            )
        )
    ),
    shiny::fluidRow(
        shiny::column(
            6,
            shiny::actionButton("previous_file", "Previous source file")
        ),
        shiny::column(6, shiny::actionButton("next_file", "Next source file"))
    ),
    shiny::fluidRow(
        shiny::column(
            3,
            shiny::actionButton("mark_advanced", "Mark visible advanced")
        ),
        shiny::column(
            3,
            shiny::actionButton("mark_regular", "Mark visible regular")
        ),
        shiny::column(
            3,
            shiny::textInput(
                "path",
                "Save JSON as",
                "tools/param-advanced-decisions.json"
            )
        ),
        shiny::column(
            3,
            shiny::textInput(
                "load_path",
                "Load JSON from",
                "tools/param-advanced-decisions.json"
            )
        )
    ),
    shiny::fluidRow(
        shiny::column(
            4,
            shiny::actionButton("save", "Save JSON", class = "btn-primary")
        ),
        shiny::column(4, shiny::actionButton("load", "Load JSON")),
        shiny::column(4, shiny::downloadButton("download", "Download JSON"))
    ),
    shiny::textOutput("summary"),
    shiny::uiOutput("table"),
    shiny::verbatimTextOutput("message")
)

server <- function(input, output, session) {
    values <- shiny::reactiveVal(setNames(params$advanced, params$key))
    file_choices <- sort(unique(params$file))
    visible <- shiny::reactive({
        x <- params
        if (input$file != "All files") {
            x <- x[x$file == input$file, , drop = FALSE]
        }
        state <- values()[x$key]
        if (input$status == "Advanced") {
            x <- x[which(state %in% TRUE), , drop = FALSE]
        }
        if (input$status == "Regular") {
            x <- x[which(state %in% FALSE), , drop = FALSE]
        }
        q <- tolower(trimws(input$search))
        if (nzchar(q)) {
            x <- x[
                grepl(
                    q,
                    tolower(paste(x$file, x$class, x$parameter)),
                    fixed = TRUE
                ),
                ,
                drop = FALSE
            ]
        }
        x
    })
    navigate <- function(direction) {
        index <- match(input$file, file_choices)
        if (is.na(index)) {
            index <- if (direction > 0) 0L else length(file_choices) + 1L
        }
        index <- ((index - 1L + direction) %% length(file_choices)) + 1L
        shiny::updateSelectInput(
            session,
            "file",
            selected = file_choices[[index]]
        )
    }
    shiny::observeEvent(input$previous_file, navigate(-1L))
    shiny::observeEvent(input$next_file, navigate(1L))
    output$summary <- shiny::renderText(sprintf(
        "Visible %d / %d | Advanced %d | Regular %d",
        nrow(visible()),
        nrow(params),
        sum(values()),
        length(values()) - sum(values())
    ))
    output$table <- shiny::renderUI({
        x <- visible()
        if (!nrow(x)) {
            return(shiny::p("No matching parameters."))
        }
        shiny::tagList(
            shiny::fluidRow(
                shiny::column(3, shiny::strong("Source file")),
                shiny::column(3, shiny::strong("Class")),
                shiny::column(3, shiny::strong("Parameter")),
                shiny::column(3, shiny::strong("Advanced"))
            ),
            lapply(seq_len(nrow(x)), function(i) {
                key <- x$key[[i]]
                shiny::fluidRow(
                    shiny::column(3, x$file[[i]]),
                    shiny::column(3, x$class[[i]]),
                    shiny::column(3, x$parameter[[i]]),
                    shiny::column(
                        3,
                        shiny::checkboxInput(
                            params$input_id[match(key, params$key)],
                            NULL,
                            value = isTRUE(values()[[key]])
                        )
                    )
                )
            })
        )
    })
    for (i in seq_len(nrow(params))) {
        local({
            row <- params[i, , drop = FALSE]
            shiny::observeEvent(input[[row$input_id]], {
                state <- isTRUE(input[[row$input_id]])
                v <- values()
                if (!identical(v[[row$key]], state)) {
                    v[[row$key]] <- state
                    values(v)
                    message(sprintf(
                        "Checkbox update: %s (%s) = %s",
                        row$key,
                        row$input_id,
                        state
                    ))
                }
            }, ignoreInit = TRUE)
        })
    }
    shiny::observeEvent(input$mark_advanced, {
        v <- values()
        v[visible()$key] <- TRUE
        values(v)
    })
    shiny::observeEvent(input$mark_regular, {
        v <- values()
        v[visible()$key] <- FALSE
        values(v)
    })
    json <- shiny::reactive({
        v <- values()
        payload <- list(
            schema_version = 1L,
            description = "Desired advanced status for parameter definitions; decisions only, not a source patch.",
            decisions = lapply(seq_len(nrow(params)), function(i) {
                list(
                    file = params$file[[i]],
                    class = params$class[[i]],
                    parameter = params$parameter[[i]],
                    advanced = unname(v[[params$key[[i]]]])
                )
            })
        )
        jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = TRUE)
    })
    output$download <- shiny::downloadHandler(
        "param-advanced-decisions.json",
        content = function(file) writeLines(json(), file, useBytes = TRUE)
    )
    message_text <- shiny::reactiveVal("")
    output$message <- shiny::renderText(message_text())
    shiny::observeEvent(input$load, {
        target <- trimws(input$load_path)
        if (!nzchar(target)) {
            message_text("Enter a JSON input path.")
            return()
        }
        result <- tryCatch(
            {
                payload <- jsonlite::fromJSON(target, simplifyVector = FALSE)
                if (
                    !is.list(payload) ||
                        !is.list(payload$decisions) ||
                        !length(payload$decisions)
                ) {
                    stop("JSON must contain a non-empty 'decisions' array.")
                }
                decisions <- payload$decisions
                required <- c("file", "class", "parameter", "advanced")
                valid <- vapply(
                    decisions,
                    function(item) {
                        is.list(item) &&
                            all(required %in% names(item)) &&
                            all(vapply(
                                item[required[1:3]],
                                function(x) {
                                    is.character(x) &&
                                        length(x) == 1L &&
                                        !is.na(x)
                                },
                                logical(1)
                            )) &&
                            is.logical(item$advanced) &&
                            length(item$advanced) == 1L &&
                            !is.na(item$advanced)
                    },
                    logical(1)
                )
                if (!all(valid)) {
                    stop(
                        "Each decision must contain file, class, parameter, and a logical advanced value."
                    )
                }
                keys <- vapply(
                    decisions,
                    function(item) {
                        paste(item$file, item$class, item$parameter, sep = "::")
                    },
                    character(1)
                )
                matched <- keys %in% names(values())
                if (!any(matched)) {
                    stop(
                        "No decisions in this file match the current parameter definitions."
                    )
                }
                v <- values()
                for (i in which(matched)) {
                    v[[keys[[i]]]] <- decisions[[i]]$advanced
                }
                values(v)
                paste(
                    "Loaded",
                    sum(matched),
                    "matching decisions from",
                    normalizePath(target, winslash = "/")
                )
            },
            error = function(e) {
                paste("Could not load JSON:", conditionMessage(e))
            }
        )
        message_text(result)
    })
    shiny::observeEvent(input$save, {
        target <- trimws(input$path)
        if (!nzchar(target)) {
            message_text("Enter an output path.")
            warning("JSON save skipped: output path is empty.", call. = FALSE)
            return()
        }
        payload <- json()
        destination <- normalizePath(target, winslash = "/", mustWork = FALSE)
        message(sprintf(
            "JSON save requested: path=%s; advanced=%d/%d",
            destination,
            sum(values()),
            length(values())
        ))
        result <- tryCatch({
            dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
            writeLines(payload, target, useBytes = TRUE)
            info <- file.info(target)
            if (is.na(info$size) || info$size == 0) {
                stop("Output file is missing or empty after write.")
            }
            message(sprintf(
                "JSON save succeeded: path=%s; bytes=%s",
                normalizePath(target, winslash = "/", mustWork = TRUE),
                info$size
            ))
            paste(
                "Saved",
                normalizePath(target, winslash = "/", mustWork = TRUE),
                sprintf("(%s bytes)", info$size)
            )
        }, error = function(e) {
            warning(sprintf(
                "JSON save failed for %s: %s",
                destination,
                conditionMessage(e)
            ), call. = FALSE)
            paste("Could not save JSON:", conditionMessage(e))
        })
        message_text(result)
    })
}

shiny::runApp(shiny::shinyApp(ui, server), launch.browser = TRUE)
