getWorkPath <- function(file = "", ...) if (nzchar(file)) file.path("test_temp", file, ...) else "test_temp"
getTestDataPath <- function() "test_data"
testFile <- function(f, ..., text = FALSE) file.path(getTestDataPath(), paste0(f, ..., if (!text) ".Rds" else ".txt", collapse = ""))
getTestAnaInfo <- function(path = patRoonData::exampleDataPath()) generateAnalysisInfo(path,
                                                                                       groups = c(rep("solvent", 3), rep("standard", 3)),
                                                                                       refs = "solvent")
getTestFGroups <- function(anaInfo = getTestAnaInfo()) groupFeatures(findFeatures(anaInfo, "openms", logPath = NULL), "openms")
getEmptyTestFGroups <- function() getTestFGroups()[, "none"]
getEmptyPLists <- function() MSPeakLists()

makeMZXMLs <- function(anaInfo)
{
    exDataFiles <- list.files(patRoonData::exampleDataPath(), "\\.mzML$", full.names = TRUE)
    convertMSFiles(anaInfo, "mzML", "mzXML", getWorkPath())
    return(getTestAnaInfo(getWorkPath()))
}

expect_file <- function(object, file, removeIfExists = TRUE)
{
    if (removeIfExists && file.exists(file))
        file.remove(file)

    act <- quasi_label(rlang::enquo(object))
    expect(file.exists(file), sprintf("failed to generate %s", file))
    invisible(act$val)
}

expect_range <- function(object, r)
{
    act <- quasi_label(rlang::enquo(object))
    act$r <- range(act$val)
    expect(act$r[1] >= r[1] && act$r[2] <= r[2],
           sprintf("range of %s is %.1f - %.1f which is outside %.1f - %.1f",
                   act$lab, act$r[1], act$r[2], r[1], r[2]))
    invisible(act$val)
}

expect_known_show <- function(object, file)
{
    act <- quasi_label(rlang::enquo(object))

    text <- capture_output_lines(show(act$val))

    # remove last line with object size as it may vary even if object remain the same
    text <- text[seq_len(length(text)-1)]

    text <- paste0(text, collapse = "")

    # based on simplified expect_known_output
    if (!file.exists(file))
    {
        warning("Creating reference output", call. = FALSE)
        cat(text, file = file)
        succeed()
    }
    else
    {
        ref <- patRoon:::readAllFile(file)
        cat(text, file = file)

        cmp <- compare(text, enc2native(ref))
        expect(cmp$equal, sprintf("show reference of %s has changed\n%s", act$lab, cmp$message))
    }

    invisible(act$val)
}

expect_plot <- function(object)
{
    tf <- tempfile()
    withr::with_png(tf, act <- quasi_label(rlang::enquo(object)))
    expect(file.exists(tf), "failed to generate plot")
    invisible(act$val)
}

makeReportMD <- function(fGroups, ...) reportMD(fGroups, getWorkPath(), openReport = FALSE, ...)
expect_reportMD <- function(object)
{
    act <- quasi_label(rlang::enquo(object))
    expect(file.exists(getWorkPath("report.html")), "failed to generate report")
    invisible(act$val)
}

expect_doppel <- function(...) vdiffr::expect_doppelganger(..., path = getOption("patRoon.path.vdiffr"))

# HACK: workaround for non imported checkmate namespace
makeExpectation <- checkmate::makeExpectation
vname <- checkmate::vname

expect_csv_file <- checkmate::makeExpectationFunction(patRoon:::checkCSVFile)

# call dollar operator with string
callDollar <- function(x, name) eval(substitute(x$NAME_ARG, list(NAME_ARG = name)))
